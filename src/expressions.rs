//! Polars expression entry points.
//!
//! Arguments arrive already validated by the Python layer (`python/polars_seq/__init__.py`),
//! so a bad `table` or a two-character `stop_symbol` raises a clean `ValueError` at expression
//! construction time rather than a `ComputeError` from deep inside the query engine.

use polars::prelude::*;
use pyo3_polars::derive::polars_expr;
use serde::Deserialize;

use crate::codon_tables::table_by_id;
use crate::translate::{reverse_complement, translate};

#[derive(Deserialize)]
pub struct TranslateKwargs {
    table: u8,
    stop_symbol: char,
    to_stop: bool,
    cds: bool,
    gap: Option<char>,
    /// If true, a sequence that BioPython would reject becomes `null` instead of aborting
    /// the whole query. See README "Differences from BioPython".
    null_on_error: bool,
}

#[polars_expr(output_type=String)]
fn translate_expr(inputs: &[Series], kwargs: TranslateKwargs) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;

    let tbl = table_by_id(kwargs.table)
        .ok_or_else(|| polars_err!(ComputeError: "unknown codon table id {}", kwargs.table))?;

    let gap: Option<u8> = kwargs.gap.map(|c| c as u8);

    let mut builder = StringChunkedBuilder::new(ca.name().clone(), ca.len());
    // One buffer for the whole column: the per-row output is built in place and copied into
    // the Arrow buffer, so translation itself allocates nothing.
    let mut buf = String::new();

    for (i, opt) in ca.iter().enumerate() {
        let Some(s) = opt else {
            builder.append_null();
            continue;
        };
        buf.clear();
        match translate(
            s.as_bytes(),
            tbl,
            kwargs.stop_symbol,
            kwargs.to_stop,
            kwargs.cds,
            gap,
            &mut buf,
        ) {
            Ok(()) => builder.append_value(&buf),
            Err(e) if kwargs.null_on_error => {
                let _ = e;
                builder.append_null();
            }
            Err(e) => {
                // Quote the row so a failure in a million-row frame is actionable.
                return Err(polars_err!(ComputeError: "{e} (row {i}, sequence '{}')", truncate(s)));
            }
        }
    }

    Ok(builder.finish().into_series())
}

#[polars_expr(output_type=String)]
fn reverse_complement_expr(inputs: &[Series]) -> PolarsResult<Series> {
    let ca = inputs[0].str()?;
    let mut builder = StringChunkedBuilder::new(ca.name().clone(), ca.len());
    let mut buf: Vec<u8> = Vec::new();

    for opt in ca.iter() {
        let Some(s) = opt else {
            builder.append_null();
            continue;
        };
        reverse_complement(s.as_bytes(), &mut buf);
        // The complement table only ever maps ASCII to ASCII and passes other bytes through
        // untouched, so a valid UTF-8 input stays valid UTF-8 -- unless it was multi-byte, in
        // which case reversing the bytes would corrupt it. Guard that case.
        match std::str::from_utf8(&buf) {
            Ok(rc) => builder.append_value(rc),
            Err(_) => {
                return Err(polars_err!(
                    ComputeError: "reverse_complement requires ASCII sequences, got '{}'", truncate(s)
                ))
            }
        }
    }

    Ok(builder.finish().into_series())
}

/// Keep error messages readable when the offending sequence is long.
fn truncate(s: &str) -> String {
    const MAX: usize = 40;
    if s.len() <= MAX {
        return s.to_string();
    }
    let cut = s
        .char_indices()
        .map(|(i, _)| i)
        .take_while(|&i| i <= MAX)
        .last()
        .unwrap_or(0);
    format!("{}...", &s[..cut])
}
