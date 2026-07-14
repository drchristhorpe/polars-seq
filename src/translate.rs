//! The translation kernel.
//!
//! A faithful port of BioPython's `Bio.Seq._translate_str`. All of the genuinely subtle
//! logic -- resolving ambiguous IUPAC codons -- has already been pre-computed into the
//! static tables in `codon_tables.rs`, so what remains here is the control flow: framing,
//! stops, gaps, and the `cds` validation rules. The rules and their ordering are load
//! bearing; see PLAN.md section 4.

use crate::codon_tables::{CodonTable, ALPHABET, INVALID, STOP};
use std::fmt;

/// Maps an input byte to its position in `ALPHABET` (0..=16), or `NOT_A_LETTER`.
///
/// Upper- and lower-case fold to the same code, which is how BioPython's leading `.upper()`
/// is implemented here without allocating an uppercased copy of every sequence.
const NOT_A_LETTER: u8 = 0xFF;

static LETTER: [u8; 256] = build_letter();

const fn build_letter() -> [u8; 256] {
    let mut t = [NOT_A_LETTER; 256];
    let mut i = 0;
    while i < 17 {
        let c = ALPHABET[i];
        t[c as usize] = i as u8;
        t[(c + 32) as usize] = i as u8; // lower-case
        i += 1;
    }
    t
}

/// Codon -> table index (base 17), or `None` if a letter is not a nucleotide code at all.
///
/// Note `None` means "not even a letter we could look up" (e.g. '-' or '?'), which is a
/// different thing from a codon that *is* looked up and turns out to be `INVALID` (e.g. XXX).
#[inline(always)]
fn codon_index(c: &[u8]) -> Option<usize> {
    let a = LETTER[c[0] as usize];
    let b = LETTER[c[1] as usize];
    let d = LETTER[c[2] as usize];
    if a > 16 || b > 16 || d > 16 {
        return None;
    }
    Some((a as usize) * 289 + (b as usize) * 17 + (d as usize))
}

#[derive(Debug)]
pub enum TranslateError {
    InvalidCodon(String),
    NotStartCodon(String),
    NotMultipleOfThree(usize),
    NotStopCodon(String),
    ExtraStop,
    NonAscii,
}

// Message text mirrors BioPython's `CodonTable.TranslationError` strings, so anyone who has
// seen these errors from BioPython recognises them immediately.
impl fmt::Display for TranslateError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidCodon(c) => write!(f, "Codon '{c}' is invalid"),
            Self::NotStartCodon(c) => write!(f, "First codon '{c}' is not a start codon"),
            Self::NotMultipleOfThree(n) => write!(f, "Sequence length {n} is not a multiple of three"),
            Self::NotStopCodon(c) => write!(f, "Final codon '{c}' is not a stop codon"),
            Self::ExtraStop => write!(f, "Extra in frame stop codon found."),
            Self::NonAscii => write!(f, "Sequence contains non-ASCII characters"),
        }
    }
}

fn upper(bytes: &[u8]) -> String {
    bytes.iter().map(|b| b.to_ascii_uppercase() as char).collect()
}

/// Translate one sequence, appending to `out` (reused across rows to avoid re-allocating).
pub fn translate(
    seq: &[u8],
    tbl: &CodonTable,
    stop_symbol: char,
    to_stop: bool,
    cds: bool,
    gap: Option<u8>,
    out: &mut String,
) -> Result<(), TranslateError> {
    // We frame by bytes, which is only equivalent to BioPython's per-character framing for
    // ASCII input. Nucleotide data always is; anything else is rejected rather than silently
    // mis-framed.
    if !seq.is_ascii() {
        return Err(TranslateError::NonAscii);
    }

    let n = seq.len();

    // `body` is the region that goes through the codon loop. Under `cds` the start and stop
    // codons are handled separately and excluded from it.
    let body: &[u8] = if cds {
        // Order matters: BioPython checks start, then length, then final stop. A sequence can
        // fail more than one of these and the message you get depends on this ordering.
        if n < 3 {
            return Err(TranslateError::NotStartCodon(upper(seq)));
        }
        let first = &seq[..3];
        match codon_index(first) {
            Some(i) if tbl.is_start_codon(i) => {}
            _ => return Err(TranslateError::NotStartCodon(upper(first))),
        }
        if n % 3 != 0 {
            return Err(TranslateError::NotMultipleOfThree(n));
        }
        let last = &seq[n - 3..];
        match codon_index(last) {
            Some(i) if tbl.is_stop_codon(i) => {}
            _ => return Err(TranslateError::NotStopCodon(upper(last))),
        }
        // The start codon is reported as Met regardless of what it actually encodes, and the
        // terminal stop is dropped.
        out.push('M');
        &seq[3..n - 3]
    } else {
        // A trailing partial codon is dropped. BioPython warns and does the same; we cannot
        // warn per-row from a parallel kernel, so we are simply quiet about it.
        &seq[..n - n % 3]
    };

    for codon in body.chunks_exact(3) {
        match codon_index(codon) {
            Some(idx) => match tbl.code[idx] {
                STOP => {
                    if cds {
                        return Err(TranslateError::ExtraStop);
                    }
                    if to_stop {
                        break;
                    }
                    out.push(stop_symbol);
                }
                // An 'X'-bearing codon that does not resolve to a single amino acid. BioPython
                // rejects it, because 'X' is not in its set of valid input letters -- even
                // though 'X' *is* a nucleotide expansion key, which is why e.g. CTX -> 'L' is
                // accepted while XXX is not.
                INVALID => return Err(TranslateError::InvalidCodon(upper(codon))),
                // An amino acid, or 'X' where the codon was too ambiguous to pin down. A
                // dual-coding codon (tables 27/28/31) lands here as an amino acid even though
                // it is also in the stop set -- which is exactly why `code` and `is_stop` are
                // stored independently.
                code => out.push(code as char),
            },
            None => {
                // The codon contains something that is not a nucleotide letter at all. The gap
                // check comes after the table lookup, matching BioPython: if `gap` is set to a
                // real nucleotide letter, the table wins and the gap rule never fires.
                if let Some(g) = gap {
                    if codon.iter().all(|&b| b.to_ascii_uppercase() == g) {
                        out.push(g as char);
                        continue;
                    }
                }
                return Err(TranslateError::InvalidCodon(upper(codon)));
            }
        }
    }

    Ok(())
}

/// IUPAC-aware complement, case preserving. Unmapped bytes are passed through unchanged,
/// which is what BioPython does.
static COMPLEMENT: [u8; 256] = build_complement();

const fn build_complement() -> [u8; 256] {
    let mut t = [0u8; 256];
    let mut i = 0;
    while i < 256 {
        t[i] = i as u8;
        i += 1;
    }
    // Ambiguity codes complement as their base sets do: R(AG)<->Y(CT), K(GT)<->M(AC),
    // B(CGT)<->V(ACG), D(AGT)<->H(ACT); S(CG) and W(AT) are self-complementary.
    let pairs: &[(u8, u8)] = &[
        (b'A', b'T'), (b'T', b'A'), (b'C', b'G'), (b'G', b'C'), (b'U', b'A'),
        (b'R', b'Y'), (b'Y', b'R'), (b'S', b'S'), (b'W', b'W'),
        (b'K', b'M'), (b'M', b'K'), (b'B', b'V'), (b'V', b'B'),
        (b'D', b'H'), (b'H', b'D'), (b'N', b'N'), (b'X', b'X'),
    ];
    let mut j = 0;
    while j < pairs.len() {
        let (from, to) = pairs[j];
        t[from as usize] = to;
        t[(from + 32) as usize] = to + 32; // lower-case
        j += 1;
    }
    t
}

pub fn reverse_complement(seq: &[u8], out: &mut Vec<u8>) {
    out.clear();
    out.extend(seq.iter().rev().map(|&b| COMPLEMENT[b as usize]));
}
