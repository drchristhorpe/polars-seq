"""Produce validation artefacts into `tmp/` so the BioPython-parity claim can be inspected.

Run:  uv run python tools/validate.py

Writes:
  tmp/01_codon_sweep.csv        all 4096 codons x 27 tables, ours vs BioPython
  tmp/02_golden_cases.csv       the documented edge cases, ours vs BioPython
  tmp/03_real_data_parity.csv   every sequence in test_data/*.parquet, ours vs BioPython
  tmp/04_real_data_output.parquet   the translated real dataset
  tmp/05_benchmark.txt          throughput vs BioPython
  tmp/00_SUMMARY.txt            read this one first
"""

from __future__ import annotations

import io
import sys
import time
import warnings
from itertools import product
from pathlib import Path

import polars as pl
from Bio.Seq import Seq

import polars_seq  # noqa: F401  -- registers .seq

warnings.simplefilter("ignore")

ROOT = Path(__file__).resolve().parent.parent
TMP = ROOT / "tmp"
TMP.mkdir(exist_ok=True)

IUPAC = "ACGTURYWSMKHBVDNX"  # 17: 'X' included -- see README, "The X-codon quirk"
ALL_TABLES = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]

summary = io.StringIO()


def say(msg: str = "") -> None:
    print(msg)
    summary.write(msg + "\n")


REFUSED = "<refused>"  # BioPython raised / we raised -- compared like any other value


def bio(seq, **kw):
    try:
        return str(Seq(seq).translate(**kw))
    except Exception:
        return REFUSED


def ours(seqs, **kw):
    """Translate with on_error='null', so a refusal is a comparable value rather than a crash."""
    out = (
        pl.DataFrame({"s": seqs}, schema={"s": pl.String})
        .select(pl.col("s").seq.translate(on_error="null", **kw))
        .to_series()
        .to_list()
    )
    return [REFUSED if v is None else v for v in out]


def ours_one(seq, **kw):
    return ours([seq], **kw)[0]


# --------------------------------------------------------------------------------------
say("=" * 78)
say("polars-seq validation report")
say("=" * 78)
say(f"python      {sys.version.split()[0]}")
say(f"polars      {pl.__version__}")
import Bio  # noqa: E402

say(f"biopython   {Bio.__version__}")
say(f"polars_seq  {polars_seq.__version__}")
say()

# --------------------------------------------------------------------------------------
# 1. Exhaustive codon sweep
# --------------------------------------------------------------------------------------
say("-" * 78)
say("1. EXHAUSTIVE CODON SWEEP -- every codon, every table")
say("-" * 78)

codons = ["".join(c) for c in product(IUPAC, repeat=3)]
rows = []
total = mismatch = 0
for table in ALL_TABLES:
    got = ours(codons, table=table)
    exp = [bio(c, table=table) for c in codons]
    for c, g, e in zip(codons, got, exp):
        total += 1
        ok = g == e
        if not ok:
            mismatch += 1
        rows.append({"table": table, "codon": c, "polars_seq": g, "biopython": e, "match": ok})

sweep = pl.DataFrame(rows)
sweep.write_csv(TMP / "01_codon_sweep.csv")
say(f"  codons compared : {total:,}  ({len(codons):,} codons x {len(ALL_TABLES)} tables)")
say(f"  mismatches      : {mismatch:,}")
say(f"  result          : {'PASS' if mismatch == 0 else 'FAIL'}")
say(f"  -> tmp/01_codon_sweep.csv")
say()

# --------------------------------------------------------------------------------------
# 2. Golden edge cases
# --------------------------------------------------------------------------------------
say("-" * 78)
say("2. DOCUMENTED EDGE CASES")
say("-" * 78)

cases = [
    ("ATGAAATTT", {}, "basic translation"),
    ("atgaaa", {}, "lower-case input is uppercased"),
    ("AUGAAA", {}, "RNA (U) accepted"),
    ("GGN", {}, "ambiguous codon, single amino acid"),
    ("TAR", {}, "ambiguous codon, all expansions are stops"),
    ("TAN", {}, "ambiguous codon, mixed stop/aa -> X"),
    ("RAT", {}, "THE TRAP: {Asp,Asn} -> B, not X"),
    ("SAA", {}, "{Glu,Gln} -> Z"),
    ("NNN", {}, "fully ambiguous -> X"),
    # --- the X-codon quirk: accepted iff the expansion hits no stop codon ---
    ("CTX", {}, "X QUIRK: expansion has no stop -> resolves, like CTN -> L"),
    ("CTTCTX", {}, "X QUIRK: the sequence the fuzz caught the bug on"),
    ("AAX", {}, "X QUIRK: several aa, no stop -> X (an *amino-acid* ambiguity code)"),
    ("AAN", {}, "  ... and the N spelling agrees"),
    ("TAX", {}, "X QUIRK: expansion hits a stop -> ERROR"),
    ("TAN", {}, "  ... whereas the N spelling degrades to X"),
    ("XXX", {}, "X QUIRK: expansion hits a stop -> ERROR"),
    ("---", {}, "whole gap codon -> one gap char"),
    ("A-A", {}, "partial gap codon -> error"),
    ("ATG---AAA", {}, "gaps interleaved"),
    ("ATGAA", {}, "trailing partial codon dropped"),
    ("AT", {}, "shorter than a codon"),
    ("", {}, "empty string"),
    ("ATGAAATAAGGG", {"to_stop": True}, "to_stop excludes the stop"),
    ("ATGAAATAAGGG", {"stop_symbol": "#"}, "custom stop symbol"),
    ("ATGCCCTAG", {"cds": True}, "valid CDS"),
    ("TTGAAATAA", {"cds": True}, "alternative start -> M"),
    ("AAACCCTAG", {"cds": True}, "bad start codon -> error"),
    ("ATGCCCTAGCCCTAG", {"cds": True}, "internal stop -> error"),
    ("ATG", {"cds": True}, "no terminal stop -> error"),
    ("ATGAGATAA", {"table": 2}, "vertebrate mito: AGA is a stop"),
    ("ATGAGATAA", {"table": "SGC1"}, "same, selected by alias"),
    ("ATGTAAAAA", {"table": 28}, "dual-coding table: TAA is an amino acid"),
    ("ATGAAATAA", {"table": 28, "cds": True}, "dual-coding codon still ends a CDS"),
]

rows = []
bad = 0
for seq, kw, note in cases:
    g, e = ours_one(seq, **kw), bio(seq, **kw)
    ok = g == e
    if not ok:
        bad += 1
    rows.append(
        {
            "sequence": seq,
            "options": str(kw) if kw else "",
            "polars_seq": g,
            "biopython": e,
            "match": ok,
            "note": note,
        }
    )
golden = pl.DataFrame(rows)
golden.write_csv(TMP / "02_golden_cases.csv")
with pl.Config(fmt_str_lengths=60, tbl_rows=-1, tbl_width_chars=200):
    say(str(golden.select("sequence", "options", "polars_seq", "biopython", "match")))
say(f"  mismatches: {bad}   result: {'PASS' if bad == 0 else 'FAIL'}")
say(f"  -> tmp/02_golden_cases.csv")
say()

# --------------------------------------------------------------------------------------
# 3 & 4. Real data
# --------------------------------------------------------------------------------------
say("-" * 78)
say("3. REAL DATA -- test_data/*.parquet")
say("-" * 78)

parquets = sorted((ROOT / "test_data").glob("*.parquet"))
if not parquets:
    say("  no parquet found in test_data/ -- skipping")
else:
    src = parquets[0]
    df = pl.read_parquet(src)
    seqs = df["sequence"].to_list()
    say(f"  file       : {src.name}")
    say(f"  rows       : {df.height:,}")
    say(f"  columns    : {df.columns}")
    lengths = df.select(pl.col("sequence").str.len_chars()).to_series()
    say(f"  seq length : min={lengths.min()} max={lengths.max()} mean={lengths.mean():.1f}")
    say()

    # Translate with both, three ways, and compare every single row.
    variants = {
        "default": {},
        "to_stop": {"to_stop": True},
        "table11_bacterial": {"table": 11},
    }
    for name, kw in variants.items():
        t0 = time.perf_counter()
        got = df.select(pl.col("sequence").seq.translate(**kw)).to_series().to_list()
        t_ours = time.perf_counter() - t0

        t0 = time.perf_counter()
        exp = [bio(s, **kw) for s in seqs]
        t_bio = time.perf_counter() - t0

        mism = [
            {"row": i, "sequence": s, "polars_seq": g, "biopython": e}
            for i, (s, g, e) in enumerate(zip(seqs, got, exp))
            if g != e
        ]
        say(
            f"  [{name:18}] rows={len(seqs):,}  mismatches={len(mism):,}  "
            f"{'PASS' if not mism else 'FAIL'}   "
            f"(polars_seq {t_ours:6.3f}s vs biopython {t_bio:7.3f}s -> {t_bio / t_ours:5.1f}x)"
        )
        if name == "default":
            pl.DataFrame(
                {
                    "sequence": seqs,
                    "polars_seq": got,
                    "biopython": exp,
                    "match": [g == e for g, e in zip(got, exp)],
                }
            ).write_csv(TMP / "03_real_data_parity.csv")
            if mism:
                pl.DataFrame(mism).write_csv(TMP / "03b_MISMATCHES.csv")

    say()
    out = df.with_columns(
        protein=pl.col("sequence").seq.translate(on_error="null"),
        protein_to_stop=pl.col("sequence").seq.translate(to_stop=True, on_error="null"),
        protein_revcomp=pl.col("sequence")
        .seq.reverse_complement()
        .seq.translate(on_error="null"),
    )
    n_refused = out["protein"].null_count()
    say(f"  sequences refused by translation: {n_refused:,} of {out.height:,}")
    out.write_parquet(TMP / "04_real_data_output.parquet")
    say(f"  -> tmp/03_real_data_parity.csv  (every row, side by side)")
    say(f"  -> tmp/04_real_data_output.parquet")
    say()
    with pl.Config(fmt_str_lengths=42, tbl_width_chars=220):
        say("  sample of the translated output:")
        say(str(out.head(8)))
    say()

    # --------------------------------------------------------------------------------
    # 5. Benchmark
    # --------------------------------------------------------------------------------
    say("-" * 78)
    say("4. BENCHMARK -- polars-seq vs BioPython via map_elements")
    say("-" * 78)

    bench = io.StringIO()
    for n in (10_000, 100_000, df.height):
        sub = df.head(n)
        subs = sub["sequence"].to_list()

        t0 = time.perf_counter()
        for _ in range(3):
            sub.select(pl.col("sequence").seq.translate(on_error="null"))
        t_ours = (time.perf_counter() - t0) / 3

        t0 = time.perf_counter()
        [str(Seq(s).translate()) for s in subs]
        t_bio = time.perf_counter() - t0

        line = (
            f"  n={n:>9,}   polars-seq {t_ours * 1000:9.1f} ms   "
            f"biopython {t_bio * 1000:10.1f} ms   speed-up {t_bio / t_ours:6.1f}x   "
            f"({n / t_ours / 1e6:.1f}M seq/s)"
        )
        say(line)
        bench.write(line + "\n")

    (TMP / "05_benchmark.txt").write_text(bench.getvalue())
    say(f"  -> tmp/05_benchmark.txt")

say()
say("=" * 78)
say("DONE")
say("=" * 78)

(TMP / "00_SUMMARY.txt").write_text(summary.getvalue())
print(f"\nSummary written to {TMP / '00_SUMMARY.txt'}")
