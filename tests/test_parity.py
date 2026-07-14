"""Differential fuzz against BioPython. This is the test the whole project rests on.

It asserts not only that we produce the same protein as BioPython, but that we *fail on
exactly the same inputs* -- a plugin that quietly accepted something BioPython rejects would
be just as wrong as one that mistranslates.

Seeded, so any failure reproduces.
"""

from __future__ import annotations

import random

import pytest

from conftest import ALL_TABLES, DUAL_CODING_TABLES, IUPAC, UNAMBIGUOUS, bio_or_error, ours_or_error

# Alphabets chosen to exercise each distinct code path in the kernel.
ALPHABETS = {
    "unambiguous": UNAMBIGUOUS,
    "iupac": IUPAC,
    "rna": "ACGU",
    "gapped": "ACGT---",
    "gappy": "ACGT-",
    "lowercase": "acgt",
    "mixed_case": "AcGtRyN",
    "with_invalid": "ACGTXZ?*",  # X, Z, ?, * are not valid nucleotide letters
    "all_gaps": "-",
    "ambiguous_heavy": "RYWSMKHBVDN",
}


def _seqs(rng: random.Random, alphabet: str, n: int) -> list[str]:
    out = []
    for _ in range(n):
        # Deliberately include lengths that are not multiples of three, and the empty string.
        length = rng.choice([0, 1, 2, 3, 4, 5, 6, 9, 12, 17, 30, 31, 99, 100])
        out.append("".join(rng.choice(alphabet) for _ in range(length)))
    return out


@pytest.mark.parametrize("alphabet_name", sorted(ALPHABETS))
def test_parity_across_alphabets(alphabet_name, rng):
    alphabet = ALPHABETS[alphabet_name]
    for seq in _seqs(rng, alphabet, 400):
        assert ours_or_error(seq) == bio_or_error(seq), f"seq={seq!r}"


@pytest.mark.parametrize("table", ALL_TABLES)
def test_parity_across_all_tables(table, rng):
    for seq in _seqs(rng, IUPAC, 150):
        assert ours_or_error(seq, table=table) == bio_or_error(seq, table=table), (
            f"table={table} seq={seq!r}"
        )


@pytest.mark.parametrize("to_stop", [False, True])
@pytest.mark.parametrize("stop_symbol", ["*", "#"])
def test_parity_stop_options(to_stop, stop_symbol, rng):
    for seq in _seqs(rng, UNAMBIGUOUS, 300):
        kw = {"to_stop": to_stop, "stop_symbol": stop_symbol}
        assert ours_or_error(seq, **kw) == bio_or_error(seq, **kw), f"seq={seq!r} {kw}"


def test_parity_cds(rng):
    """cds=True has four rejection rules whose *ordering* is observable. Fuzz it hard."""
    seqs = []
    # Mostly-valid CDSs, so we actually exercise the success path and the internal-stop rule
    # rather than just bouncing off the start-codon check.
    starts = ["ATG", "TTG", "CTG", "GTG", "AAA", "ATT"]
    stops = ["TAA", "TAG", "TGA", "AAA", "CCC"]
    for _ in range(2000):
        body = "".join(rng.choice(UNAMBIGUOUS) for _ in range(3 * rng.randint(0, 8)))
        seq = rng.choice(starts) + body + rng.choice(stops)
        if rng.random() < 0.15:  # sometimes break the frame
            seq += rng.choice(UNAMBIGUOUS)
        seqs.append(seq)
    seqs += ["", "A", "AT", "ATG", "ATGTAA", "TAA"]

    for seq in seqs:
        assert ours_or_error(seq, cds=True) == bio_or_error(seq, cds=True), f"seq={seq!r}"


def test_parity_gap_options(rng):
    for gap in ["-", ".", None]:
        for seq in _seqs(rng, "ACGT-.", 300):
            assert ours_or_error(seq, gap=gap) == bio_or_error(seq, gap=gap), (
                f"seq={seq!r} gap={gap!r}"
            )


@pytest.mark.parametrize("table", DUAL_CODING_TABLES)
def test_parity_dual_coding_tables(table, rng):
    """Tables where a codon is both a stop and an amino acid."""
    for seq in _seqs(rng, UNAMBIGUOUS, 300):
        assert ours_or_error(seq, table=table) == bio_or_error(seq, table=table), (
            f"table={table} seq={seq!r}"
        )
        # ... and they must still work as CDS terminators.
        assert ours_or_error(seq, table=table, cds=True) == bio_or_error(seq, table=table, cds=True)


def test_parity_combinatorial(rng):
    """Cross-product of the option space on a shared pool of sequences."""
    seqs = _seqs(rng, IUPAC, 40) + _seqs(rng, "ACGT-", 40) + _seqs(rng, UNAMBIGUOUS, 40)
    for table in (1, 2, 6, 11, 16, 25, 33):
        for to_stop in (False, True):
            for cds in (False, True):
                for gap in ("-", None):
                    for stop_symbol in ("*", "@"):
                        kw = dict(
                            table=table,
                            to_stop=to_stop,
                            cds=cds,
                            gap=gap,
                            stop_symbol=stop_symbol,
                        )
                        for seq in seqs:
                            assert ours_or_error(seq, **kw) == bio_or_error(seq, **kw), (
                                f"seq={seq!r} {kw}"
                            )


def test_parity_reverse_complement(rng):
    from Bio.Seq import Seq

    import polars as pl

    seqs = _seqs(rng, IUPAC, 200) + _seqs(rng, "acgtryn", 100)
    got = (
        pl.DataFrame({"s": seqs}, schema={"s": pl.String})
        .select(pl.col("s").seq.reverse_complement())
        .to_series()
        .to_list()
    )
    expected = [str(Seq(s).reverse_complement()) for s in seqs]
    assert got == expected
