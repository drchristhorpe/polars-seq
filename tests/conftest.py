from __future__ import annotations

import warnings

import polars as pl
import pytest

import polars_seq  # noqa: F401  -- registers the .seq namespace

# BioPython warns about partial codons and dual-coding tables. We are comparing behaviour,
# not warnings, so silence them across the suite.
warnings.filterwarnings("ignore")

# The 17 letters that can appear in a codon BioPython will look up. 'X' is the odd one out:
# it expands to GATC like 'N', but is *not* in BioPython's set of valid input letters, so an
# X-bearing codon is accepted only when it still resolves to one amino acid (CTX -> L) and is
# rejected otherwise (XXX -> invalid). See README, "The X-codon quirk".
IUPAC = "ACGTURYWSMKHBVDNX"
# The 16 that are valid input letters in their own right.
IUPAC_VALID = "ACGTURYWSMKHBVDN"
UNAMBIGUOUS = "ACGT"

# All 27 NCBI tables.
ALL_TABLES = [1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 13, 14, 15, 16, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33]
DUAL_CODING_TABLES = [27, 28, 31]


def ours(seqs, **kwargs) -> list[str | None]:
    """Translate via the plugin."""
    return (
        pl.DataFrame({"s": seqs}, schema={"s": pl.String})
        .select(pl.col("s").seq.translate(**kwargs))
        .to_series()
        .to_list()
    )


def theirs(seqs, **kwargs) -> list[str | None]:
    """Translate via BioPython, row by row."""
    from Bio.Seq import Seq

    out = []
    for s in seqs:
        out.append(None if s is None else str(Seq(s).translate(**kwargs)))
    return out


def bio_or_error(seq: str, **kwargs):
    """BioPython's answer, or the sentinel ``ERROR`` if it refuses."""
    from Bio.Seq import Seq

    try:
        return str(Seq(seq).translate(**kwargs))
    except Exception:
        return "ERROR"


def ours_or_error(seq: str, **kwargs):
    """Our answer, or the sentinel ``ERROR`` if we refuse."""
    try:
        return ours([seq], **kwargs)[0]
    except Exception:
        return "ERROR"


@pytest.fixture(scope="session")
def rng():
    import random

    return random.Random(20260714)
