"""Integration with the Polars engine itself, as opposed to the translation semantics."""

from __future__ import annotations

import polars as pl
import pytest

import polars_seq
from conftest import ours


def test_lazy_frame():
    lf = pl.LazyFrame({"dna": ["ATGAAA", "ATGTTT"]})
    got = lf.with_columns(protein=pl.col("dna").seq.translate()).collect()
    assert got["protein"].to_list() == ["MK", "MF"]


def test_streaming_engine():
    lf = pl.LazyFrame({"dna": ["ATGAAA"] * 1000})
    got = lf.with_columns(protein=pl.col("dna").seq.translate()).collect(engine="streaming")
    assert got["protein"].to_list() == ["MK"] * 1000


def test_empty_frame():
    df = pl.DataFrame({"dna": []}, schema={"dna": pl.String})
    got = df.select(pl.col("dna").seq.translate())
    assert got.height == 0
    assert got.schema["dna"] == pl.String


def test_all_null_column():
    assert ours([None, None]) == [None, None]


def test_output_dtype_is_string():
    df = pl.DataFrame({"dna": ["ATG"]})
    assert df.select(pl.col("dna").seq.translate()).schema["dna"] == pl.String


def test_chunked_series():
    """A Series assembled from several chunks must translate identically to a contiguous one."""
    a = pl.Series("dna", ["ATGAAA", "ATGTTT"])
    b = pl.Series("dna", ["ATGCCC", None])
    chunked = a.append(b)
    assert chunked.n_chunks() > 1

    got = pl.DataFrame({"dna": chunked}).select(pl.col("dna").seq.translate())
    assert got.to_series().to_list() == ["MK", "MF", "MP", None]


def test_sliced_series():
    """A slice does not start at offset 0 in the underlying buffer -- easy thing to get wrong."""
    df = pl.DataFrame({"dna": ["XXX", "ATGAAA", "ATGTTT", "XXX"]})
    got = df.slice(1, 2).select(pl.col("dna").seq.translate())
    assert got.to_series().to_list() == ["MK", "MF"]


def test_group_by_agg():
    df = pl.DataFrame({"g": ["a", "a", "b"], "dna": ["ATG", "AAA", "TTT"]})
    got = (
        df.group_by("g")
        .agg(pl.col("dna").seq.translate().alias("prot"))
        .sort("g")
    )
    assert got["prot"].to_list() == [["M", "K"], ["F"]]


def test_large_frame_multithreaded():
    n = 200_000
    df = pl.DataFrame({"dna": ["ATGAAATTTGGGTAA"] * n})
    got = df.with_columns(protein=pl.col("dna").seq.translate())
    assert got["protein"].n_unique() == 1
    assert got["protein"][0] == "MKFG*"
    assert got.height == n


def test_expression_reuse():
    """The same expression object used twice in one query."""
    e = pl.col("dna").seq.translate()
    df = pl.DataFrame({"dna": ["ATGAAA"]})
    got = df.select(a=e, b=e)
    assert got.row(0) == ("MK", "MK")


def test_chained_with_other_expressions():
    df = pl.DataFrame({"dna": ["atgaaataaggg"]})
    got = df.select(
        pl.col("dna").str.to_uppercase().seq.translate(to_stop=True).str.len_chars().alias("n")
    )
    assert got.item() == 2


def test_error_message_includes_row_index():
    df = pl.DataFrame({"dna": ["ATG", "ATG", "XXX"]})
    with pytest.raises(pl.exceptions.ComputeError, match=r"row 2"):
        df.select(pl.col("dna").seq.translate())


def test_non_string_column_raises():
    df = pl.DataFrame({"dna": [1, 2, 3]})
    with pytest.raises(Exception):
        df.select(pl.col("dna").seq.translate())


def test_codon_tables_helper():
    tables = polars_seq.codon_tables()
    assert tables.height == 27
    assert tables.filter(pl.col("id") == 1)["name"].item() == "Standard"
    assert tables.filter(pl.col("id") == 28)["dual_coding"].item().to_list() == ["TAA", "TAG", "TGA"]


def test_six_frame_translation():
    """The use case reverse_complement exists for."""
    df = pl.DataFrame({"dna": ["ATGAAATTTGGGCCC"]})
    got = df.select(
        **{
            f"fwd{i}": pl.col("dna").str.slice(i).seq.translate() for i in range(3)
        },
        **{
            f"rev{i}": pl.col("dna").seq.reverse_complement().str.slice(i).seq.translate()
            for i in range(3)
        },
    )
    assert got["fwd0"].item() == "MKFGP"
    assert got.width == 6
    assert all(got[c].item() is not None for c in got.columns)
