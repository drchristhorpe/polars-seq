"""Golden tests: every rule and trap from PLAN.md section 4, asserted explicitly.

The parity fuzz in test_parity.py would catch a regression in any of these too, but it would
report it as "sequence 4192 differs". These name the rule that broke.
"""

from __future__ import annotations

import polars as pl
import pytest

from conftest import ours, ours_or_error


class TestBasics:
    def test_simple(self):
        assert ours(["ATGAAATTT"]) == ["MKF"]

    def test_empty_string(self):
        assert ours([""]) == [""]

    def test_null_in_null_out(self):
        assert ours([None, "ATG", None]) == [None, "M", None]

    def test_lowercase_is_uppercased(self):
        assert ours(["atgaaa"]) == ["MK"]

    def test_mixed_case(self):
        assert ours(["AtGaAa"]) == ["MK"]

    def test_rna_u_accepted(self):
        assert ours(["AUGAAA"]) == ["MK"]

    def test_stop_codon(self):
        assert ours(["TGA"]) == ["*"]


class TestAmbiguity:
    """The subtle part. See PLAN.md 4.2."""

    def test_unresolvable_but_unique_amino_acid(self):
        # All four GGx codons are glycine, so GGN is unambiguously G.
        assert ours(["GGN"]) == ["G"]

    def test_all_expansions_are_stops(self):
        # TAR -> {TAA, TAG}, both stops.
        assert ours(["TAR"]) == ["*"]

    def test_mixed_stop_and_amino_acid_is_X(self):
        # TAN -> {TAA*, TAC=Y, TAG*, TAT=Y}: unresolvable.
        assert ours(["TAN"]) == ["X"]

    def test_nnn_is_X(self):
        assert ours(["NNN"]) == ["X"]

    def test_asx_ambiguity_code(self):
        # THE trap: RAT -> {GAT=Asp, AAT=Asn} -> 'B' (Asx), *not* 'X'.
        assert ours(["RAT"]) == ["B"]

    def test_glx_ambiguity_code(self):
        # SAA -> {CAA=Gln, GAA=Glu} -> 'Z' (Glx).
        assert ours(["SAA"]) == ["Z"]

    def test_ambiguous_but_single_amino_acid(self):
        assert ours(["AAY"]) == ["N"]  # AAT, AAC both Asn


class TestTheXQuirk:
    """'X' expands to GATC like 'N', but is NOT a valid *input* letter.

    So an X-bearing codon is accepted precisely when it still resolves to a single amino acid,
    and rejected otherwise. It can never produce 'X' or a stop. This is the sharpest edge in
    the whole of `Seq.translate()`; see README, "The X-codon quirk".
    """

    def test_x_codon_resolving_to_one_amino_acid_is_accepted(self):
        # All four CTx codons are Leu, so the X does not actually cost us any information.
        assert ours(["CTX"]) == ["L"]
        assert ours(["CTTCTX"]) == ["LL"]

    def test_x_codon_that_does_not_resolve_is_invalid(self):
        # XXX spans every codon -> amino acids and stops -> cannot resolve -> rejected,
        # NOT reported as 'X'. Contrast with NNN, which *is* 'X'.
        assert ours_or_error("XXX") == "ERROR"
        assert ours(["NNN"]) == ["X"]

    def test_x_codon_mixing_stop_and_amino_acid_is_invalid(self):
        # TAX -> {TAA*, TAC=Y, TAG*, TAT=Y}. Same expansion as TAN, but TAN -> 'X' and
        # TAX -> error, purely because of which letter spelled the ambiguity.
        assert ours_or_error("TAX") == "ERROR"
        assert ours(["TAN"]) == ["X"]

    def test_x_and_n_agree_whenever_the_codon_resolves(self):
        assert ours(["CTX"]) == ours(["CTN"]) == ["L"]
        assert ours(["GGX"]) == ours(["GGN"]) == ["G"]


class TestAlphabet:
    def test_junk_letter_is_invalid(self):
        assert ours_or_error("ATZ") == "ERROR"

    def test_digit_is_invalid(self):
        assert ours_or_error("AT1") == "ERROR"

    def test_error_message_names_the_codon(self):
        with pytest.raises(pl.exceptions.ComputeError, match="Codon 'XXX' is invalid"):
            ours(["XXX"])

    def test_error_message_is_uppercased(self):
        with pytest.raises(pl.exceptions.ComputeError, match="Codon 'ATZ' is invalid"):
            ours(["atz"])


class TestGaps:
    def test_full_gap_codon(self):
        assert ours(["---"]) == ["-"]

    def test_multiple_gap_codons(self):
        assert ours(["------"]) == ["--"]

    def test_gap_mixed_with_bases_is_invalid(self):
        assert ours_or_error("A-A") == "ERROR"
        assert ours_or_error("--A") == "ERROR"

    def test_gaps_interleaved_with_codons(self):
        assert ours(["ATG---AAA"]) == ["M-K"]

    def test_custom_gap_char(self):
        assert ours(["ATG...AAA"], gap=".") == ["M.K"]

    def test_gap_none_makes_gaps_invalid(self):
        assert ours_or_error("---", gap=None) == "ERROR"

    def test_table_wins_over_gap(self):
        # If gap is set to a real nucleotide letter, the codon table still wins.
        assert ours(["AAA"], gap="A") == ["K"]


class TestPartialCodons:
    """Trailing partial codons are silently dropped. PLAN.md 4.5."""

    def test_one_trailing_base(self):
        assert ours(["ATGA"]) == ["M"]

    def test_two_trailing_bases(self):
        assert ours(["ATGAA"]) == ["M"]

    def test_shorter_than_a_codon(self):
        assert ours(["AT"]) == [""]


class TestStopHandling:
    def test_to_stop_excludes_the_stop(self):
        assert ours(["ATGAAATAAGGG"], to_stop=True) == ["MK"]

    def test_to_stop_with_no_stop_present(self):
        assert ours(["ATGAAAGGG"], to_stop=True) == ["MKG"]

    def test_custom_stop_symbol(self):
        assert ours(["ATGAAATAAGGG"], stop_symbol="#") == ["MK#G"]

    def test_stop_symbol_ignored_when_to_stop(self):
        assert ours(["ATGAAATAAGGG"], stop_symbol="#", to_stop=True) == ["MK"]


class TestCds:
    """PLAN.md 4.7. The ordering of these checks is itself part of the behaviour."""

    def test_valid_cds(self):
        assert ours(["ATGCCCTAG"], cds=True) == ["MP"]

    def test_terminal_stop_is_excluded(self):
        assert ours(["ATGTAA"], cds=True) == ["M"]

    def test_alternative_start_becomes_M(self):
        # TTG encodes Leu, but as a start codon under cds it is reported as Met.
        assert ours(["TTGAAATAA"], cds=True) == ["MK"]
        assert ours(["TTGAAATAA"], cds=False) == ["LK*"]

    def test_bad_start(self):
        with pytest.raises(pl.exceptions.ComputeError, match="First codon 'AAA' is not a start codon"):
            ours(["AAACCCTAG"], cds=True)

    def test_missing_stop(self):
        with pytest.raises(pl.exceptions.ComputeError, match="Final codon 'ATG' is not a stop codon"):
            ours(["ATG"], cds=True)

    def test_internal_stop(self):
        with pytest.raises(pl.exceptions.ComputeError, match="Extra in frame stop codon found"):
            ours(["ATGCCCTAGCCCTAG"], cds=True)

    def test_length_not_multiple_of_three(self):
        with pytest.raises(pl.exceptions.ComputeError, match="not a multiple of three"):
            ours(["ATGCCCTAGA"], cds=True)

    def test_start_is_checked_before_length(self):
        # 'AAAC' fails both rules; BioPython reports the start-codon one.
        with pytest.raises(pl.exceptions.ComputeError, match="is not a start codon"):
            ours(["AAAC"], cds=True)

    def test_to_stop_is_ignored_under_cds(self):
        assert ours(["ATGCCCTAG"], cds=True, to_stop=True) == ["MP"]


class TestTables:
    def test_by_id(self):
        assert ours(["ATGAGATAA"], table=2) == ["M**"]  # AGA is a stop in vertebrate mito

    def test_by_name(self):
        assert ours(["ATGAGATAA"], table="Vertebrate Mitochondrial") == ["M**"]

    def test_by_alias(self):
        assert ours(["ATGAGATAA"], table="SGC1") == ["M**"]

    def test_name_is_case_insensitive(self):
        assert ours(["ATGAGATAA"], table="vertebrate mitochondrial") == ["M**"]

    def test_stringified_id(self):
        assert ours(["ATGAGATAA"], table="2") == ["M**"]

    def test_standard_is_the_default(self):
        assert ours(["ATGAGATAA"]) == ours(["ATGAGATAA"], table=1) == ["MR*"]

    def test_unknown_id_raises_valueerror_not_computeerror(self):
        with pytest.raises(ValueError, match="unknown codon table id"):
            ours(["ATG"], table=99)

    def test_unknown_name(self):
        with pytest.raises(ValueError, match="unknown codon table name"):
            ours(["ATG"], table="Klingon Mitochondrial")


class TestDualCodingTables:
    """Tables 27/28/31: codons that are both a stop and an amino acid. PLAN.md 4.8."""

    def test_dual_coding_codon_translates_as_amino_acid(self):
        # Table 28: TAA/TAG/TGA are all also amino acids.
        assert ours(["ATGTAAAAA"], table=28) == ["MQK"]

    def test_to_stop_is_rejected(self):
        with pytest.raises(ValueError, match="cannot use to_stop=True"):
            ours(["ATGTAA"], table=28, to_stop=True)

    def test_warns(self):
        from polars_seq import AmbiguousStopCodonWarning

        with pytest.warns(AmbiguousStopCodonWarning):
            pl.col("s").seq.translate(table=28)

    def test_dual_coding_codon_still_valid_as_cds_terminator(self):
        # It is an amino acid mid-sequence, yet still a legal terminal stop.
        assert ours(["ATGAAATAA"], table=28, cds=True) == ["MK"]


class TestArgumentValidation:
    def test_multi_char_stop_symbol(self):
        with pytest.raises(ValueError, match="single character"):
            ours(["ATG"], stop_symbol="**")

    def test_multi_char_gap(self):
        with pytest.raises(ValueError, match="single character"):
            ours(["ATG"], gap="--")

    def test_non_ascii_stop_symbol(self):
        with pytest.raises(ValueError, match="ASCII"):
            ours(["ATG"], stop_symbol="†")

    def test_bad_on_error(self):
        with pytest.raises(ValueError, match="on_error"):
            ours(["ATG"], on_error="explode")


class TestOnError:
    """Our one behavioural extension over BioPython."""

    def test_null_on_error_isolates_the_bad_row(self):
        assert ours(["ATG", "XXX", "AAA"], on_error="null") == ["M", None, "K"]

    def test_raise_is_the_default(self):
        with pytest.raises(pl.exceptions.ComputeError):
            ours(["ATG", "XXX", "AAA"])

    def test_null_on_error_with_cds(self):
        assert ours(["ATGTAA", "AAATAA"], cds=True, on_error="null") == ["M", None]


class TestReverseComplement:
    def test_simple(self):
        assert _rc(["ATGC"]) == ["GCAT"]

    def test_iupac_codes(self):
        assert _rc(["RYKMBVDHNSW"]) == ["WSNDHBVKMRY"]

    def test_case_preserved(self):
        assert _rc(["atgc"]) == ["gcat"]

    def test_null(self):
        assert _rc([None]) == [None]

    def test_composes_with_translate(self):
        # revcomp("TTTCAT") == "ATGAAA" -> "MK"
        df = pl.DataFrame({"s": ["TTTCAT"]})
        got = df.select(pl.col("s").seq.reverse_complement().seq.translate()).item()
        assert got == "MK"


def _rc(seqs):
    return (
        pl.DataFrame({"s": seqs}, schema={"s": pl.String})
        .select(pl.col("s").seq.reverse_complement())
        .to_series()
        .to_list()
    )
