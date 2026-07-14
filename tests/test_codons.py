"""Exhaustive codon sweep: all 4913 codons x all 27 NCBI tables, against BioPython.

This is cheap (132,651 comparisons, a few seconds) and *total* rather than sampled, so the
generated lookup tables are proven correct rather than spot-checked. If the codegen ever
drifts -- a new BioPython, a botched edit to codon_tables.rs -- this catches it immediately.

We compare with ``on_error="null"`` so that a codon BioPython *refuses* is also part of the
comparison: getting the right answer on the codons it accepts is only half the requirement.
"""

from __future__ import annotations

from itertools import product

import pytest

from conftest import ALL_TABLES, DUAL_CODING_TABLES, IUPAC, IUPAC_VALID, ours

ALL_CODONS = ["".join(c) for c in product(IUPAC, repeat=3)]
X_CODONS = [c for c in ALL_CODONS if "X" in c]


def _biopython(codons, **kw):
    """BioPython's answer per codon, with a refusal represented as ``None``."""
    from Bio.Seq import Seq

    out = []
    for c in codons:
        try:
            out.append(str(Seq(c).translate(**kw)))
        except Exception:
            out.append(None)
    return out


@pytest.mark.parametrize("table", ALL_TABLES)
def test_every_codon_matches_biopython(table):
    expected = _biopython(ALL_CODONS, table=table)
    got = ours(ALL_CODONS, table=table, on_error="null")

    mismatches = [(c, e, g) for c, e, g in zip(ALL_CODONS, expected, got) if e != g]
    assert not mismatches, (
        f"table {table}: {len(mismatches)}/{len(ALL_CODONS)} codons differ from BioPython; "
        f"first 10: {mismatches[:10]}"
    )


@pytest.mark.parametrize("table", ALL_TABLES)
def test_every_codon_matches_biopython_to_stop(table):
    """Same sweep with to_stop, which exercises the stop-detection path."""
    if table in DUAL_CODING_TABLES:
        pytest.skip("to_stop is rejected for dual-coding tables, as in BioPython")

    expected = _biopython(ALL_CODONS, table=table, to_stop=True)
    got = ours(ALL_CODONS, table=table, to_stop=True, on_error="null")
    assert expected == got


@pytest.mark.parametrize("table", ALL_TABLES)
def test_x_bearing_codons(table):
    """The X quirk, swept exhaustively: 17^3 - 16^3 = 817 codons contain an X."""
    expected = _biopython(X_CODONS, table=table)
    got = ours(X_CODONS, table=table, on_error="null")
    assert expected == got

    accepted = [(c, g) for c, g in zip(X_CODONS, got) if g is not None]
    assert accepted, "expected at least some X-bearing codons to resolve"
    for codon, aa in accepted:
        # An accepted X codon can be an amino acid (CTX -> L) or even an ambiguity code
        # (AAX -> X, since AAA/AAG are Lys and AAC/AAT are Asn). What it can never be is a
        # stop -- if the expansion reached a stop codon, the codon would have been refused.
        assert aa != "*", f"{codon} -> {aa!r}: an accepted X codon must never be a stop"


@pytest.mark.parametrize("table", ALL_TABLES)
def test_x_codon_is_accepted_exactly_when_its_expansion_has_no_stop(table):
    """The rule, stated exactly.

    An X-bearing codon is accepted if and only if none of the concrete codons it expands to is
    a stop. That is the whole quirk: 'X' resolves through the expansion table like 'N', but is
    not a letter BioPython will fall back on when the expansion fails -- and the expansion only
    fails when it runs into a stop.
    """
    from Bio.Data import CodonTable

    unamb = CodonTable.unambiguous_dna_by_id[table]
    stops, forward = set(unamb.stop_codons), unamb.forward_table

    got = ours(X_CODONS, table=table, on_error="null")
    for codon, out in zip(X_CODONS, got):
        concrete = [
            "".join(p) for p in product(*(_EXPAND[letter] for letter in codon))
        ]
        # In dual-coding tables the forward table wins, so such a codon is not a stop here.
        hits_a_stop = any(c in stops and c not in forward for c in concrete)
        assert (out is None) == hits_a_stop, (
            f"table {table} codon {codon}: accepted={out is not None} "
            f"but expansion {'does' if hits_a_stop else 'does not'} hit a stop"
        )


def test_x_codons_agree_with_n_codons_wherever_they_are_accepted():
    """Wherever an X codon is accepted, it means exactly what its N spelling means."""
    with_x = [c for c in X_CODONS if set(c) <= set("ACGTX")]
    as_n = [c.replace("X", "N") for c in with_x]

    got_x = ours(with_x, on_error="null")
    got_n = ours(as_n, on_error="null")

    for cx, gx, gn in zip(with_x, got_x, got_n):
        if gx is not None:
            assert gx == gn, f"{cx} -> {gx!r} but its N-spelling gives {gn!r}"

    # ... and where it is refused, the N spelling is nonetheless fine. That asymmetry IS the bug.
    refused = [(cx, gn) for cx, gx, gn in zip(with_x, got_x, got_n) if gx is None]
    assert refused, "expected some X codons to be refused"
    assert all(gn is not None for _, gn in refused)


_EXPAND = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "W": "AT", "S": "CG", "M": "AC", "K": "GT",
    "H": "ACT", "B": "CGT", "V": "ACG", "D": "AGT", "N": "ACGT", "X": "ACGT",
}


def test_valid_input_letters_are_exactly_the_sixteen():
    """Any single letter repeated three times: we must accept exactly what BioPython accepts."""
    from Bio.Seq import Seq

    for byte in range(32, 127):
        ch = chr(byte)
        codon = ch * 3
        try:
            str(Seq(codon).translate(gap=None))
            bio_ok = True
        except Exception:
            bio_ok = False

        we_ok = ours([codon], gap=None, on_error="null")[0] is not None

        assert bio_ok == we_ok, f"codon {codon!r}: BioPython ok={bio_ok}, ours ok={we_ok}"
        if bio_ok:
            assert ch.upper() in IUPAC_VALID, f"{ch!r} accepted but not a valid input letter"

    # 'X' is the letter that makes this interesting: XXX is refused by both.
    assert ours(["XXX"], on_error="null") == [None]
