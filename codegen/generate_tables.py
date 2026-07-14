"""Generate `src/codon_tables.rs` and `python/polars_seq/_tables.py` from BioPython.

Rather than re-implement BioPython's ambiguous-codon resolver in Rust (fiddly, and the
source of the classic `RAT -> X` bug when the right answer is `B`), we use BioPython as
an *oracle at build time*: enumerate every codon over the 17 IUPAC nucleotide letters,
ask BioPython what it does, and bake the answers into a static Rust table.

The generated files are committed, so BioPython is a build-/test-time dependency only --
never a runtime one.

Run with:  uv run --with biopython python codegen/generate_tables.py
"""

from __future__ import annotations

import warnings
from pathlib import Path

import Bio
from Bio.Data import CodonTable
from Bio.Seq import _translate_str

warnings.simplefilter("ignore")  # dual-coding tables warn; we handle them explicitly

# The 17 letters that can appear in a codon BioPython is willing to look up. This is the union
# of the nucleotide-expansion keys, which is *not* the same as the set of letters BioPython
# considers "valid input" -- 'X' is a key (it expands to GATC, like 'N') but is NOT in
# `_ambiguous_dna_letters | _ambiguous_rna_letters`.
#
# The consequence is genuinely surprising, and is why this alphabet has 17 letters and not 16:
#   CTX -> 'L'      (all four CTx codons are Leu, so the codon resolves despite the X)
#   XXX -> invalid  (does not resolve; the fallback validity check then rejects the 'X')
# The rule, exactly: an X-bearing codon is accepted iff none of the codons it expands to is a
# stop. When accepted it means whatever the 'N' spelling means (so AAX -> 'X', an *amino-acid*
# ambiguity code). It can never come out as a stop. See PLAN.md section 4.3.
ALPHABET = "ACGTURYWSMKHBVDNX"
assert len(ALPHABET) == len(set(ALPHABET)) == 17

VALID_INPUT_LETTERS = set("ACGTURYWSMKHBVDN")  # note: no 'X'

# Sentinels stored in the `code` array, chosen to not collide with any amino-acid byte.
STOP = 0  # a true stop codon -> emit the user's stop_symbol
INVALID = 1  # BioPython refuses this codon -> raise "Codon '...' is invalid"

N_CODONS = 17**3  # 4913
N_WORDS = (N_CODONS + 63) // 64  # 77

ROOT = Path(__file__).resolve().parent.parent


def codon_index(codon: str) -> int:
    """Map a 3-letter codon to 0..4912 (base 17)."""
    a, b, c = (ALPHABET.index(x) for x in codon)
    return a * 289 + b * 17 + c


def build_table(table_id: int) -> dict:
    """Resolve all 4913 codons for one NCBI table, using BioPython as the oracle.

    We do not reimplement the resolver -- we *ask the real one*, via `_translate_str`, which
    is the exact function `Seq.translate()` calls. Whatever it says, including which codons it
    refuses, is what we bake in.
    """
    amb = CodonTable.ambiguous_generic_by_id[table_id]
    stop_codons = set(amb.stop_codons)
    start_codons = set(amb.start_codons)

    code = [0] * N_CODONS
    is_stop = [False] * N_CODONS
    is_start = [False] * N_CODONS

    for a in ALPHABET:
        for b in ALPHABET:
            for c in ALPHABET:
                codon = a + b + c
                idx = codon_index(codon)

                try:
                    # The oracle. '*' is unambiguous here: no amino-acid code is '*'.
                    out = _translate_str(codon, amb, stop_symbol="*")
                except CodonTable.TranslationError:
                    value = INVALID
                else:
                    assert len(out) == 1, f"{codon!r} -> {out!r}"
                    value = STOP if out == "*" else ord(out)

                code[idx] = value
                is_stop[idx] = codon in stop_codons
                is_start[idx] = codon in start_codons

    # The kernel decides "is this an in-frame stop?" from `code == STOP`, whereas BioPython
    # asks `codon in table.stop_codons`. Those two must agree, or `cds`'s internal-stop rule
    # would fire in the wrong places.
    for a in ALPHABET:
        for b in ALPHABET:
            for c in ALPHABET:
                idx = codon_index(a + b + c)
                if (code[idx] == STOP) != is_stop[idx]:
                    # Tolerated in one direction only: a dual-coding codon is in stop_codons
                    # but translates as an amino acid.
                    if not (is_stop[idx] and code[idx] not in (STOP, INVALID)):
                        raise AssertionError(
                            f"table {table_id}: codon {a + b + c!r} disagrees: "
                            f"code={code[idx]} is_stop={is_stop[idx]}"
                        )

    unamb = CodonTable.unambiguous_dna_by_id[table_id]
    dual = sorted(c for c in unamb.stop_codons if c in unamb.forward_table)

    return {
        "id": table_id,
        "names": [n for n in unamb.names if n],
        "code": code,
        "is_stop": is_stop,
        "is_start": is_start,
        "dual_coding": dual,
    }


def to_bitset(flags: list[bool]) -> list[int]:
    """Pack N_CODONS bools into N_WORDS u64s."""
    words = [0] * N_WORDS
    for i, f in enumerate(flags):
        if f:
            words[i >> 6] |= 1 << (i & 63)
    return words


def rust_byte_string(code: list[int]) -> str:
    """Render the code array as a compact Rust byte-string literal."""
    out = []
    for v in code:
        ch = chr(v)
        if v < 32 or v >= 127:
            out.append(f"\\x{v:02x}")
        elif ch == "\\":
            out.append("\\\\")
        elif ch == '"':
            out.append('\\"')
        else:
            out.append(ch)
    # Wrap for readability without changing the literal: in Rust a trailing `\` before a newline
    # swallows the newline and the following indentation. Crucially, wrap on *token* boundaries
    # -- splitting mid-escape would turn `\x01` into `\x0` + `1`, which is not valid Rust.
    lines, cur = [], ""
    for tok in out:
        if len(cur) + len(tok) > 96:
            lines.append(cur)
            cur = ""
        cur += tok
    if cur:
        lines.append(cur)
    return 'b"' + "\\\n        ".join(lines) + '"'


def main() -> None:
    ids = sorted(CodonTable.unambiguous_dna_by_id)
    tables = [build_table(i) for i in ids]
    print(f"Verified {len(tables)} tables x {N_CODONS} codons against BioPython {Bio.__version__}")

    # ---------------- Rust ----------------
    r = []
    r.append("// @generated by codegen/generate_tables.py -- DO NOT EDIT BY HAND.")
    r.append(f"// Source of truth: BioPython {Bio.__version__} (NCBI genetic codes).")
    r.append("//")
    r.append("// Every entry was produced by *calling* `Bio.Seq._translate_str` -- the very function")
    r.append("// `Seq.translate()` uses -- on the codon, so the ambiguity rules are BioPython's own")
    r.append("// rather than a reimplementation of them.")
    r.append("")
    r.append("/// The 17 letters that may appear in a codon.")
    r.append("///")
    r.append("/// 'X' is in here but is NOT one of the letters BioPython treats as valid *input*:")
    r.append("/// it expands to GATC (like 'N'), so `CTX` resolves to 'L' and is accepted, but")
    r.append("/// `XXX` does not resolve and is then rejected as an invalid codon.")
    r.append(f'pub const ALPHABET: &[u8; 17] = b"{ALPHABET}";')
    r.append("")
    r.append("/// A true stop codon -> emit the caller's `stop_symbol`.")
    r.append("pub const STOP: u8 = 0;")
    r.append("/// A codon BioPython refuses -> raise \"Codon '...' is invalid\".")
    r.append("pub const INVALID: u8 = 1;")
    r.append("")
    r.append("pub struct CodonTable {")
    r.append("    pub id: u8,")
    r.append("    /// Amino-acid byte per codon index, or `STOP`, or `INVALID`.")
    r.append(f"    pub code: &'static [u8; {N_CODONS}],")
    r.append("    /// Membership of the ambiguity-expanded stop set. Deliberately independent of")
    r.append("    /// `code`: in dual-coding tables (27/28/31) a codon can translate as an amino")
    r.append("    /// acid and still be a legal terminal stop under `cds`.")
    r.append(f"    pub is_stop: &'static [u64; {N_WORDS}],")
    r.append("    /// Membership of the ambiguity-expanded start set (used by `cds`).")
    r.append(f"    pub is_start: &'static [u64; {N_WORDS}],")
    r.append("}")
    r.append("")
    r.append("// The list of dual-coding codons per table lives in `polars_seq/_tables.py`, not")
    r.append("// here: it is only needed to reject `to_stop` and to warn, both of which happen in")
    r.append("// the Python layer once per expression rather than once per row.")
    r.append("")
    r.append("impl CodonTable {")
    r.append("    #[inline(always)]")
    r.append("    pub fn is_stop_codon(&self, idx: usize) -> bool {")
    r.append("        self.is_stop[idx >> 6] & (1u64 << (idx & 63)) != 0")
    r.append("    }")
    r.append("    #[inline(always)]")
    r.append("    pub fn is_start_codon(&self, idx: usize) -> bool {")
    r.append("        self.is_start[idx >> 6] & (1u64 << (idx & 63)) != 0")
    r.append("    }")
    r.append("}")
    r.append("")

    # `#[rustfmt::skip]` on every generated static is load bearing, not cosmetic. `cargo fmt`
    # follows `mod` declarations, so formatting src/lib.rs reaches in here and reflows the
    # bitsets to one entry per line -- after which this file no longer matches what the codegen
    # produces, and the `codegen-is-current` CI job fails on a diff nobody made. Skipping fmt
    # keeps the generated output the single canonical form.
    for t in tables:
        i = t["id"]
        r.append(f"// --- Table {i}: {', '.join(t['names'])} ---")
        r.append("#[rustfmt::skip]")
        r.append(f"static CODE_{i}: &[u8; {N_CODONS}] = {rust_byte_string(t['code'])};")
        r.append("#[rustfmt::skip]")
        r.append(f"static IS_STOP_{i}: &[u64; {N_WORDS}] = &{to_bitset(t['is_stop'])!r};")
        r.append("#[rustfmt::skip]")
        r.append(f"static IS_START_{i}: &[u64; {N_WORDS}] = &{to_bitset(t['is_start'])!r};")
        r.append("")

    r.append("#[rustfmt::skip]")
    r.append(f"pub static TABLES: [CodonTable; {len(tables)}] = [")
    for t in tables:
        i = t["id"]
        r.append(
            f"    CodonTable {{ id: {i}, code: CODE_{i}, is_stop: IS_STOP_{i}, "
            f"is_start: IS_START_{i} }},"
        )
    r.append("];")
    r.append("")
    r.append("/// Look up a table by its NCBI id.")
    r.append("pub fn table_by_id(id: u8) -> Option<&'static CodonTable> {")
    r.append("    TABLES.iter().find(|t| t.id == id)")
    r.append("}")
    r.append("")

    rust_path = ROOT / "src" / "codon_tables.rs"
    rust_path.parent.mkdir(parents=True, exist_ok=True)
    rust_path.write_text("\n".join(r))
    print(f"wrote {rust_path.relative_to(ROOT)} ({rust_path.stat().st_size / 1024:.0f} KB)")

    # ---------------- Python ----------------
    p = []
    p.append('"""NCBI genetic-code table names and ids.')
    p.append("")
    p.append("@generated by codegen/generate_tables.py -- DO NOT EDIT BY HAND.")
    p.append(f'Source of truth: BioPython {Bio.__version__}.')
    p.append('"""')
    p.append("")
    p.append("from __future__ import annotations")
    p.append("")
    p.append(f"TABLE_IDS: frozenset[int] = frozenset({sorted(t['id'] for t in tables)!r})")
    p.append("")
    p.append('"""Every accepted table name/alias, lower-cased, mapped to its NCBI id."""')
    p.append("NAME_TO_ID: dict[str, int] = {")
    for t in tables:
        for n in t["names"]:
            p.append(f"    {n.lower()!r}: {t['id']},")
    p.append("}")
    p.append("")
    p.append('"""Primary display name per id."""')
    p.append("ID_TO_NAME: dict[int, str] = {")
    for t in tables:
        p.append(f"    {t['id']}: {t['names'][0]!r},")
    p.append("}")
    p.append("")
    p.append('"""Tables with codons that are both a stop and an amino acid (NCBI 27/28/31)."""')
    p.append("DUAL_CODING: dict[int, tuple[str, ...]] = {")
    for t in tables:
        if t["dual_coding"]:
            p.append(f"    {t['id']}: {tuple(t['dual_coding'])!r},")
    p.append("}")
    p.append("")

    py_path = ROOT / "python" / "polars_seq" / "_tables.py"
    py_path.parent.mkdir(parents=True, exist_ok=True)
    py_path.write_text("\n".join(p))
    print(f"wrote {py_path.relative_to(ROOT)}")

    n_dual = sum(1 for t in tables if t["dual_coding"])
    print(f"  tables: {ids}")
    print(f"  dual-coding tables: {[t['id'] for t in tables if t['dual_coding']]} ({n_dual})")


if __name__ == "__main__":
    main()
