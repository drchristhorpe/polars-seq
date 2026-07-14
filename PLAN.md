# PLAN — `polars-seq`: a Polars plugin for DNA → protein translation

A native (Rust) Polars expression plugin that translates nucleotide sequences to amino-acid
sequences, **bit-for-bit compatible with BioPython's `Seq.translate()`**.

---

## 1. Goal & success criterion

Expose a Polars expression namespace so that this works on a `DataFrame`/`LazyFrame` column:

```python
import polars as pl
import polars_seq  # registers the `.seq` namespace

df.with_columns(protein=pl.col("dna").seq.translate())
```

**Success criterion (the whole point):** for any input sequence and any parameter combination,

```python
df.select(pl.col("dna").seq.translate(**kw)).to_series().to_list()
==
[str(Bio.Seq.Seq(s).translate(**kw)) for s in dna]
```

This is not an aspiration — it is enforced by a randomised differential test against BioPython
across all 27 NCBI tables (§7). Anything that does not match BioPython is a bug.

Non-goal: inventing a "better" translation semantics. Where BioPython is quirky, we are quirky.

---

## 2. Why a Rust plugin (and not a Python UDF)

The obvious implementation is `pl.col("dna").map_elements(lambda s: str(Seq(s).translate()))`.
That is correct but pays Python-interpreter overhead *per row* and holds the GIL, so it does not
parallelise and it defeats the point of using Polars at all.

A Polars **expression plugin** (`pyo3-polars`) compiles to a shared library that Polars calls
directly on the Arrow buffers. It runs inside the engine: parallel over chunks, no GIL, no
per-row Python object churn, and it composes with the lazy optimiser. That is the canonical
meaning of "Polars plugin", and it is what we build.

Expected shape of the win: a tight `&str → String` loop over contiguous UTF-8 buffers, with an
O(1) array lookup per codon. Benchmark vs. BioPython is part of the deliverable (§9), so the
claim gets measured rather than asserted.

---

## 3. Toolchain (verified before planning, not assumed)

| Component | Version | Status |
|---|---|---|
| Python | **3.14.3** | uv-managed; user's stated target |
| polars (py) | 1.42.1 | verified importable on 3.14.3 |
| polars (rust) | 0.54.4 | pinned to match |
| pyo3-polars | 0.27 | depends on polars 0.54.4 + pyo3 0.28 |
| maturin | 1.14 | build backend |
| Rust | stable (installed via rustup) | `cc` present for linking |
| BioPython | 1.78 | **test oracle + codegen source only** — not a runtime dependency |

Rust was not installed at the start; it has been installed. Python 3.14 + polars was verified by
actually resolving and importing it, not by reading a compatibility table.

---

## 4. The BioPython semantics we must reproduce

Derived by reading `Bio.Seq._translate_str`, `Bio.Data.CodonTable.AmbiguousForwardTable`, and
`list_possible_proteins` in the installed 1.78 source, then confirming each rule empirically.
**This section is the specification.** Every claim below was executed against BioPython.

### 4.1 Signature

`translate(table="Standard", stop_symbol="*", to_stop=False, cds=False, gap="-")`

### 4.2 The ambiguity resolver (the subtle part)

`Seq.translate()` uses `CodonTable.ambiguous_generic_by_id[...]`, whose forward table resolves
*ambiguous* codons. This is where a naive implementation goes wrong. For a codon:

1. Expand each IUPAC nucleotide letter to its concrete bases and form the cartesian product.
2. Look each concrete codon up in the table's forward table.
   - **all concrete codons are stops** → the codon is a **stop** (`TAR` → `*`).
   - **mixed stops and amino acids** → not resolvable → **`X`** (`TAN` → `X`).
   - **exactly one distinct amino acid** → that amino acid (`GGN` → `G`).
   - **several distinct amino acids** → find every extended-IUPAC *protein* letter whose
     expansion is a superset of the possible set, and pick the least-ambiguous one
     (ties broken alphabetically).

Step 2's last rule is the trap: **`RAT` → `B`, not `X`.** `GAT`=Asp, `AAT`=Asn, and `B` is the
IUPAC code for "Asp or Asn". Likewise `SAA` → `Z` (Glu or Gln). An implementation that emits `X`
for every unresolvable codon passes casual testing and is silently wrong here.

The mapping is `Bio.Data.IUPACData.extended_protein_values` (`B: ND`, `Z: QE`, `J: IL`,
`X`: the 20 standard residues, plus identity entries).

### 4.3 Alphabet, and the `X` trap

> **REVISED DURING EXECUTION.** The rule below is the corrected one. My first pass got this
> wrong, shipped a 16-letter alphabet, and the differential fuzz (§7.2) caught it on `CTTCTX`.
> Recording the correction here because it is the single subtlest thing in the whole project.

The generic table works internally in **RNA space**: its nucleotide-expansion dict maps `T → U`,
and has **17** keys — including `X`, which expands to `GATC` exactly like `N`. But
`_translate_str`'s *validity* check uses `_ambiguous_dna_letters + _ambiguous_rna_letters` —
only **16 letters, and `X` is not among them**:

```
A C G T U R Y W S M K H B V D N          (valid input letters)
A C G T U R Y W S M K H B V D N X        (letters the expansion table can resolve)
```

Translation tries the expansion **first**, and reaches the validity check only when the expansion
fails — which happens exactly when it runs into a **stop codon**. Hence the rule, verified
exhaustively (817 X-bearing codons × 27 tables, zero counterexamples):

> **An `X`-bearing codon is accepted iff none of the codons it expands to is a stop.** When
> accepted, it means precisely what the `N` spelling means.

| expansion | `N` spelling | `X` spelling |
|---|---|---|
| one amino acid, no stops | `CTN` → `L` | `CTX` → **`L`** |
| several amino acids, no stops | `AAN` → `X` | `AAX` → **`X`** |
| amino acids **and** stops | `TAN` → `X` | `TAX` → **error** |
| everything | `NNN` → `X` | `XXX` → **error** |

Note `AAX` → `X`: an accepted X-codon *can* output `X`, but there `X` is an **amino-acid**
ambiguity code (Lys-or-Asn), not a nucleotide. It can never output a stop.

So the lookup is built over **17 letters (17³ = 4913 codons)**, and each entry is one of: an
amino-acid byte, `STOP`, or `INVALID` (the codons BioPython refuses). Anything containing a
character outside those 17 is invalid input, handled separately (and is where `gap` is checked).

Input is uppercased first, so `atgaaa` → `MK`. `U` is accepted, so RNA works: `AUGAAA` → `MK`.

### 4.4 Gaps

A codon that is exactly `gap*3` translates to a single `gap` character. The gap check happens
*after* the validity check, so a **partial** gap codon is an error:

- `------` → `--`
- `A-A` → `TranslationError: Codon 'A-A' is invalid`

### 4.5 Partial codons

Trailing bases that do not complete a codon are **silently dropped** (BioPython warns, then
translates `floor(n/3)` codons anyway):

- `ATGAA` → `M`, `ATGA` → `M`, `AT` → `''`, `''` → `''`

### 4.6 `to_stop` and `stop_symbol`

- `to_stop=True` breaks at the first in-frame stop and **excludes** it: `ATGAAATAAGGG` → `MK`
- `stop_symbol="#"`: `ATGAAATAAGGG` → `MK#G`

### 4.7 `cds=True`

Validated as a complete coding sequence, in this order:
1. first codon ∈ (ambiguity-expanded) **start** codons, else `First codon 'AAA' is not a start codon`
2. `len % 3 == 0`, else `Sequence length N is not a multiple of three`
3. last codon ∈ **stop** codons, else `Final codon 'ATG' is not a stop codon`
4. the start codon is force-translated to **`M`** regardless of what it encodes (`TTGAAATAA` → `MK`)
5. the terminal stop is **excluded** from the output
6. any *internal* in-frame stop → `Extra in frame stop codon found.`
7. `to_stop` is ignored under `cds=True`

Verified: `ATGCCCTAG` → `MP`, `ATGTAA` → `M`, `ATG` → error (no stop).

### 4.8 Dual-coding tables (27, 28, 31)

Tables `27` (Karyorelict), `28` (Condylostoma) and `31` (Blastocrithidia) contain codons that are
**both** a stop and an amino acid (e.g. table 28: `TAA`, `TAG`, `TGA`). BioPython:

- looks the forward table up **first**, so such codons translate as the **amino acid**;
- still keeps them in `stop_codons`, so they remain valid *terminal* codons under `cds=True`;
- **raises `ValueError` if `to_stop=True`** on such a table;
- emits a `BiopythonWarning` otherwise.

So "is an amino acid" and "is a stop codon" must be stored as **independent** facts per codon, not
as one enum. This falls out of the design in §5.

### 4.9 Tables

27 NCBI tables: ids `1..6, 9..16, 21..31, 33` (exactly: 1,2,3,4,5,6,9,10,11,12,13,14,15,16,21,22,
23,24,25,26,27,28,29,30,31,32,33). Selectable by **id** (`table=2`) or by **name**
(`table="Vertebrate Mitochondrial"`, incl. aliases like `SGC1`). Verified both give `M**` for
`ATGAGATAA` under table 2 (`AGA` is a stop in vertebrate mitochondria).

---

## 5. Design: precompute the resolver, don't reimplement it

The single most important design decision.

Re-implementing §4.2's resolver in Rust by hand is exactly the kind of fiddly logic that produces
the `RAT → X` class of bug. Instead:

> **Use BioPython itself as the oracle at code-generation time.** For each of the 27 tables,
> enumerate all 4096 codons over the 16-letter alphabet, ask BioPython what each one does, and
> emit the answers as a static Rust table.

`codegen/generate_tables.py` produces `src/codon_tables.rs`, which is **committed** to the repo.
BioPython is therefore a *build-time and test-time* dependency only — the shipped wheel has zero
Python dependencies beyond polars, and the runtime does no ambiguity reasoning at all.

Per table we emit three parallel 4096-entry arrays, indexed by
`idx = c1*256 + c2*16 + c3` (each letter mapped to 0..15):

| array | meaning |
|---|---|
| `code: [u8; 4096]` | amino-acid byte (`b'K'`…), or `0` = stop-only, or `1` = unresolvable (→ `pos_stop`, `X`) |
| `is_stop: [u64; 64]` (bitset) | membership of the ambiguity-expanded stop-codon set |
| `is_start: [u64; 64]` (bitset) | membership of the ambiguity-expanded start-codon set |

`code` and `is_stop` are deliberately independent — that is precisely what §4.8 requires.

The codegen script also **asserts** its own output against BioPython before writing, so a bad
generation fails loudly rather than silently shipping.

Translation then reduces to: uppercase byte → 4-bit code → array index → one `u8`. No hashing, no
allocation per codon, branch-predictable.

---

## 6. API

```python
pl.col("dna").seq.translate(
    table=1,             # int id or str name ("Standard", "SGC1", ...)
    stop_symbol="*",
    to_stop=False,
    cds=False,
    gap="-",
    on_error="raise",    # extension: "raise" | "null"
)
```

Faithful to BioPython, with two deliberate, documented deviations that exist because a
*dataframe* is not a *single sequence*:

1. **`on_error="null"`** — in a 10-million-row frame, one malformed sequence killing the entire
   query is usually the wrong behaviour. `on_error="null"` maps the offending row to `null` and
   lets the rest through. Default remains `"raise"` (BioPython's behaviour); the error message
   quotes the row index and the offending codon.
2. **No per-row warnings.** BioPython warns on partial codons; emitting a Python warning per row
   from inside a parallel Rust kernel is not viable. Behaviour is identical (§4.5), just quiet.
   Documented in the README.

`null` in → `null` out.

Also shipped, because it is ~20 lines and makes six-frame translation possible (the immediate
next thing anyone asks for): `pl.col("dna").seq.reverse_complement()`. Everything else
(`gc_content`, ORF finding, …) is explicitly out of scope for v0.1.0.

Validation of `table` / `stop_symbol` / `gap` happens in the **Python** layer, once per expression,
so bad arguments raise a clean `ValueError` immediately instead of a `ComputeError` from inside
the engine.

---

## 7. Testing strategy

The parity claim in §1 is only worth as much as the test that enforces it.

1. **Golden unit tests** — every example in §4, asserted explicitly. These encode the traps
   (`RAT`→`B`, `XXX`→invalid, `A-A`→invalid, `TTGAAATAA`+cds→`MK`) so a future refactor cannot
   quietly regress them.
2. **Differential fuzz vs. BioPython — the load-bearing test.** Generate thousands of random
   sequences (varying length incl. non-multiples of 3, over plain ACGT / IUPAC-ambiguous /
   gapped / RNA / lowercase / invalid-character alphabets), cross-producted with all 27 tables and
   the `to_stop` / `cds` / `stop_symbol` / `gap` options. Assert our output equals BioPython's —
   **and that our errors occur exactly where BioPython's do.** Seeded, so failures reproduce.
3. **Exhaustive codon sweep** — all 4096 codons × all 27 tables against BioPython. This is cheap
   (110k comparisons) and total: it makes the lookup tables provably correct rather than sampled.
4. **Polars integration** — nulls, empty frames, chunked/sliced series, `LazyFrame` + `.collect()`,
   `group_by().agg()`, and multi-threaded execution over a large frame.
5. **Benchmark** — vs. BioPython `map_elements`, to substantiate §2.

## 8. Layout

```
polars_seq/
├── PLAN.md  CHANGELOG.md  README.md
├── Cargo.toml            # rust: polars 0.54, pyo3-polars 0.27
├── pyproject.toml        # maturin backend, requires-python >=3.10, target 3.14
├── rust-toolchain.toml
├── codegen/generate_tables.py     # BioPython -> Rust; self-asserting
├── src/
│   ├── lib.rs
│   ├── codon_tables.rs   # GENERATED, committed
│   ├── translate.rs      # the kernel
│   └── expressions.rs    # #[polars_expr] entry points
├── python/polars_seq/
│   ├── __init__.py       # registers the `.seq` namespace
│   ├── _tables.py        # name/alias -> id resolution, arg validation
│   └── py.typed
└── tests/
    ├── test_golden.py  test_parity.py  test_codons.py  test_polars.py
    └── bench.py
```

## 9. Deliverables

- Working, `maturin`-built plugin importable on Python 3.14.
- `CHANGELOG.md`, `README.md` (usage-focused).
- **Validation artefacts written to a tmp directory for review** (final step): the exhaustive
  4096×27 codon sweep result, the fuzz-parity report, a side-by-side
  BioPython-vs-`polars_seq` CSV, and the benchmark numbers — so the parity claim can be
  inspected rather than taken on trust.

## 10. Risks

| Risk | Mitigation |
|---|---|
| Ambiguity resolver subtly wrong | Not hand-written — generated from BioPython and exhaustively verified (§5, §7.3) |
| pyo3-polars ↔ py-polars ABI mismatch | Versions pinned to a known-compatible set (§3); integration tests would fail immediately |
| Dual-coding tables mishandled | Independent `code`/`is_stop` arrays (§5); covered by the 27-table sweep |
| BioPython 1.78 (2020) tables stale vs. current NCBI | Tables 1–33 are stable; codegen is re-runnable against any newer BioPython |
| Silent divergence appearing later | Parity fuzz is a permanent test, not a one-off script |
```
