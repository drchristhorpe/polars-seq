# Changelog

All notable changes to `polars-seq` are documented here.

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [0.1.0] — 2026-07-14

First release. A native Polars expression plugin that translates DNA/RNA to protein with the
exact semantics of BioPython's `Seq.translate()`.

### Added

- **`pl.col(...).seq.translate()`** — a Polars expression plugin written in Rust
  (`pyo3-polars`), so translation runs inside the query engine: parallel across chunks, no GIL,
  no per-row Python objects, and composable with the lazy optimiser and streaming engine.
  Supports the full BioPython argument set — `table`, `stop_symbol`, `to_stop`, `cds`, `gap` —
  with identical behaviour.
- **All 27 NCBI genetic codes**, selectable by id (`table=2`), by name
  (`table="Vertebrate Mitochondrial"`), or by alias (`table="SGC1"`), case-insensitively.
- **Full IUPAC ambiguity support**, resolved the way BioPython resolves it: `GGN` → `G`,
  `TAR` → `*`, `TAN` → `X`, and — the one people get wrong — `RAT` → `B` (Asx) and `SAA` → `Z`
  (Glx), rather than collapsing every unresolvable codon to `X`.
- **`cds=True` validation**: start-codon check (with the start reported as `M` whatever it
  encodes), length-divisible-by-three check, terminal-stop check (excluded from the output), and
  internal-stop rejection — applied in BioPython's order, which is observable through the error
  message you get when a sequence violates more than one rule.
- **Dual-coding tables (27, 28, 31)** handled correctly: codons that are simultaneously a stop
  and an amino acid translate as the amino acid, yet remain valid CDS terminators. `to_stop` is
  rejected for these tables, and a warning is emitted once per expression — matching BioPython.
- **`pl.col(...).seq.reverse_complement()`** — IUPAC-aware, case-preserving; makes six-frame
  translation a one-liner.
- **`polars_seq.codon_tables()`** — the supported genetic codes as a DataFrame.
- **`on_error="null"`** — map invalid sequences to `null` instead of aborting the query. The
  default (`"raise"`) matches BioPython, and reports the row index and offending codon.

### Correctness

- Codon lookup tables are **generated from BioPython** (`codegen/generate_tables.py`) by calling
  the real `_translate_str` on every codon, rather than by re-implementing its ambiguity
  resolver in Rust. BioPython is a build- and test-time dependency only, never a runtime one.
- **Exhaustive verification**: all 4913 codons × 27 tables (132,651 comparisons) are checked
  against BioPython at codegen time *and again* in the test suite, including which codons
  BioPython *refuses*.
- **Differential fuzz** against BioPython across ten alphabets and the full option space,
  asserting both that we produce the same protein and that we fail on exactly the same inputs.

### Notes

#### The X-codon quirk

Discovered while fuzzing, and worth recording because it is genuinely counter-intuitive:

> **`CTX` translates to `L`, but `XXX` is an error.**

BioPython holds two sets of nucleotide letters that disagree with each other. The *expansion*
table has 17 keys including `X` (which expands to `GATC`, exactly like `N`); the *validity* check
has only 16, and `X` is not one of them. Expansion is tried first and only fails when it runs
into a stop codon — at which point the validity check runs, spots the `X`, and rejects the codon.

The exact rule, verified exhaustively against BioPython (817 X-bearing codons × 27 tables, zero
counterexamples):

> An `X`-bearing codon is accepted **iff none of the codons it expands to is a stop**. When
> accepted, it means exactly what the `N` spelling means.

So `X` is a perfect synonym for `N` until the expansion touches a stop, where `N` degrades to `X`
and `X` becomes a hard error: `CTN`/`CTX` → `L`, `AAN`/`AAX` → `X`, but `TAN` → `X` while
`TAX` → error, and `NNN` → `X` while `XXX` → error. (An accepted `X` codon *can* come out as `X`
— as an **amino-acid** ambiguity code, e.g. `AAX` = Lys-or-Asn. It can never come out as a stop.)

The first implementation here treated `X` as invalid everywhere. That passed every hand-written
test and was caught only by the differential fuzz, on the sequence `CTTCTX`. The alphabet was
widened from 16 to 17 letters, a distinct `INVALID` outcome was added to the generated tables,
and the X-codon sweep is now a permanent test. This is BioPython's real, current behaviour (1.78
and 1.87 agree), so it is matched, not "fixed". See the README section of the same name.

#### Deliberate differences from BioPython

- `on_error="null"` has no BioPython equivalent; it exists because one malformed row should not
  have to abort a million-row query. The default remains BioPython's raise-on-error.
- No per-row warnings. BioPython emits a `BiopythonWarning` for a trailing partial codon; doing
  that per row from a parallel Rust kernel is not viable, so we are quiet. The behaviour is
  unchanged — the partial codon is dropped either way.

[0.1.0]: https://github.com/christhorpe/polars-seq/releases/tag/v0.1.0
