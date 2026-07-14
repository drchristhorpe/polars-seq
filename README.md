# polars-seq

[![PyPI](https://img.shields.io/pypi/v/polars-seq)](https://pypi.org/project/polars-seq/)
[![Python](https://img.shields.io/pypi/pyversions/polars-seq)](https://pypi.org/project/polars-seq/)
[![CI](https://github.com/drchristhorpe/polars-seq/actions/workflows/ci.yml/badge.svg)](https://github.com/drchristhorpe/polars-seq/actions/workflows/ci.yml)
[![Licence](https://img.shields.io/pypi/l/polars-seq)](LICENSE)

Translate DNA/RNA to protein **inside Polars**, at Rust speed, with **exactly** the semantics of
BioPython's `Seq.translate()`.

```python
import polars as pl
import polars_seq  # noqa: F401 -- importing registers the .seq namespace

df = pl.DataFrame({"dna": ["ATGAAATTTTAA", "ATGGGCCCCTGA"]})

df.with_columns(protein=pl.col("dna").seq.translate())
# ┌──────────────┬─────────┐
# │ dna          ┆ protein │
# ╞══════════════╪═════════╡
# │ ATGAAATTTTAA ┆ MKF*    │
# │ ATGGGCCCCTGA ┆ MGP*    │
# └──────────────┴─────────┘
```

It is a native Polars **expression plugin** (Rust, via `pyo3-polars`), not a Python UDF: it runs
inside the query engine, parallelises across threads, holds no GIL, and composes with the lazy
optimiser and the streaming engine.

---

## Install

```bash
uv add polars-seq
```

or

```bash
pip install polars-seq
```

That's it. Pre-built wheels are published for Linux (x86_64, aarch64), macOS (Intel and Apple
Silicon) and Windows (x64), so **you do not need Rust** — it only gets compiled if you are on a
platform without a wheel, or you are building from source.

Requires Python ≥ 3.10 and polars ≥ 1.0. The extension is built against the stable ABI (`abi3`),
so one wheel serves 3.10 through 3.14 and beyond.

Check it works:

```bash
uv run --with polars-seq --with polars python -c "
import polars as pl, polars_seq
print(pl.DataFrame({'dna': ['ATGAAATTTTAA']}).with_columns(p=pl.col('dna').seq.translate()))
"
```

---

## Building from source (uv)

Only needed if you want to hack on it. You will need a **Rust toolchain**.

### 1. Prerequisites

**uv** — if you don't have it:

```bash
curl -LsSf https://astral.sh/uv/install.sh | sh
```

**Rust** — `rustup` gives you `cargo` and `rustc`. The `stable` channel is pinned in
`rust-toolchain.toml`, so rustup will fetch the right one automatically:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
source "$HOME/.cargo/env"     # or restart your shell
cargo --version               # should print something
```

You also need a C linker (`cc`). On Debian/Ubuntu: `sudo apt install build-essential`.

### 2. Build and install

```bash
git clone https://github.com/drchristhorpe/polars-seq
cd polars-seq

uv sync          # creates .venv on Python 3.14, compiles the Rust extension, installs everything
```

`uv sync` reads `.python-version` (3.14) and `pyproject.toml`, downloads the interpreter if you
don't have it, builds the crate through the `maturin` backend, and installs `polars-seq` into
`.venv` along with the dev dependencies (pytest, biopython).

The first build compiles all of `polars`'s Rust dependencies and takes a few minutes. Later
builds are incremental and take seconds.

Check it works:

```bash
uv run python -c "
import polars as pl, polars_seq
print(pl.DataFrame({'dna': ['ATGAAATTTTAA']}).with_columns(p=pl.col('dna').seq.translate()))
"
```

Run the tests (this is also the BioPython parity check):

```bash
uv run pytest -q
```

### 3. Rebuilding after you change the Rust

**This is the one thing that will bite you.** `uv sync` caches the built wheel and does **not**
notice that you edited `src/*.rs` — you will keep importing the old binary and wonder why your
change did nothing. To actually rebuild, use `maturin develop`, which compiles in place:

```bash
uv run maturin develop --uv            # debug build, fast to compile
uv run maturin develop --uv --release  # optimised; use this for benchmarking
```

Or force uv to rebuild from scratch:

```bash
uv sync --reinstall-package polars-seq --no-cache
```

Pure-Python changes under `python/polars_seq/` need no rebuild at all — the install is editable.

### 4. Regenerating the codon tables

`src/codon_tables.rs` is generated from BioPython and committed. You only need to regenerate it
if you bump BioPython or change the alphabet:

```bash
uv run python codegen/generate_tables.py   # verifies 27 tables x 4913 codons, then writes
uv run maturin develop --uv                # rebuild so the new tables are compiled in
```

### 5. Building a wheel to install elsewhere

```bash
uv build                       # -> dist/polars_seq-0.1.0-*.whl
uv pip install dist/*.whl      # on any machine with Python >= 3.10; no Rust needed
```

### Troubleshooting

**`Both VIRTUAL_ENV and CONDA_PREFIX are set`** — maturin refuses to guess which environment you
mean. If you have conda on your `PATH` (an active `base` env is enough), either `conda deactivate`
first, or just unset it for the one command:

```bash
env -u CONDA_PREFIX uv run maturin develop --uv
```

**`cargo: command not found`** — `~/.cargo/bin` isn't on your `PATH`. `source "$HOME/.cargo/env"`.

**Your Rust change had no effect** — see §3. `uv sync` served you a cached wheel; use
`maturin develop`.

**Python version** — requires ≥ 3.10; developed and tested on **3.14**. The extension is built
against the stable ABI (`abi3`), so one wheel works across versions.

---

## Usage

### `.seq.translate()`

```python
pl.col("dna").seq.translate(
    table=1,           # NCBI id, or name/alias: "Standard", "SGC0", "Vertebrate Mitochondrial"
    stop_symbol="*",
    to_stop=False,     # stop at the first in-frame stop, excluding it
    cds=False,         # validate as a complete coding sequence
    gap="-",
    on_error="raise",  # or "null" -- see below
)
```

Every argument means what it means in BioPython.

```python
df.with_columns(
    protein   = pl.col("dna").seq.translate(),
    orf       = pl.col("dna").seq.translate(to_stop=True),
    mito      = pl.col("dna").seq.translate(table=2),
    bacterial = pl.col("dna").seq.translate(table="Bacterial"),
)
```

`cds=True` validates the sequence as a complete CDS: it must start with a start codon (reported
as `M` whatever it actually encodes), have a length divisible by three, end with a stop codon
(dropped from the output), and contain no internal stop.

```python
pl.col("dna").seq.translate(cds=True)   # "TTGAAATAA" -> "MK"  (TTG is a start codon)
```

Nulls pass through as nulls. A trailing partial codon is dropped, as in BioPython.

### `.seq.reverse_complement()`

IUPAC-aware and case-preserving. Six-frame translation is then just:

```python
df.select(
    **{f"fwd{i}": pl.col("dna").str.slice(i).seq.translate() for i in range(3)},
    **{f"rev{i}": pl.col("dna").seq.reverse_complement().str.slice(i).seq.translate()
       for i in range(3)},
)
```

### `polars_seq.codon_tables()`

All 27 NCBI genetic codes as a DataFrame — ids, names, aliases, and which are dual-coding.

---

## Ambiguity codes are handled properly

Ambiguous IUPAC nucleotides resolve the way BioPython resolves them, which is more subtle than
"anything unclear becomes `X`":

| codon | → | why |
|---|---|---|
| `GGN` | `G` | all four `GGx` codons are Gly, so there is no ambiguity to report |
| `TAR` | `*` | `R`=A/G, and both `TAA` and `TAG` are stops |
| `TAN` | `X` | expands to two stops *and* two Tyr — genuinely unresolvable |
| `RAT` | **`B`** | `GAT`=Asp, `AAT`=Asn → **`B`** is the IUPAC code for "Asp or Asn" |
| `SAA` | **`Z`** | `CAA`=Gln, `GAA`=Glu → **`Z`** = "Glu or Gln" |

`RAT` → `B` is the one that catches people out. An implementation that emits `X` for every
unresolvable codon looks correct until someone runs it on ambiguous data.

---

## The X-codon quirk

While building this I ran into a genuinely surprising corner of `Seq.translate()`. It is worth
knowing about whether or not you use this library.

**`CTX` translates to `L`. `XXX` is an error.**

Both contain `X`. Here is why they differ.

BioPython keeps *two* different sets of nucleotide letters, and they disagree with each other:

1. the **expansion table**, which says what each letter can stand for. It has **17** keys, and
   `X` is one of them — it expands to `GATC`, exactly like `N`.
2. the **validity check**, which decides whether a codon it *failed* to translate is at least
   made of legal characters. That set is `_ambiguous_dna_letters | _ambiguous_rna_letters` —
   **16** letters, and `X` is **not** among them.

Translation consults the expansion table **first**, and falls back to the validity check only
when the expansion fails. And the expansion fails in exactly one situation: **when it runs into a
stop codon.** So:

- `CTX` → expands to `CTA/CTC/CTG/CTT` → all four are Leucine, no stops → resolves to **`L`**.
  The validity check is never reached and the `X` costs nothing.
- `XXX` → expands to all 64 codons → amino acids *and stops* → the expansion gives up → *now*
  the validity check runs, sees the `X`, and rejects the codon: `Codon 'XXX' is invalid`.

The precise rule, which we verified exhaustively against BioPython (817 X-bearing codons × 27
tables, zero counterexamples):

> **An `X`-bearing codon is accepted if and only if none of the codons it expands to is a stop.**
> When accepted, it means exactly what the `N` spelling means.

So `X` is a perfect synonym for `N` — right up to the moment the expansion touches a stop codon,
where `N` degrades gracefully and `X` becomes a hard error:

| expansion | `N` spelling | `X` spelling |
|---|---|---|
| one amino acid, no stops | `CTN` → `L` | `CTX` → `L` |
| several amino acids, no stops | `AAN` → `X` | `AAX` → `X` |
| amino acids **and stops** | `TAN` → `X` | `TAX` → **error** |
| everything | `NNN` → `X` | `XXX` → **error** |

Note the second row: an accepted `X` codon *can* come out as `X` — there `X` is an **amino-acid**
ambiguity code (Lys-or-Asn), not a nucleotide. What it can never be is a stop.

This is easy to get wrong in both directions: treat `X` as invalid everywhere and you break
`CTX` → `L`; treat it as a plain synonym for `N` and you wrongly accept `TAX` and `XXX`. I got it
wrong the first way, and only the differential fuzz caught it — on the sequence `CTTCTX`. Every
hand-written test I had passed. `polars-seq` now reproduces the real behaviour and re-sweeps all
817 X-bearing codons × 27 tables against BioPython on each test run.

*This is BioPython's actual, current behaviour — 1.78 and 1.87 agree — so it is a quirk to match,
not a bug to route around. It falls out of the two letter-sets having drifted apart.*

---

## Differences from BioPython

Two, both deliberate, because a DataFrame is not a single sequence.

**1. `on_error="null"`.** BioPython raises on a malformed sequence, and so do we, by default. But
in a million-row frame, one bad sequence aborting the whole query is usually not what you want:

```python
df.with_columns(protein=pl.col("dna").seq.translate(on_error="null"))
# invalid sequences become null; everything else still translates
```

The default, `on_error="raise"`, reports the row index and the offending codon.

**2. No per-row warnings.** BioPython emits a `BiopythonWarning` for a trailing partial codon.
Warning once per row from inside a parallel Rust kernel is not viable, so we are quiet about it.
The *behaviour* is identical — the partial codon is dropped either way.

Dual-coding tables (27/28/31) do still warn, once, when the expression is built — and `to_stop`
with those tables is still rejected, exactly as BioPython rejects it.

---

## Correctness

The parity claim is enforced, not asserted. `uv run pytest` runs:

- **an exhaustive codon sweep** — all 4913 codons × all 27 NCBI tables = **132,651 comparisons**
  against BioPython, including which codons it *refuses*;
- **differential fuzz** — thousands of random sequences over ten alphabets (unambiguous, IUPAC,
  RNA, gapped, lower-case, invalid-character, …) crossed with the full option space, asserting we
  produce the same protein *and* fail on exactly the same inputs;
- **golden tests** for every documented rule and trap;
- **Polars integration** — nulls, empty frames, chunked and sliced Series, lazy, streaming,
  `group_by().agg()`.

The lookup tables are **generated from BioPython** (`codegen/generate_tables.py`) by calling the
real `_translate_str` on every codon, rather than by re-implementing its resolver in Rust — which
is exactly the sort of code that produces `RAT` → `X` bugs. BioPython is a build- and test-time
dependency only; it is not needed at runtime.

---

## Performance

A Rust kernel with one array lookup per codon, no per-row allocation, parallel across chunks.

`uv run python tools/validate.py` reproduces the numbers on your own machine and data, and writes
a full row-by-row BioPython comparison to `tmp/` so you can inspect the output rather than trust
it.

---

## Layout

```
codegen/generate_tables.py   BioPython -> Rust lookup tables (self-verifying)
src/translate.rs             the kernel: framing, stops, gaps, cds rules
src/codon_tables.rs          GENERATED -- do not edit
src/expressions.rs           Polars expression entry points
python/polars_seq/           the .seq namespace and argument validation
tests/                       golden, exhaustive, differential-fuzz, integration
tools/validate.py            writes validation + benchmark artefacts to tmp/
```

## Citing — please cite BioPython

**If you find this tool useful, cite BioPython.** This library would not exist without it.

BioPython is not merely an inspiration here, it is the *specification*: the codon tables shipped
in this package were generated by calling BioPython's own resolver on every codon, its semantics
are what the test suite asserts against, and every subtlety `polars-seq` gets right — `RAT` → `B`,
the dual-coding tables, the X-codon quirk — is right because BioPython worked it out first. All
`polars-seq` adds is speed and a Polars namespace.

> Cock, P.J.A., Antao, T., Chang, J.T., Chapman, B.A., Cox, C.J., Dalke, A., Friedberg, I.,
> Hamelryck, T., Kauff, F., Wilczynski, B. and de Hoon, M.J.L. (2009)
> **Biopython: freely available Python tools for computational molecular biology and
> bioinformatics.** *Bioinformatics* **25**(11), 1422–1423.
> <https://doi.org/10.1093/bioinformatics/btp163> · PMID: 19304878

```bibtex
@article{Cock2009Biopython,
  author  = {Cock, Peter J. A. and Antao, Tiago and Chang, Jeffrey T. and Chapman, Brad A.
             and Cox, Cymon J. and Dalke, Andrew and Friedberg, Iddo and Hamelryck, Thomas
             and Kauff, Frank and Wilczynski, Bartek and de Hoon, Michiel J. L.},
  title   = {Biopython: freely available {P}ython tools for computational molecular biology
             and bioinformatics},
  journal = {Bioinformatics},
  year    = {2009},
  volume  = {25},
  number  = {11},
  pages   = {1422--1423},
  doi     = {10.1093/bioinformatics/btp163},
  pmid    = {19304878}
}
```

## Licence

MIT. BioPython is distributed under the [Biopython License Agreement / BSD 3-Clause
](https://github.com/biopython/biopython/blob/master/LICENSE.rst); it is used here at build and
test time only and is not redistributed as part of this package.
