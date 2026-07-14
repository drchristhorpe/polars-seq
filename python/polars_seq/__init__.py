"""polars-seq -- translate DNA/RNA to protein inside Polars.

Importing this module registers a ``.seq`` namespace on Polars expressions::

    import polars as pl
    import polars_seq  # noqa: F401  (import for the side effect)

    df.with_columns(protein=pl.col("dna").seq.translate())

The semantics are those of BioPython's ``Seq.translate()``, reproduced exactly -- including
its handling of ambiguous IUPAC codons, gaps, partial codons and ``cds`` validation. See the
README for the (two, deliberate) differences.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import polars as pl
from polars.plugins import register_plugin_function

from ._tables import DUAL_CODING, ID_TO_NAME, NAME_TO_ID, TABLE_IDS

if TYPE_CHECKING:
    from polars._typing import IntoExpr

__all__ = ["SeqNameSpace", "codon_tables", "__version__"]
__version__ = "0.1.1"

_PLUGIN_PATH = Path(__file__).parent


class AmbiguousStopCodonWarning(UserWarning):
    """Raised for NCBI tables 27/28/31, where a codon is both a stop and an amino acid.

    Mirrors the ``BiopythonWarning`` that ``Seq.translate()`` emits for these tables.
    """


def codon_tables() -> pl.DataFrame:
    """Return the NCBI genetic-code tables this plugin supports, as a DataFrame."""
    return pl.DataFrame(
        {
            "id": sorted(TABLE_IDS),
            "name": [ID_TO_NAME[i] for i in sorted(TABLE_IDS)],
            "aliases": [
                sorted(n for n, tid in NAME_TO_ID.items() if tid == i and n != ID_TO_NAME[i].lower())
                for i in sorted(TABLE_IDS)
            ],
            "dual_coding": [list(DUAL_CODING.get(i, ())) for i in sorted(TABLE_IDS)],
        }
    )


def _resolve_table(table: int | str) -> int:
    """Accept an NCBI id (1), a stringified id ("1"), or a name/alias ("Standard", "SGC0")."""
    # BioPython tries int() first, so table="1" is the id, not a name.
    try:
        table_id = int(table)
    except (TypeError, ValueError):
        pass
    else:
        if table_id not in TABLE_IDS:
            raise ValueError(
                f"unknown codon table id {table_id!r}; valid ids are {sorted(TABLE_IDS)}"
            )
        return table_id

    if not isinstance(table, str):
        raise TypeError(f"table must be an int id or a str name, got {type(table).__name__}")

    try:
        return NAME_TO_ID[table.strip().lower()]
    except KeyError:
        raise ValueError(
            f"unknown codon table name {table!r}; "
            f"call polars_seq.codon_tables() to list the {len(TABLE_IDS)} available tables"
        ) from None


def _one_ascii_char(value: str, argname: str) -> str:
    if not isinstance(value, str):
        raise TypeError(f"{argname} must be a single-character string, got {type(value).__name__}")
    if len(value) != 1:
        raise ValueError(f"{argname} must be a single character, got {value!r}")
    if not value.isascii():
        raise ValueError(f"{argname} must be ASCII, got {value!r}")
    return value


@pl.api.register_expr_namespace("seq")
class SeqNameSpace:
    """The ``.seq`` expression namespace."""

    def __init__(self, expr: pl.Expr) -> None:
        self._expr = expr

    def translate(
        self,
        table: int | str = 1,
        stop_symbol: str = "*",
        to_stop: bool = False,
        cds: bool = False,
        gap: str | None = "-",
        *,
        on_error: Literal["raise", "null"] = "raise",
    ) -> pl.Expr:
        """Translate a nucleotide column to a protein column.

        Arguments match ``Bio.Seq.Seq.translate``:

        table
            NCBI genetic code, as an id (``1``) or a name/alias (``"Standard"``, ``"SGC0"``,
            ``"Vertebrate Mitochondrial"``). Case-insensitive.
        stop_symbol
            Character emitted for an in-frame stop codon.
        to_stop
            Stop at the first in-frame stop codon and exclude it from the output.
        cds
            Validate the sequence as a complete coding sequence: it must begin with a start
            codon (reported as ``M`` whatever it encodes), have a length that is a multiple of
            three, end with a stop codon (excluded from the output), and contain no internal
            stop. Sequences that fail are errors.
        gap
            Gap character. A whole codon of gaps translates to one gap character; a partial
            one (``"A--"``) is an error. Pass ``None`` to disable.
        on_error
            ``"raise"`` (default, BioPython's behaviour) aborts the query on the first invalid
            sequence. ``"null"`` maps invalid sequences to ``null`` and translates the rest --
            usually what you want on real, messy data.

        Null inputs produce null outputs.
        """
        table_id = _resolve_table(table)
        stop_symbol = _one_ascii_char(stop_symbol, "stop_symbol")
        if gap is not None:
            gap = _one_ascii_char(gap, "gap")

        if on_error not in ("raise", "null"):
            raise ValueError(f"on_error must be 'raise' or 'null', got {on_error!r}")

        # Tables 27/28/31 have codons that are both a stop and an amino acid. BioPython refuses
        # to combine those with to_stop (there is no single right answer) and warns otherwise;
        # we do the same, once per expression rather than once per row.
        if table_id in DUAL_CODING:
            dual = DUAL_CODING[table_id]
            if to_stop:
                raise ValueError(
                    f"cannot use to_stop=True with table {table_id} "
                    f"({ID_TO_NAME[table_id]}): {', '.join(dual)} can be both a stop codon "
                    f"and an amino acid"
                )
            warnings.warn(
                f"table {table_id} ({ID_TO_NAME[table_id]}) contains {len(dual)} codon(s) "
                f"({', '.join(dual)}) coding for both STOP and an amino acid; "
                f"they will be translated as the amino acid",
                AmbiguousStopCodonWarning,
                stacklevel=2,
            )

        return register_plugin_function(
            plugin_path=_PLUGIN_PATH,
            function_name="translate_expr",
            args=self._expr,
            kwargs={
                "table": table_id,
                "stop_symbol": stop_symbol,
                "to_stop": to_stop,
                "cds": cds,
                "gap": gap,
                "null_on_error": on_error == "null",
            },
            is_elementwise=True,
        )

    def reverse_complement(self) -> pl.Expr:
        """Reverse-complement a nucleotide column, IUPAC-aware and case-preserving.

        Combine with :meth:`translate` for the reverse strand::

            df.with_columns(rev_protein=pl.col("dna").seq.reverse_complement().seq.translate())
        """
        return register_plugin_function(
            plugin_path=_PLUGIN_PATH,
            function_name="reverse_complement_expr",
            args=self._expr,
            is_elementwise=True,
        )


def translate(expr: IntoExpr, /, **kwargs) -> pl.Expr:
    """Function form of :meth:`SeqNameSpace.translate`, for when you prefer it."""
    return pl.col(expr).seq.translate(**kwargs) if isinstance(expr, str) else expr.seq.translate(**kwargs)
