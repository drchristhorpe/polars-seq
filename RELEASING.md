# Releasing `polars-seq` to PyPI

Publishing is automated: **push a `v*` tag and GitHub Actions builds every wheel, smoke-tests
one, and publishes to PyPI.** There is no API token to create, store, or rotate — it uses PyPI's
Trusted Publishing (OIDC).

The one-time setup below only has to be done once, before the first release.

---

## One-time setup (do this before v0.1.0)

### 1. Register the Trusted Publisher on PyPI

The project does not exist on PyPI yet, so create a **pending** publisher — PyPI will create the
project automatically on first upload.

1. Go to <https://pypi.org/manage/account/publishing/>
2. Under **"Add a new pending publisher"**, fill in exactly:

   | field | value |
   |---|---|
   | PyPI Project Name | `polars-seq` |
   | Owner | `drchristhorpe` |
   | Repository name | `polars-seq` |
   | Workflow name | `release.yml` |
   | Environment name | `pypi` |

3. Save.

> The **Environment name must be `pypi`** — it has to match `environment: name: pypi` in
> `.github/workflows/release.yml`, or PyPI will reject the upload.

The name `polars-seq` was available at the time of writing. If someone has taken it since, change
`name` in `pyproject.toml` and the values above together.

### 2. Create the `pypi` environment on GitHub

1. Go to <https://github.com/drchristhorpe/polars-seq/settings/environments>
2. **New environment** → name it `pypi` → Configure.
3. Optionally add yourself as a **required reviewer**, so every publish waits for you to click
   approve. Recommended: it turns an accidental tag push into a prompt rather than a release.

Nothing else is needed. Do **not** add a `PYPI_API_TOKEN` secret — trusted publishing does not
use one, and adding one would just be a credential sitting there to leak.

---

## Cutting a release

1. **Bump the version in `pyproject.toml`** (and `__version__` in
   `python/polars_seq/__init__.py`, and `Cargo.toml` — all three should agree).
2. **Update `CHANGELOG.md`**: move the changes under a new version heading with today's date.
3. Commit, then tag and push:

   ```bash
   git commit -am "Release v0.1.0"
   git tag v0.1.0
   git push origin master --tags
   ```

4. Watch <https://github.com/drchristhorpe/polars-seq/actions>. The release workflow will:
   - build wheels for Linux (x86_64, aarch64), macOS (x86_64, arm64) and Windows (x64);
   - build the sdist;
   - install the Linux wheel on a Python it was *not* built against and check it still agrees
     with BioPython (including the `XXX` rejection);
   - publish everything to PyPI.

5. Confirm it landed: <https://pypi.org/project/polars-seq/>

```bash
uv run --no-project --with polars-seq --with polars python -c "
import polars as pl, polars_seq
print(pl.DataFrame({'dna': ['ATGAAATTTTAA']}).with_columns(p=pl.col('dna').seq.translate()))"
```

---

## Notes

- **One wheel per platform, not per Python.** The extension is built against the stable ABI
  (`abi3-py310`), so a single wheel covers Python 3.10 through 3.14 and every future 3.x. This is
  why the wheel matrix has five entries rather than twenty-five.
- **The sdist is buildable.** `src/codon_tables.rs` is generated from BioPython but *committed*,
  so anyone building from source needs only Rust — not BioPython. CI enforces that the committed
  tables still match what BioPython produces (`codegen-is-current` job).
- **Test PyPI first, if you want.** Add a second pending publisher on
  <https://test.pypi.org/manage/account/publishing/> with environment `testpypi`, then add a
  matching job with `repository-url: https://test.pypi.org/legacy/`.
- **Yanking**, if a release turns out to be broken:
  <https://pypi.org/manage/project/polars-seq/releases/> → *Options* → *Yank*. Yanking hides it
  from resolvers without breaking pinned installs. You cannot re-upload the same version number,
  so bump the patch and release again.
