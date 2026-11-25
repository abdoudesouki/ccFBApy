# Building and publishing ccFBApy

The package is defined in `pyproject.toml` (setuptools). The code is in the folder `ccfbapy` (imported as
`ccfbapy`, published as `ccFBApy`). The data shipped with it are in `ccfbapy/data` (SBML and CSV files only; the
package does not need R or `rdata`). They are made by `tools/make_package_data.py`:

```
pip install rdata                        # only for this tool, to read the .RData files
python tools/make_package_data.py        # or --no-rdata to skip the export of the R data sets
```

It copies the CSV tables of `data/` and the files of `R/sybilccFBA/inst/extdata`, compresses the SBML models of
`models/` (`ccfbapy/data/bigg/`), and exports the data sets of `R/sybilccFBA/data/*.RData` to SBML and CSV
(`ccfbapy/data/sybilccFBA/`), checking each export against the R data set. Run it again when these files change.

The sdist contains the package, its data and the tests (`MANIFEST.in`); the R sources, the uncompressed SBML
models, the notebooks and the documentation are not in it.

## 1. Version and tests

1. Set the version in `ccfbapy/__init__.py` (`__version__`; `pyproject.toml` reads it from there).
2. Run the tests:

   ```
   pip install -e ".[all,test]"
   pytest
   ```

## 2. Build and check

```
pip install build twine
rm -rf dist build
python -m build                 # dist/ccfbapy-<version>.tar.gz and dist/ccfbapy-<version>-py3-none-any.whl
twine check dist/*
```

Test the wheel in a new environment, outside the repository (so that the folder `ccfbapy` of the repository is
not imported):

```
python -m venv /tmp/venv && . /tmp/venv/bin/activate
pip install "dist/ccfbapy-<version>-py3-none-any.whl[test]"
cp -r tests /tmp/tt && cd /tmp/tt && pytest tests
```

## 3. Publish

You need an account on [PyPI](https://pypi.org) (and, for a test, on [TestPyPI](https://test.pypi.org)) and an
API token (Account settings, API tokens). The first upload creates the project `ccFBApy`; the name was free
when this package was made.

1. Optional, a test on TestPyPI:

   ```
   twine upload --repository testpypi dist/*
   pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ ccFBApy
   ```

2. PyPI:

   ```
   twine upload dist/*
   ```

   Twine asks for the user name `__token__` and the token as password (or use `~/.pypirc`, or the environment
   variables `TWINE_USERNAME=__token__` and `TWINE_PASSWORD=<token>`).

3. Check: `pip install ccFBApy` in a new environment, and `python -c "import ccfbapy; print(ccfbapy.__version__)"`.

A version that is on PyPI cannot be uploaded again: increase the version for each release.

Before the first release, you may add the address of the repository to `pyproject.toml`, so that PyPI links to
it and to the documentation:

```toml
[project.urls]
Homepage = "https://github.com/<user>/<repository>"
Documentation = "https://github.com/<user>/<repository>/tree/main/docs"
```
