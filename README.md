# JITs Benchmarks
Here I will try to explain how can you reproduce the results I got for ACO presentations 2022.

## Requirements
* `gcc`
* `bash`
* `julia` (version >= 1.5)
* `python` (version >= 3.6)

One needs also the necessary packages for `julia` and `python` (TOML parser for `C` is provided by https://github.com/cktan/tomlc99).

Please run the following on `julia`:
```julia
using Pkg
Pkg.add("StaticArrays")
Pkg.add("TOML")
```

Please run the following on your favourite terminal:
```bash
python -m pip install --user numpy numba
```

## Steps to get the results
1) Enter any of the folders.
2) Execute `make`.
3) Adapt `template.toml` to your needs. DO NOT CHANGE THE LINE WHICH SAYS "PATATA".
4) Execute `bash run.sh`.

Once this is done, the results will be in:
* `ccfile` for `gcc`.
* `cofile` for `gcc -O3`.
* `jlfile` for `julia`.
* `pyfile` for `python + numpy + numba`.
* `ppfile` for `python + numpy`.

If you want some graphs, consider running (`matplotlib` and `pandas` is needed)
```bash
python plots.py
```
in each folder.