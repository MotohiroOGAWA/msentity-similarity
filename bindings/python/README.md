# msentity-similarity (Python)

Python binding of msentity-similarity. Guide: [docs/bindings.md](../../docs/bindings.md#python).

```python
import msentity as ms
import msentity_similarity as mss

ds = ms.load("example.msp")
hits = mss.library_search(ds, ds, threshold=0.7)
hits["cosine"]       # [1.0, 0.7597143650054932, 0.7597143650054932, 1.0]
hits.save("hits.mssim")
```

| | |
|---|---|
| Requires | Python ≥ 3.10, NumPy, the `msentity` package; pyarrow + pandas for `to_arrow()` / `to_pandas()` |
| Native library | `MSSIM_C_LIBRARY`, else `msentity_similarity/lib/`, else `<repo>/build/` |
| Files | `msentity_similarity/__init__.py` (API), `msentity_similarity/_native.py` (ctypes declarations) |

Use `entropy` for weighted spectral entropy similarity, or
`unweighted_entropy` for the unweighted score. To return both columns
alongside cosine for spectra with matching IDs:

```python
scores = mss.similarity_by_key(
    experimental, predicted,
    methods=["cosine", "reverse_cosine", "entropy", "unweighted_entropy"],
    **mss.parse_mz_tolerance("10ppm"),
)
weighted = scores["entropy"]
unweighted = scores["unweighted_entropy"]
```

Both methods also work with `pair_scores`, `all_pairs`, and `library_search`.
Use `"0.02Da"` for Da-only matching or `"0.02Da,10ppm"` to enable either
rule. Keep `intensity_exponent=1.0` to use the original relative intensities.
Close peaks in the same spectrum stay separate. See
[the calculation guide](../../docs/similarity.md#6-spectral-entropy-similarity-entropy-and-unweighted_entropy)
for peak pairing, preprocessing, and weighting details.

Development (from this directory, with msentity checked out next to this repository):

```console
cmake --build ../../build
PYTHONPATH=.:../../../msentity/bindings/python \
MSENTITY_C_LIBRARY=../../build/_deps/msentity-build/libmsentity_c.so \
  python -m unittest discover -s tests -v
```
