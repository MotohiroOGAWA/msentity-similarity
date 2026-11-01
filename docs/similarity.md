# How scores are calculated

msentity-similarity compares spectra by **cosine similarity of matched peaks**:
peaks of the two spectra are paired when their m/z values agree within a
tolerance in Da, in ppm, or either. This page defines the calculation precisely
and works through it on real data, so results can be checked by hand or
reproduced by another implementation.

## 1. Peak matching

1. A negative intensity is replaced by 0, and each intensity is raised to the
   power `intensity_exponent` (default 1).
2. Two peaks, one of each spectrum, **can pair** when

   ```text
   |mz1 − mz2| ≤ mz_tolerance_da                                  (Da rule)
   |mz1 − mz2| / ((mz1 + mz2) / 2) × 10⁶ ≤ mz_tolerance_ppm       (ppm rule)
   ```

   holds for at least one rule. A tolerance of 0 turns its rule off; at least
   one must be positive. The defaults are `mz_tolerance_da = 0.01` and
   `mz_tolerance_ppm = 0` (Da only). Both bounds are inclusive. The ppm rule
   divides by the mean of the two m/z values, so it gives the same result
   whichever spectrum comes first.
3. **Each peak pairs at most once.** Of all the pairs that can pair, the one
   with the smallest `|mz1 − mz2|` is taken first, then the next smallest whose
   peaks are both still free, and so on. Equal differences are taken in the
   original peak order of the first spectrum, then of the second.
4. Peaks with intensity 0 (after step 1) never pair; they could not change a
   score but would otherwise take a partner away from a real peak.

Examples with `mz_tolerance_da = 0.01`:

| first | second | pairs |
|---|---|---|
| 100.009 | 100.011 | 100.009–100.011 (difference 0.002) |
| 100.000 | 99.995, 100.002 | 100.000–100.002 (the closer one) |
| 100.000, 100.005 | 100.003 | 100.005–100.003 only; 100.000 stays unmatched |

Examples of the two rules for 1000 vs 1000.008 (0.008 Da, 8.0 ppm) and 100 vs
100.008 (0.008 Da, 80 ppm):

| `mz_tolerance_da` | `mz_tolerance_ppm` | 1000 / 1000.008 | 100 / 100.008 |
|---|---|---|---|
| 0.005 | 0 | no | no |
| 0 | 10 | yes | no |
| 0.005 | 10 | yes (ppm) | no |
| 0.01 | 10 | yes | yes (Da) |

Combining both rules allows, for example, a fixed tolerance at low m/z and a
relative one at high m/z.

For the two spectra of [examples/example.msp](examples/example.msp), with the defaults:

| m/z | spectrum 0 | spectrum 1 |
|---|---|---|
| 91.0542 | 245 | 99 |
| 125.0233 | 1000 | 281 |
| 185.0961 | 349 | — |
| 246.1125 | 735 | 1000 |

## 2. Cosine (`method = "cosine"`)

```text
            Σ  a · r                  (over the matched peak pairs)
cosine = ─────────────────────
           ‖a‖ · ‖r‖                  (‖x‖ = sqrt of Σ x² over all peaks of x)
```

For spectrum 0 (a) vs spectrum 1 (r):

```text
a·r = 245·99 + 1000·281 + 735·1000 = 1 040 255
‖a‖ = sqrt(245² + 1000² + 349² + 735²) = 1312.2694
‖r‖ = sqrt(99² + 281² + 1000²)         = 1043.4376
cosine = 1 040 255 / (1312.2694 · 1043.4376) = 0.75971
```

The score is symmetric: spectrum 1 vs spectrum 0 is also 0.75971.

## 3. Reverse cosine (`method = "reverse_cosine"`)

The same, except that the **query** norm counts only the query peaks that are
matched to a reference peak. Extra query peaks (for example impurities or
fragments missing from the library spectrum) then do not lower the score:

```text
                    a·r
reverse_cosine = ────────────────
                  ‖a matched‖ · ‖r‖
```

Spectrum 0 as query vs spectrum 1 as reference: the 185.0961 peak has no partner,
so

```text
‖a matched‖ = sqrt(245² + 1000² + 735²) = 1265.0099
reverse_cosine = 1 040 255 / (1265.0099 · 1043.4376) = 0.78810
```

Spectrum 1 as query vs spectrum 0 as reference: every peak of spectrum 1 is
matched, so the score equals the cosine, 0.75971. Reverse cosine is therefore
**not** symmetric.

## 4. Edge cases and output

* An empty spectrum, or one whose intensities are all 0, scores 0 against anything.
* Scores are clipped to [0, 1] (to absorb rounding) and returned as 32-bit floats.
* `mz_tolerance_da` and `mz_tolerance_ppm` must be finite and ≥ 0 with at least
  one positive, and `intensity_exponent` must be positive; otherwise the call
  fails with an invalid-argument error.
* `intensity_exponent = 0.5` (square root) reduces the weight of dominant peaks.
  For the example spectra it raises the cosine from 0.75971 to 0.86077.

## 5. Calculation modes

| Mode | Pairs compared | Result rows |
|---|---|---|
| `pair_scores` | the pairs you list: (ds1[i1[k]], ds2[i2[k]]) | one score per pair (array) |
| `all_pairs_matrix` | every spectrum of ds1 with every spectrum of ds2 | len(ds1) × len(ds2) matrix |
| `library_search` | every query × every reference spectrum | one row per pair with score ≥ threshold, ordered by query then reference |
| `similarity_by_key` | spectra whose key values are equal (`key1` in ds1 = `key2` in ds2) | one row per matching key, in ds1 order |

Notes on `similarity_by_key`:

* Each key value must occur at most once in each dataset; duplicates are an error
  (the message lists up to five of them).
* Rows whose key is missing (null) are skipped.
* Integer and floating keys compare by value (`3` matches `3.0`); strings compare
  exactly; a string never matches a number.

Each spectrum is sorted by m/z once per call; a pair is then scored by one
sliding-window pass over both peak lists (cost about the number of peaks plus
the number of candidate pairs) without allocating memory.
`library_search` and `all_pairs_matrix` use several threads when there is
enough work (from 2048 pairs). The order of the output does not depend on the
number of threads.
