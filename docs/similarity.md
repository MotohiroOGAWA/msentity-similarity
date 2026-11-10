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
   `mz_tolerance_ppm = 10` (0.01 Da or 10 ppm; the ppm rule is the wider one
   above m/z 1000). Both bounds are inclusive. The ppm rule
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

The two tolerances can also be written as one text, as `mssim-tool
--mz-tolerance` takes them: `<value>Da` and/or `<value>ppm` joined by a comma
in either order, units case-insensitive and spaces allowed (`0.01Da,10ppm`,
`10 PPM, 0.01da`, `10ppm`). A unit left out is 0. `parse_mz_tolerance` reads
this text and `format_mz_tolerance` writes its canonical form (`0.01Da,10ppm`;
a zero term is left out), in C++, the C API, and every binding.

For the two spectra of [examples/example.msp](examples/example.msp), with the defaults:

| m/z | spectrum 0 | spectrum 1 |
|---|---|---|
| 91.0542 | 245 | 99 |
| 125.0233 | 1000 | 281 |
| 185.0961 | 349 | — |
| 246.1125 | 735 | 1000 |

## 2. Cosine (`cosine`)

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

## 3. Reverse cosine (`reverse_cosine`)

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

## 4. Modified cosine (`modified_cosine`)

Fragments that keep the part of the molecule where two related compounds
differ appear shifted by the difference of their precursor m/z. Modified
cosine therefore also pairs such peaks: with `shift = precursor1 − precursor2`,
a peak at mz1 can pair with a peak at mz2 directly (§1) **or** when mz1 and
mz2 + shift match under the same tolerance rule. Direct and shifted pairs are
taken together, smallest error (|mz1 − mz2| or |mz1 − (mz2 + shift)|) first,
each peak at most once. The score is then the cosine of these pairs, with the
norms over all peaks:

```text
                    Σ  a · r            (over the direct and shifted pairs)
modified_cosine = ─────────────
                    ‖a‖ · ‖r‖
```

With equal precursor m/z (as for the two example spectra, both 288.1225) it
equals the cosine, 0.75971.

## 5. BONANZA (`bonanza`)

Peaks are paired as for modified cosine. Each spectrum's intensities are
divided by their sum (p = a / Σa, q = r / Σr). Of these shares, p_u and q_v
denote those of the **unpaired** peaks (peaks without a partner) of the first
and the second spectrum; p·q runs over the pairs:

```text
                Σ p·q                         (p·q: over the pairs)
bonanza = ──────────────────────────────────
           Σ p·q + Σ p_u² + Σ q_v²            (p_u, q_v: unpaired peaks only)
```

Unpaired peaks lower the score through their squared share instead of through
the norms. For the example spectra (sums 2329 and 1380), the 185.0961 peak of
spectrum 0 is its only unpaired peak, and every peak of spectrum 1 has a
partner, so Σ q_v² is 0:

```text
Σ p·q   = (245·99 + 1000·281 + 735·1000) / (2329 · 1380) = 0.323662
Σ p_u²  = (349 / 2329)²                                  = 0.022455
Σ q_v²  = 0                                              (no unpaired peak in spectrum 1)
bonanza = 0.323662 / (0.323662 + 0.022455 + 0)           = 0.93512
```

### Precursor m/z

Modified cosine and BONANZA need the precursor m/z of both spectra. It is
read from the metadata column `precursor_column`, or, when that is empty (the
default), from the first of `PrecursorMZ`, `precursor_mz`, `PEPMASS`,
`ParentMZ` (names compared case-insensitively) that holds a positive number
for the spectrum. Text values are read up to the first space (`301.14 1000`
in MGF PEPMASS gives 301.14). A named column that does not exist is an error.
When either spectrum has no precursor m/z, these two scores are **missing**:
null in result tables, NaN from `pair_scores` / `all_pairs_matrix`; the other
methods are still calculated for that pair.

## 6. Several methods at once

`methods` selects one or more of the four; result tables get one score column
per method, named as above and always in the order `cosine`,
`reverse_cosine`, `modified_cosine`, `bonanza`. `library_search` keeps a pair
when **any** selected method scores at least the threshold (missing scores
never do). `pair_scores` and `all_pairs_matrix` return one score per pair and
need exactly one method.

## 7. Edge cases and output

* An empty spectrum, or one whose intensities are all 0, scores 0 against anything.
* Scores are clipped to [0, 1] (to absorb rounding) and returned as 32-bit floats.
* `mz_tolerance_da` and `mz_tolerance_ppm` must be finite and ≥ 0 with at least
  one positive, and `intensity_exponent` must be positive; otherwise the call
  fails with an invalid-argument error.
* `intensity_exponent = 0.5` (square root) reduces the weight of dominant peaks.
  For the example spectra it raises the cosine from 0.75971 to 0.86077.

## 8. Peak intensity coverage

Result tables (`library_search`, `similarity_by_key`) also report, for each
side, how much of the spectrum the matched peaks explain:

```text
                        Σ  x²   (over the matched peaks of x)
peak_intensity_coverage = ─────────
                        Σ  x²   (over all peaks of x)
```

where x are the intensities after step 1 (raised to `intensity_exponent`, so
x² = I^(2·exponent)). With `intensity_exponent = 0.5` this is simply the
matched share of the summed raw intensity. A spectrum without positive
intensity has coverage 0. The columns are `peak_intensity_coverage1` (first /
query spectrum) and `peak_intensity_coverage2` (second / reference spectrum).

By the Cauchy–Schwarz inequality the score cannot exceed what the coverages
allow:

```text
cosine          ≤ sqrt(peak_intensity_coverage1 · peak_intensity_coverage2) ≤ sqrt(min of the two)
reverse_cosine  ≤ sqrt(peak_intensity_coverage2)
```

The coverages come from the direct pairing (§1), the one cosine and reverse
cosine use. Modified cosine pairs more peaks (the shifted ones too), so it can
exceed sqrt(peak_intensity_coverage). When modified cosine or BONANZA is
calculated, the table also has `modified_peak_intensity_coverage1` /
`modified_peak_intensity_coverage2`, the same shares under their pairing (direct and
shifted pairs; null without a precursor m/z), and the same argument gives

```text
modified_cosine ≤ sqrt(modified_peak_intensity_coverage1 · modified_peak_intensity_coverage2)
```

BONANZA uses the same pairing, so these columns show how much of each spectrum
its pairs cover, but it normalizes differently (sums instead of norms), so no
such bound holds for it.

So if spectrum 1 is a measured spectrum and spectrum 2 a prediction, a cosine
above sqrt(peak_intensity_coverage1) is impossible however good the predicted
intensities are; equality needs every peak of the other side matched and the
matched intensities proportional.

For the example spectra (spectrum 0 vs spectrum 1, exponent 1), the 185.0961
peak of spectrum 0 is unmatched:

```text
peak_intensity_coverage1 = (245² + 1000² + 735²) / (245² + 1000² + 349² + 735²) = 0.92927
peak_intensity_coverage2 = 1
cosine 0.75971 ≤ sqrt(0.92927 · 1) = 0.96399
```

## 9. Calculation modes

| Mode | Pairs compared | Result rows |
|---|---|---|
| `pair_scores` | the pairs you list: (ds1[i1[k]], ds2[i2[k]]) | one score per pair (array) |
| `all_pairs_matrix` | every spectrum of ds1 with every spectrum of ds2 | len(ds1) × len(ds2) matrix |
| `library_search` | every query × every reference spectrum | one row per pair with score ≥ threshold, ordered by query then reference; with peak intensity coverages |
| `similarity_by_key` | spectra whose key values are equal (`key1` in ds1 = `key2` in ds2) | one row per matching key, in ds1 order; with peak intensity coverages |

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
