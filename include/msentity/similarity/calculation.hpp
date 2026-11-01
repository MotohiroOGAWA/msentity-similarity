#pragma once

#include <arrow/table.h>

#include <cstdint>
#include <functional>
#include <memory>
#include <string>
#include <vector>

#include <msentity/msentity.hpp>

namespace msentity::similarity {

// Tolerance-based spectrum similarity (docs/similarity.md).
//
// Negative intensities become 0 and intensities are raised to
// intensity_exponent. Peaks of the two spectra are then paired one-to-one:
// two peaks can pair when
//   |mz1 - mz2| <= mz_tolerance_da   or
//   |mz1 - mz2| / ((mz1 + mz2) / 2) * 1e6 <= mz_tolerance_ppm,
// and pairs are taken in order of increasing |mz1 - mz2| (ties: original peak
// order), each peak at most once.
//   cosine          dot / sqrt(|a|² |b|²)
//   reverse_cosine  dot / sqrt(|a matched|² |b|²) — peaks of the first
//                   (query) spectrum without a partner in the second do not
//                   count, so noise in the query is tolerated.
// Scores are clipped to [0, 1]; a zero denominator gives 0.
enum class Method { Cosine, ReverseCosine };

// "cosine" / "reverse_cosine".
std::string method_name(Method method);
Method parse_method(const std::string& name);

// An m/z tolerance as text: "<value>Da" and/or "<value>ppm" joined by a comma,
// in either order, e.g. "0.01Da,10ppm", "10ppm", "10 PPM, 0.01da". Units are
// case-insensitive, spaces are allowed, and each unit may appear once; a unit
// left out is 0 (off). Values must be finite and >= 0, at least one positive.
struct MzTolerance {
  double da = 0.0;
  double ppm = 0.0;
};
MzTolerance parse_mz_tolerance(const std::string& text);
// Canonical text: "0.01Da,10ppm"; a zero term is left out ("10ppm").
std::string format_mz_tolerance(double da, double ppm);

struct Options {
  Method method = Method::Cosine;
  double mz_tolerance_da = 0.01;     // m/z tolerance in Da (>= 0; 0 = Da rule off)
  double mz_tolerance_ppm = 10.0;    // m/z tolerance in ppm (>= 0; 0 = ppm rule off);
                                     // at least one of the two must be positive
  double intensity_exponent = 1.0;   // e.g. 0.5 for square-root intensities (> 0)
  int threads = 0;                   // worker threads; 0 = all hardware threads
};

// (completed pairs, total pairs)
using ProgressCallback = std::function<void(int64_t, int64_t)>;

// Score of ds1[index1[k]] vs ds2[index2[k]] for every k (positions in the
// dataset views).
std::vector<float> pair_scores(const MSDataset& ds1, const std::vector<int64_t>& index1,
                               const MSDataset& ds2, const std::vector<int64_t>& index2,
                               const Options& options = {});

// Every spectrum of ds1 against every spectrum of ds2; row-major
// (size(ds1) × size(ds2)).
std::vector<float> all_pairs_matrix(const MSDataset& ds1, const MSDataset& ds2,
                                    const Options& options = {});

// Pair spectra whose metadata keys are equal (ds1[key1] == ds2[key2]) and
// score them. Missing keys are skipped; a key value that occurs twice in one
// dataset is an error. Columns: key1 [, key2 if different], index1, index2,
// similarity (float32); rows in ds1 order.
std::shared_ptr<arrow::Table> similarity_by_key(const MSDataset& ds1, const MSDataset& ds2,
                                                const std::string& key1 = "SpecID",
                                                const std::string& key2 = "SpecID",
                                                const Options& options = {});

// Score every query spectrum against every reference spectrum and keep the
// pairs with score >= threshold. Columns: index1 (query), index2
// (reference), similarity (float32); rows ordered by index1, index2.
std::shared_ptr<arrow::Table> library_search(const MSDataset& query, const MSDataset& reference,
                                             double threshold = 0.8,
                                             const Options& options = {},
                                             const ProgressCallback& progress = {});

}  // namespace msentity::similarity
