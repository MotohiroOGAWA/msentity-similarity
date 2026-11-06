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
//   cosine           dot / sqrt(|a|² |b|²)
//   reverse_cosine   dot / sqrt(|a matched|² |b|²) — peaks of the first
//                    (query) spectrum without a partner in the second do not
//                    count, so noise in the query is tolerated.
// The precursor methods also pair peaks whose m/z differ by the precursor
// m/z difference (mz1 - mz2 ≈ precursor1 - precursor2, the same tolerance
// rule applied to mz1 and mz2 + shift), taking direct and shifted pairs
// together in order of increasing error:
//   modified_cosine  dot / sqrt(|a|² |b|²) over these pairs
//   bonanza          with p = a / Σa and q = b / Σb, and p_u, q_v those of
//                    the unpaired peaks of each spectrum:
//                    Σ p·q (pairs) / (Σ p·q + Σ p_u² + Σ q_v²)
// They need the precursor m/z of both spectra; without it their score is
// missing (null in result tables, NaN from pair_scores).
// Scores are clipped to [0, 1]; a zero denominator gives 0.
//
// Result tables also report the intensity coverage of each side, from the
// direct pairing (as for cosine):
//   intensity_coverage = Σ_matched I^(2·exponent) / Σ_all I^(2·exponent),
// the share of the spectrum's squared norm carried by its matched peaks (0
// for a spectrum without positive intensity). By Cauchy–Schwarz,
//   cosine         <= sqrt(intensity_coverage1 · intensity_coverage2)
//   reverse_cosine <= sqrt(intensity_coverage2).
// With modified_cosine or bonanza, the tables also have
// modified_intensity_coverage1/2: the same shares under their pairing (direct
// and shifted pairs; null without a precursor m/z), and
//   modified_cosine <= sqrt(modified_intensity_coverage1 · modified_intensity_coverage2).
// No such bound holds for bonanza.
enum class Method { Cosine, ReverseCosine, ModifiedCosine, Bonanza };

// Every method, in the order of result columns.
const std::vector<Method>& all_methods();
// "cosine", "reverse_cosine", "modified_cosine", "bonanza" (also the names of
// the result columns).
std::string method_name(Method method);
Method parse_method(const std::string& name);
// True for the methods that need a precursor m/z (modified_cosine, bonanza).
bool needs_precursor(Method method);
// Methods from comma-separated names ("cosine,bonanza"; spaces allowed).
std::vector<Method> parse_methods(const std::string& names);

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

// Metadata columns tried, in order, when precursor_column is "" (names
// compared case-insensitively; the first with a positive number is used).
const std::vector<std::string>& precursor_columns();

struct Options {
  // Scores to calculate (at least one). Result tables have one column per
  // method, in the order of all_methods() whatever the order here.
  std::vector<Method> methods = {Method::Cosine};
  double mz_tolerance_da = 0.01;     // m/z tolerance in Da (>= 0; 0 = Da rule off)
  double mz_tolerance_ppm = 10.0;    // m/z tolerance in ppm (>= 0; 0 = ppm rule off);
                                     // at least one of the two must be positive
  double intensity_exponent = 1.0;   // e.g. 0.5 for square-root intensities (> 0)
  // Metadata column with the precursor m/z (modified_cosine, bonanza); ""
  // tries precursor_columns(). Numbers, or text starting with a number
  // ("301.1 1000" in MGF PEPMASS), are read.
  std::string precursor_column;
  int threads = 0;                   // worker threads; 0 = all hardware threads
};

// The methods of options without duplicates, in the order of all_methods().
std::vector<Method> selected_methods(const Options& options);

// (completed pairs, total pairs)
using ProgressCallback = std::function<void(int64_t, int64_t)>;

// Score of ds1[index1[k]] vs ds2[index2[k]] for every k (positions in the
// dataset views). options must select exactly one method.
std::vector<float> pair_scores(const MSDataset& ds1, const std::vector<int64_t>& index1,
                               const MSDataset& ds2, const std::vector<int64_t>& index2,
                               const Options& options = {});

// Every spectrum of ds1 against every spectrum of ds2; row-major
// (size(ds1) × size(ds2)). options must select exactly one method.
std::vector<float> all_pairs_matrix(const MSDataset& ds1, const MSDataset& ds2,
                                    const Options& options = {});

// Pair spectra whose metadata keys are equal (ds1[key1] == ds2[key2]) and
// score them. Missing keys are skipped; a key value that occurs twice in one
// dataset is an error. Columns: key1 [, key2 if different], index1, index2,
// one float32 column per selected method (named by method_name),
// intensity_coverage1, intensity_coverage2, and with modified_cosine or bonanza
// modified_intensity_coverage1, modified_intensity_coverage2 (float32); rows
// in ds1 order.
std::shared_ptr<arrow::Table> similarity_by_key(const MSDataset& ds1, const MSDataset& ds2,
                                                const std::string& key1 = "SpecID",
                                                const std::string& key2 = "SpecID",
                                                const Options& options = {});

// Score every query spectrum against every reference spectrum and keep the
// pairs where any selected method scores >= threshold. Columns: index1
// (query), index2 (reference), one float32 column per selected method,
// the coverage columns (as similarity_by_key); rows ordered by index1, index2.
std::shared_ptr<arrow::Table> library_search(const MSDataset& query, const MSDataset& reference,
                                             double threshold = 0.8,
                                             const Options& options = {},
                                             const ProgressCallback& progress = {});

}  // namespace msentity::similarity
