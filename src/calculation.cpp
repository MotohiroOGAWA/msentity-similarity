#include "msentity/similarity/calculation.hpp"

#include <algorithm>
#include <atomic>
#include <cctype>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <limits>
#include <locale>
#include <sstream>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>

#include "internal/table.hpp"

namespace msentity::similarity {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// A peak that can be matched; peaks with zero intensity (after clipping) or a
// non-finite m/z cannot contribute to a score and are left out.
struct Peak {
  double mz;
  double intensity;
  int64_t original_index;
};

// Peaks sorted by m/z (original order among equal m/z), the squared norm and
// the sum of all intensities, and the precursor m/z (NaN when unknown).
struct PreparedSpectrum {
  std::vector<Peak> peaks;
  double norm = 0.0;
  double total = 0.0;
  double unlisted_square = 0.0;  // Σ intensity² of positive peaks left out (no finite m/z)
  double precursor = kNaN;
};

void validate_mz_tolerance(double da, double ppm) {
  if (!std::isfinite(da) || da < 0) {
    throw std::invalid_argument("mz_tolerance_da must be finite and >= 0");
  }
  if (!std::isfinite(ppm) || ppm < 0) {
    throw std::invalid_argument("mz_tolerance_ppm must be finite and >= 0");
  }
  if (da == 0 && ppm == 0) {
    throw std::invalid_argument("at least one of mz_tolerance_da and mz_tolerance_ppm must be positive");
  }
}

void validate(const Options& options) {
  validate_mz_tolerance(options.mz_tolerance_da, options.mz_tolerance_ppm);
  if (!std::isfinite(options.intensity_exponent) || options.intensity_exponent <= 0) {
    throw std::invalid_argument("intensity_exponent must be positive");
  }
  if (options.threads < 0) throw std::invalid_argument("threads must be >= 0");
  if (options.methods.empty()) throw std::invalid_argument("methods must not be empty");
}

bool any_precursor_method(const std::vector<Method>& methods) {
  return std::any_of(methods.begin(), methods.end(), needs_precursor);
}

PreparedSpectrum prepare_spectrum(const Spectrum& spectrum, const Options& options) {
  PreparedSpectrum result;
  result.peaks.reserve(static_cast<std::size_t>(spectrum.size()));
  for (int64_t i = 0; i < spectrum.size(); ++i) {
    double intensity = std::max(0.0, spectrum.intensity(i));
    if (options.intensity_exponent != 1.0) {
      intensity = std::pow(intensity, options.intensity_exponent);
    }
    result.norm += intensity * intensity;
    result.total += intensity;
    double mz = spectrum.mz(i);
    if (intensity > 0 && std::isfinite(mz)) {
      result.peaks.push_back({mz, intensity, i});
    } else {
      result.unlisted_square += intensity * intensity;
    }
  }
  std::stable_sort(result.peaks.begin(), result.peaks.end(),
                   [](const Peak& a, const Peak& b) { return a.mz < b.mz; });
  return result;
}

std::string lower(std::string text) {
  std::transform(text.begin(), text.end(), text.begin(),
                 [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
  return text;
}

// A positive precursor m/z from a metadata value, or NaN.
double precursor_value(const Value& value) {
  double mz = kNaN;
  if (const auto* d = std::get_if<double>(&value)) mz = *d;
  if (const auto* i = std::get_if<int64_t>(&value)) mz = static_cast<double>(*i);
  if (const auto* s = std::get_if<std::string>(&value)) {
    // Leading number of the text ("301.1 1000"), read in the C locale.
    std::istringstream in(*s);
    in.imbue(std::locale::classic());
    double parsed;
    if (in >> parsed) mz = parsed;
  }
  return std::isfinite(mz) && mz > 0 ? mz : kNaN;
}

// Precursor m/z of every row of dataset (NaN when unknown): from
// options.precursor_column, or the first of precursor_columns() that has one.
std::vector<double> precursors(const MSDataset& dataset, const Options& options) {
  std::vector<double> result(static_cast<std::size_t>(dataset.size()), kNaN);
  std::vector<std::string> names;
  if (!options.precursor_column.empty()) {
    if (!dataset.has_column(options.precursor_column)) {
      throw KeyError("precursor column not found: " + options.precursor_column);
    }
    names.push_back(options.precursor_column);
  } else {
    for (const auto& wanted : precursor_columns()) {
      for (const auto& column : dataset.columns()) {
        if (lower(column) == lower(wanted)) {
          names.push_back(column);
          break;
        }
      }
    }
  }
  for (const auto& name : names) {
    auto values = dataset.column_values(name);
    for (std::size_t i = 0; i < result.size(); ++i) {
      if (std::isnan(result[i])) result[i] = precursor_value(values[i]);
    }
  }
  return result;
}

// The m/z match rule: |mz1 - mz2| <= Da tolerance, or the difference relative
// to the mean m/z <= ppm tolerance. Bounds are inclusive; a few ulps of slack
// keep e.g. 100.00 vs 100.01 at 0.01 Da inside despite rounding.
class Tolerance {
 public:
  explicit Tolerance(const Options& options)
      : da_(options.mz_tolerance_da), ppm_(options.mz_tolerance_ppm * 1e-6) {}

  bool matches(double mz1, double mz2) const {
    double difference = std::abs(mz1 - mz2);
    double slack = kSlack * std::max(std::abs(mz1), std::abs(mz2));
    if (da_ > 0 && difference <= da_ + slack) return true;
    return ppm_ > 0 && difference <= ppm_ * (mz1 + mz2) / 2.0 + slack;
  }

  // [lower, upper] containing every m/z that matches mz; both bounds grow
  // with mz, so a sorted spectrum can be scanned with a sliding window.
  std::pair<double, double> window(double mz) const {
    double lower = mz - da_, upper = mz + da_;
    if (ppm_ > 0 && mz > 0) {
      double r = ppm_ / 2.0;  // |mz - x| <= r (mz + x)
      lower = std::min(lower, mz * (1.0 - r) / (1.0 + r));
      upper = r >= 1.0 ? std::numeric_limits<double>::infinity()
                       : std::max(upper, mz * (1.0 + r) / (1.0 - r));
    }
    double pad = 1e-12 * (std::abs(mz) + 1.0);
    return {lower - pad, upper + pad};
  }

 private:
  static constexpr double kSlack = 8 * std::numeric_limits<double>::epsilon();
  double da_;
  double ppm_;
};

struct Candidate {
  double difference;
  uint32_t first;   // positions in the sorted peaks
  uint32_t second;
};

// Scores of a pair (NaN: not calculated or not available) and the intensity
// coverage of each side under the direct pairing.
struct Scored {
  double cosine = kNaN;
  double reverse_cosine = kNaN;
  double modified_cosine = kNaN;
  double bonanza = kNaN;
  double coverage1 = 0.0;
  double coverage2 = 0.0;

  double get(Method method) const {
    switch (method) {
      case Method::Cosine: return cosine;
      case Method::ReverseCosine: return reverse_cosine;
      case Method::ModifiedCosine: return modified_cosine;
      case Method::Bonanza: return bonanza;
    }
    return kNaN;
  }
};

double clip_unit(double value) {
  if (!(value >= 0.0)) value = std::isnan(value) ? value : 0.0;  // keep NaN, clip negatives
  if (value > 1.0) value = 1.0;
  return value;
}

// Candidate pairs of a with b shifted by offset (a.mz ≈ b.mz + offset),
// appended to candidates; difference = |a.mz - (b.mz + offset)|.
void collect(const PreparedSpectrum& a, const PreparedSpectrum& b, double offset,
             const Tolerance& tolerance, std::vector<Candidate>& candidates) {
  if (a.peaks.empty() || b.peaks.empty()) return;
  if (tolerance.window(a.peaks.back().mz).second < b.peaks.front().mz + offset ||
      tolerance.window(a.peaks.front().mz).first > b.peaks.back().mz + offset) {
    return;  // m/z ranges too far apart
  }
  std::size_t start = 0;
  for (std::size_t i = 0; i < a.peaks.size(); ++i) {
    double mz = a.peaks[i].mz;
    auto [low, high] = tolerance.window(mz);
    while (start < b.peaks.size() && b.peaks[start].mz + offset < low) ++start;
    for (std::size_t j = start; j < b.peaks.size() && b.peaks[j].mz + offset <= high; ++j) {
      double shifted = b.peaks[j].mz + offset;
      if (tolerance.matches(mz, shifted)) {
        candidates.push_back({std::abs(mz - shifted), static_cast<uint32_t>(i),
                              static_cast<uint32_t>(j)});
      }
    }
  }
}

// Sums over one-to-one pairs chosen from candidates.
struct Pairing {
  double dot = 0.0;
  double matched_norm_a = 0.0, matched_norm_b = 0.0;
  double bonanza_unmatched = 0.0;                        // Σ p² + Σ q² of unmatched peaks
  double bonanza_dot = 0.0;                              // Σ p·q
};

// Pairs are taken in order of increasing m/z difference (ties: original peak
// order of the first, then the second spectrum), skipping pairs whose peaks
// are already used.
Pairing pair_peaks(const PreparedSpectrum& a, const PreparedSpectrum& b,
                   std::vector<Candidate>& candidates, bool bonanza) {
  thread_local std::vector<char> used_a, used_b;
  std::sort(candidates.begin(), candidates.end(), [&](const Candidate& x, const Candidate& y) {
    if (x.difference != y.difference) return x.difference < y.difference;
    if (x.first != y.first) return a.peaks[x.first].original_index < a.peaks[y.first].original_index;
    return b.peaks[x.second].original_index < b.peaks[y.second].original_index;
  });
  used_a.assign(a.peaks.size(), 0);
  used_b.assign(b.peaks.size(), 0);
  Pairing result;
  for (const auto& candidate : candidates) {
    if (used_a[candidate.first] || used_b[candidate.second]) continue;
    used_a[candidate.first] = used_b[candidate.second] = 1;
    double x = a.peaks[candidate.first].intensity;
    double y = b.peaks[candidate.second].intensity;
    result.dot += x * y;
    result.matched_norm_a += x * x;
    result.matched_norm_b += y * y;
  }
  if (bonanza && a.total > 0 && b.total > 0) {
    for (std::size_t i = 0; i < a.peaks.size(); ++i) {
      double p = a.peaks[i].intensity / a.total;
      if (!used_a[i]) result.bonanza_unmatched += p * p;
    }
    for (std::size_t j = 0; j < b.peaks.size(); ++j) {
      double q = b.peaks[j].intensity / b.total;
      if (!used_b[j]) result.bonanza_unmatched += q * q;
    }
    result.bonanza_unmatched += a.unlisted_square / (a.total * a.total) +
                                b.unlisted_square / (b.total * b.total);
    result.bonanza_dot = result.dot / (a.total * b.total);
  }
  return result;
}

// Scores of the selected methods for one pair.
Scored score_pair(const PreparedSpectrum& a, const PreparedSpectrum& b, const Tolerance& tolerance,
                  const std::vector<Method>& methods) {
  // Per-thread buffers: scoring a pair allocates nothing once they have grown.
  thread_local std::vector<Candidate> candidates;
  auto wanted = [&](Method m) { return std::find(methods.begin(), methods.end(), m) != methods.end(); };
  Scored result;

  // Direct pairing: cosine, reverse cosine, and the coverages.
  candidates.clear();
  collect(a, b, 0.0, tolerance, candidates);
  Pairing direct = pair_peaks(a, b, candidates, false);
  double denominator = std::sqrt(a.norm * b.norm);
  if (wanted(Method::Cosine)) result.cosine = clip_unit(denominator > 0 ? direct.dot / denominator : 0.0);
  if (wanted(Method::ReverseCosine)) {
    double reverse = std::sqrt(direct.matched_norm_a * b.norm);
    result.reverse_cosine = clip_unit(reverse > 0 ? direct.dot / reverse : 0.0);
  }
  result.coverage1 = clip_unit(a.norm > 0 ? direct.matched_norm_a / a.norm : 0.0);
  result.coverage2 = clip_unit(b.norm > 0 ? direct.matched_norm_b / b.norm : 0.0);

  // Direct and precursor-shifted pairing: modified cosine and BONANZA.
  bool modified = wanted(Method::ModifiedCosine), bonanza = wanted(Method::Bonanza);
  if ((modified || bonanza) && std::isfinite(a.precursor) && std::isfinite(b.precursor)) {
    candidates.clear();
    collect(a, b, 0.0, tolerance, candidates);
    double shift = a.precursor - b.precursor;
    if (shift != 0.0) collect(a, b, shift, tolerance, candidates);
    Pairing shifted = pair_peaks(a, b, candidates, bonanza);
    if (modified) result.modified_cosine = clip_unit(denominator > 0 ? shifted.dot / denominator : 0.0);
    if (bonanza) {
      double total = shifted.bonanza_dot + shifted.bonanza_unmatched;
      result.bonanza = clip_unit(total > 0 ? shifted.bonanza_dot / total : 0.0);
    }
  }
  return result;
}

std::vector<PreparedSpectrum> prepare_all(const MSDataset& dataset, const Options& options,
                                          const std::vector<Method>& methods) {
  std::vector<PreparedSpectrum> prepared(static_cast<std::size_t>(dataset.size()));
  std::vector<double> precursor;
  if (any_precursor_method(methods)) precursor = precursors(dataset, options);
  for (int64_t i = 0; i < dataset.size(); ++i) {
    auto& spectrum = prepared[static_cast<std::size_t>(i)];
    spectrum = prepare_spectrum(dataset.peaks().spectrum(i), options);
    if (!precursor.empty()) spectrum.precursor = precursor[static_cast<std::size_t>(i)];
  }
  return prepared;
}

unsigned thread_count(const Options& options, int64_t work) {
  unsigned threads = options.threads > 0 ? static_cast<unsigned>(options.threads)
                                         : std::max(1u, std::thread::hardware_concurrency());
  if (work < 2048) threads = 1;  // not worth starting threads
  return threads;
}

// Run body(begin, end) over [0, n) split into contiguous ranges, one per thread.
template <typename Body>
void parallel_ranges(int64_t n, unsigned threads, Body&& body) {
  if (threads <= 1 || n <= 1) {
    body(int64_t{0}, n);
    return;
  }
  std::vector<std::thread> workers;
  int64_t chunk = (n + threads - 1) / threads;
  for (unsigned t = 0; t < threads; ++t) {
    int64_t begin = static_cast<int64_t>(t) * chunk;
    int64_t end = std::min(n, begin + chunk);
    if (begin >= end) break;
    workers.emplace_back([&body, begin, end] { body(begin, end); });
  }
  for (auto& worker : workers) worker.join();
}

void check_positions(const std::vector<int64_t>& index, int64_t size, const char* name) {
  for (int64_t position : index) {
    if (position < 0 || position >= size) {
      throw std::out_of_range(std::string(name) + " contains an invalid position: " +
                              std::to_string(position));
    }
  }
}

// Key of a metadata value for equality: numbers (and booleans) compare by
// value across integer/float types, text by content.
std::optional<std::string> match_key(const Value& value) {
  if (is_null(value)) return std::nullopt;
  if (const auto* d = std::get_if<double>(&value)) {
    if (std::isnan(*d)) return std::nullopt;
    return "n:" + format_float(*d);
  }
  if (const auto* i = std::get_if<int64_t>(&value)) return "n:" + format_float(static_cast<double>(*i));
  if (const auto* b = std::get_if<bool>(&value)) return "n:" + format_float(*b ? 1.0 : 0.0);
  return "s:" + std::get<std::string>(value);
}

std::shared_ptr<arrow::Array> int64_array(const std::vector<int64_t>& values) {
  arrow::Int64Builder builder;
  internal::check(builder.AppendValues(values));
  return internal::unwrap(builder.Finish());
}

// float32 values; NaN becomes null.
std::shared_ptr<arrow::Array> float_array(const std::vector<double>& values) {
  arrow::FloatBuilder builder;
  internal::check(builder.Reserve(static_cast<int64_t>(values.size())));
  for (double value : values) {
    if (std::isnan(value)) {
      builder.UnsafeAppendNull();
    } else {
      builder.UnsafeAppend(static_cast<float>(value));
    }
  }
  return internal::unwrap(builder.Finish());
}

// Result columns of scored pairs: one per method, then the coverages.
void add_score_columns(const std::vector<Scored>& scored, const std::vector<Method>& methods,
                       std::vector<std::shared_ptr<arrow::Field>>& fields,
                       std::vector<std::shared_ptr<arrow::Array>>& arrays) {
  auto add = [&](const std::string& name, const std::vector<double>& values) {
    auto array = float_array(values);
    fields.push_back(arrow::field(name, arrow::float32()));
    arrays.push_back(array);
  };
  for (Method method : methods) {
    std::vector<double> values;
    values.reserve(scored.size());
    for (const auto& s : scored) values.push_back(s.get(method));
    add(method_name(method), values);
  }
  std::vector<double> coverage1, coverage2;
  for (const auto& s : scored) {
    coverage1.push_back(s.coverage1);
    coverage2.push_back(s.coverage2);
  }
  add("intensity_coverage1", coverage1);
  add("intensity_coverage2", coverage2);
}

// Only method of options, for the single-score functions.
std::vector<Method> single_method(const Options& options) {
  auto methods = selected_methods(options);
  if (methods.size() != 1) {
    throw std::invalid_argument("pair_scores and all_pairs_matrix need exactly one method");
  }
  return methods;
}

}  // namespace

const std::vector<Method>& all_methods() {
  static const std::vector<Method> methods = {Method::Cosine, Method::ReverseCosine,
                                              Method::ModifiedCosine, Method::Bonanza};
  return methods;
}

std::string method_name(Method method) {
  switch (method) {
    case Method::Cosine: return "cosine";
    case Method::ReverseCosine: return "reverse_cosine";
    case Method::ModifiedCosine: return "modified_cosine";
    case Method::Bonanza: return "bonanza";
  }
  throw std::invalid_argument("unknown method");
}

Method parse_method(const std::string& name) {
  for (Method method : all_methods()) {
    if (method_name(method) == name) return method;
  }
  throw std::invalid_argument(
      "method must be 'cosine', 'reverse_cosine', 'modified_cosine', or 'bonanza', got '" + name + "'");
}

bool needs_precursor(Method method) {
  return method == Method::ModifiedCosine || method == Method::Bonanza;
}

const std::vector<std::string>& precursor_columns() {
  static const std::vector<std::string> names = {"PrecursorMZ", "precursor_mz", "PEPMASS", "ParentMZ"};
  return names;
}

std::vector<Method> selected_methods(const Options& options) {
  if (options.methods.empty()) throw std::invalid_argument("methods must not be empty");
  std::vector<Method> result;
  for (Method method : all_methods()) {
    if (std::find(options.methods.begin(), options.methods.end(), method) != options.methods.end()) {
      result.push_back(method);
    }
  }
  return result;
}

namespace {

std::string trim(const std::string& text) {
  auto space = [](unsigned char c) { return std::isspace(c) != 0; };
  auto begin = std::find_if_not(text.begin(), text.end(), space);
  auto end = std::find_if_not(text.rbegin(), std::string::const_reverse_iterator(begin), space).base();
  return std::string(begin, end);
}

// A decimal number filling all of text, independent of the global locale.
bool read_number(const std::string& text, double& value) {
  std::istringstream in(text);
  in.imbue(std::locale::classic());
  in >> value;
  return !in.fail() && in.peek() == std::char_traits<char>::eof();
}

// Shortest of 15 or 17 significant digits that reads back as value.
std::string write_number(double value) {
  std::string text;
  for (int precision : {15, 17}) {
    std::ostringstream out;
    out.imbue(std::locale::classic());
    out << std::setprecision(precision) << value;
    text = out.str();
    double back;
    if (read_number(text, back) && back == value) break;
  }
  return text;
}

}  // namespace

MzTolerance parse_mz_tolerance(const std::string& text) {
  auto invalid = [&] {
    return std::invalid_argument("m/z tolerance must look like '0.01Da', '10ppm', or '0.01Da,10ppm', got '" +
                                 text + "'");
  };
  MzTolerance result;
  bool has_da = false, has_ppm = false;
  std::size_t start = 0;
  while (true) {
    std::size_t comma = text.find(',', start);
    std::string term = trim(text.substr(start, comma == std::string::npos ? comma : comma - start));
    std::size_t split = term.size();
    while (split > 0 && std::isalpha(static_cast<unsigned char>(term[split - 1]))) --split;
    std::string unit = term.substr(split);
    std::transform(unit.begin(), unit.end(), unit.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    double value;
    if (!read_number(trim(term.substr(0, split)), value)) throw invalid();
    if (unit == "da" && !has_da) {
      has_da = true;
      result.da = value;
    } else if (unit == "ppm" && !has_ppm) {
      has_ppm = true;
      result.ppm = value;
    } else {
      throw invalid();
    }
    if (comma == std::string::npos) break;
    start = comma + 1;
  }
  validate_mz_tolerance(result.da, result.ppm);
  return result;
}

std::string format_mz_tolerance(double da, double ppm) {
  validate_mz_tolerance(da, ppm);
  std::string text;
  if (da > 0) text = write_number(da) + "Da";
  if (ppm > 0) text += (text.empty() ? "" : ",") + write_number(ppm) + "ppm";
  return text;
}

std::vector<Method> parse_methods(const std::string& names) {
  std::vector<Method> methods;
  std::size_t start = 0;
  while (true) {
    std::size_t comma = names.find(',', start);
    std::string name = trim(names.substr(start, comma == std::string::npos ? comma : comma - start));
    Method method = parse_method(name);
    if (std::find(methods.begin(), methods.end(), method) == methods.end()) methods.push_back(method);
    if (comma == std::string::npos) break;
    start = comma + 1;
  }
  return methods;
}

std::vector<float> pair_scores(const MSDataset& ds1, const std::vector<int64_t>& index1,
                               const MSDataset& ds2, const std::vector<int64_t>& index2,
                               const Options& options) {
  validate(options);
  const auto methods = single_method(options);
  if (index1.size() != index2.size()) {
    throw std::invalid_argument("index1 and index2 must have the same length");
  }
  check_positions(index1, ds1.size(), "index1");
  check_positions(index2, ds2.size(), "index2");

  // Prepare only the spectra that are used.
  auto prepare_used = [&](const MSDataset& dataset, const std::vector<int64_t>& index) {
    std::vector<double> precursor;
    if (any_precursor_method(methods)) precursor = precursors(dataset, options);
    std::unordered_map<int64_t, PreparedSpectrum> prepared;
    for (int64_t position : index) {
      if (prepared.count(position) == 0) {
        auto spectrum = prepare_spectrum(dataset.peaks().spectrum(position), options);
        if (!precursor.empty()) spectrum.precursor = precursor[static_cast<std::size_t>(position)];
        prepared.emplace(position, std::move(spectrum));
      }
    }
    return prepared;
  };
  auto prepared1 = prepare_used(ds1, index1);
  auto prepared2 = prepare_used(ds2, index2);
  const Tolerance tolerance(options);

  const int64_t n = static_cast<int64_t>(index1.size());
  std::vector<float> scores(static_cast<std::size_t>(n));
  parallel_ranges(n, thread_count(options, n), [&](int64_t begin, int64_t end) {
    for (int64_t k = begin; k < end; ++k) {
      auto u = static_cast<std::size_t>(k);
      scores[u] = static_cast<float>(
          score_pair(prepared1.at(index1[u]), prepared2.at(index2[u]), tolerance, methods)
              .get(methods.front()));
    }
  });
  return scores;
}

std::vector<float> all_pairs_matrix(const MSDataset& ds1, const MSDataset& ds2,
                                    const Options& options) {
  validate(options);
  const auto methods = single_method(options);
  auto prepared1 = prepare_all(ds1, options, methods);
  auto prepared2 = prepare_all(ds2, options, methods);
  const Tolerance tolerance(options);
  const int64_t n1 = ds1.size(), n2 = ds2.size();
  std::vector<float> matrix(static_cast<std::size_t>(n1 * n2));
  parallel_ranges(n1, thread_count(options, n1 * n2), [&](int64_t begin, int64_t end) {
    for (int64_t i = begin; i < end; ++i) {
      for (int64_t j = 0; j < n2; ++j) {
        matrix[static_cast<std::size_t>(i * n2 + j)] = static_cast<float>(
            score_pair(prepared1[static_cast<std::size_t>(i)], prepared2[static_cast<std::size_t>(j)],
                       tolerance, methods)
                .get(methods.front()));
      }
    }
  });
  return matrix;
}

std::shared_ptr<arrow::Table> similarity_by_key(const MSDataset& ds1, const MSDataset& ds2,
                                                const std::string& key1,
                                                const std::string& key2,
                                                const Options& options) {
  validate(options);
  const auto methods = selected_methods(options);
  if (!ds1.has_column(key1)) throw KeyError("ds1 does not contain key column: " + key1);
  if (!ds2.has_column(key2)) throw KeyError("ds2 does not contain key column: " + key2);
  std::vector<std::string> reserved = {"index1", "index2", "similarity", "intensity_coverage1",
                                       "intensity_coverage2"};
  for (Method method : all_methods()) reserved.push_back(method_name(method));
  for (const auto& key : {key1, key2}) {
    if (std::find(reserved.begin(), reserved.end(), key) != reserved.end()) {
      std::string list;
      for (const auto& name : reserved) list += (list.empty() ? "" : ", ") + name;
      throw std::invalid_argument("Key column names must not collide with result columns: " + list);
    }
  }
  auto values1 = ds1.column_values(key1);
  auto values2 = ds2.column_values(key2);

  // Non-missing keys must be unique on each side.
  auto unique_keys = [](const std::vector<Value>& values, const std::string& side,
                        const std::string& key) {
    std::unordered_map<std::string, int64_t> positions;
    std::vector<std::string> duplicates;
    for (std::size_t i = 0; i < values.size(); ++i) {
      auto k = match_key(values[i]);
      if (!k) continue;
      auto [it, inserted] = positions.emplace(*k, static_cast<int64_t>(i));
      if (!inserted) {
        std::string shown = to_display_string(values[i]);
        if (std::find(duplicates.begin(), duplicates.end(), shown) == duplicates.end() &&
            duplicates.size() < 5) {
          duplicates.push_back(shown);
        }
      }
    }
    if (!duplicates.empty()) {
      std::string examples;
      for (const auto& d : duplicates) examples += (examples.empty() ? "'" : ", '") + d + "'";
      throw std::invalid_argument(side + " key column '" + key +
                                  "' contains duplicate values; one-to-one matching is "
                                  "ambiguous; examples: [" + examples + "]");
    }
    return positions;
  };
  unique_keys(values1, "ds1", key1);
  auto right = unique_keys(values2, "ds2", key2);

  std::vector<int64_t> index1, index2;
  std::vector<Value> matched1, matched2;
  for (std::size_t i = 0; i < values1.size(); ++i) {
    auto k = match_key(values1[i]);
    if (!k) continue;
    auto it = right.find(*k);
    if (it == right.end()) continue;
    index1.push_back(static_cast<int64_t>(i));
    index2.push_back(it->second);
    matched1.push_back(values1[i]);
    matched2.push_back(values2[static_cast<std::size_t>(it->second)]);
  }

  // Score the pairs, preparing each used spectrum once.
  std::vector<double> precursor1, precursor2;
  if (any_precursor_method(methods)) {
    precursor1 = precursors(ds1, options);
    precursor2 = precursors(ds2, options);
  }
  const Tolerance tolerance(options);
  const int64_t n = static_cast<int64_t>(index1.size());
  std::vector<Scored> scored(static_cast<std::size_t>(n));
  parallel_ranges(n, thread_count(options, n), [&](int64_t begin, int64_t end) {
    for (int64_t k = begin; k < end; ++k) {
      auto u = static_cast<std::size_t>(k);
      auto a = prepare_spectrum(ds1.peaks().spectrum(index1[u]), options);
      auto b = prepare_spectrum(ds2.peaks().spectrum(index2[u]), options);
      if (!precursor1.empty()) {
        a.precursor = precursor1[static_cast<std::size_t>(index1[u])];
        b.precursor = precursor2[static_cast<std::size_t>(index2[u])];
      }
      scored[u] = score_pair(a, b, tolerance, methods);
    }
  });

  std::vector<std::shared_ptr<arrow::Field>> fields;
  std::vector<std::shared_ptr<arrow::Array>> arrays;
  auto add = [&](const std::string& name, std::shared_ptr<arrow::Array> array) {
    fields.push_back(arrow::field(name, array->type()));
    arrays.push_back(std::move(array));
  };
  add(key1, internal::values_to_array(matched1));
  if (key2 != key1) add(key2, internal::values_to_array(matched2));
  add("index1", int64_array(index1));
  add("index2", int64_array(index2));
  add_score_columns(scored, methods, fields, arrays);
  return arrow::Table::Make(arrow::schema(fields), arrays, n);
}

std::shared_ptr<arrow::Table> library_search(const MSDataset& query, const MSDataset& reference,
                                             double threshold, const Options& options,
                                             const ProgressCallback& progress) {
  validate(options);
  const auto methods = selected_methods(options);
  if (!std::isfinite(threshold) || threshold < 0 || threshold > 1) {
    throw std::invalid_argument("threshold must be a finite value between 0 and 1");
  }
  const int64_t n_query = query.size(), n_reference = reference.size();
  const int64_t total = n_query * n_reference;
  auto prepared_query = prepare_all(query, options, methods);
  auto prepared_reference = prepare_all(reference, options, methods);
  const Tolerance tolerance(options);
  // Scores are stored as float32: compare them as stored.
  auto kept = [&](const Scored& s) {
    for (Method method : methods) {
      double value = s.get(method);
      if (!std::isnan(value) && static_cast<float>(value) >= threshold) return true;
    }
    return false;
  };

  // Hits of each query row, filled in parallel, joined in order.
  std::vector<std::vector<std::pair<int64_t, Scored>>> hits(static_cast<std::size_t>(n_query));
  const unsigned threads = thread_count(options, total);
  const int64_t batch = std::max<int64_t>(1, 4096 / std::max<int64_t>(1, n_reference) + 1) * threads;
  if (progress) progress(0, total);
  for (int64_t start = 0; start < n_query; start += batch) {
    int64_t stop = std::min(n_query, start + batch);
    parallel_ranges(stop - start, threads, [&](int64_t begin, int64_t end) {
      for (int64_t q = start + begin; q < start + end; ++q) {
        auto& row = hits[static_cast<std::size_t>(q)];
        for (int64_t r = 0; r < n_reference; ++r) {
          Scored s = score_pair(prepared_query[static_cast<std::size_t>(q)],
                                prepared_reference[static_cast<std::size_t>(r)], tolerance, methods);
          if (kept(s)) row.emplace_back(r, s);
        }
      }
    });
    if (progress) progress(stop * n_reference, total);
  }

  std::vector<int64_t> index1, index2;
  std::vector<Scored> scored;
  for (int64_t q = 0; q < n_query; ++q) {
    for (const auto& [r, s] : hits[static_cast<std::size_t>(q)]) {
      index1.push_back(q);
      index2.push_back(r);
      scored.push_back(s);
    }
  }
  std::vector<std::shared_ptr<arrow::Field>> fields = {arrow::field("index1", arrow::int64()),
                                                       arrow::field("index2", arrow::int64())};
  std::vector<std::shared_ptr<arrow::Array>> arrays = {int64_array(index1), int64_array(index2)};
  add_score_columns(scored, methods, fields, arrays);
  return arrow::Table::Make(arrow::schema(fields), arrays, static_cast<int64_t>(scored.size()));
}

}  // namespace msentity::similarity
