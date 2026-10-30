#include "msentity/similarity/calculation.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>

#include "internal/table.hpp"

namespace msentity::similarity {

namespace {

// A peak that can be matched; peaks with zero intensity (after clipping) or a
// non-finite m/z cannot contribute to a score and are left out.
struct Peak {
  double mz;
  double intensity;
  int64_t original_index;
};

// Peaks sorted by m/z (original order among equal m/z) and the squared norm
// of all intensities.
struct PreparedSpectrum {
  std::vector<Peak> peaks;
  double norm = 0.0;
};

void validate(const Options& options) {
  if (!std::isfinite(options.mz_tolerance_da) || options.mz_tolerance_da < 0) {
    throw std::invalid_argument("mz_tolerance_da must be finite and >= 0");
  }
  if (!std::isfinite(options.mz_tolerance_ppm) || options.mz_tolerance_ppm < 0) {
    throw std::invalid_argument("mz_tolerance_ppm must be finite and >= 0");
  }
  if (options.mz_tolerance_da == 0 && options.mz_tolerance_ppm == 0) {
    throw std::invalid_argument("at least one of mz_tolerance_da and mz_tolerance_ppm must be positive");
  }
  if (!std::isfinite(options.intensity_exponent) || options.intensity_exponent <= 0) {
    throw std::invalid_argument("intensity_exponent must be positive");
  }
  if (options.threads < 0) throw std::invalid_argument("threads must be >= 0");
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
    double mz = spectrum.mz(i);
    if (intensity > 0 && std::isfinite(mz)) result.peaks.push_back({mz, intensity, i});
  }
  std::stable_sort(result.peaks.begin(), result.peaks.end(),
                   [](const Peak& a, const Peak& b) { return a.mz < b.mz; });
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

// Peaks are paired one-to-one: candidate pairs within tolerance are taken in
// order of increasing m/z difference (ties: original peak order of the first,
// then the second spectrum), skipping pairs whose peaks are already used.
float score(const PreparedSpectrum& a, const PreparedSpectrum& b, const Tolerance& tolerance,
            Method method) {
  // Per-thread buffers: scoring a pair allocates nothing once they have grown.
  thread_local std::vector<Candidate> candidates;
  thread_local std::vector<char> used_a, used_b;
  candidates.clear();

  if (a.peaks.empty() || b.peaks.empty()) return 0.0f;
  if (tolerance.window(a.peaks.back().mz).second < b.peaks.front().mz ||
      tolerance.window(a.peaks.front().mz).first > b.peaks.back().mz) {
    return 0.0f;  // m/z ranges too far apart
  }
  std::size_t start = 0;
  for (std::size_t i = 0; i < a.peaks.size(); ++i) {
    double mz = a.peaks[i].mz;
    auto [lower, upper] = tolerance.window(mz);
    while (start < b.peaks.size() && b.peaks[start].mz < lower) ++start;
    for (std::size_t j = start; j < b.peaks.size() && b.peaks[j].mz <= upper; ++j) {
      if (tolerance.matches(mz, b.peaks[j].mz)) {
        candidates.push_back({std::abs(mz - b.peaks[j].mz), static_cast<uint32_t>(i),
                              static_cast<uint32_t>(j)});
      }
    }
  }
  if (candidates.empty()) return 0.0f;
  std::sort(candidates.begin(), candidates.end(), [&](const Candidate& x, const Candidate& y) {
    if (x.difference != y.difference) return x.difference < y.difference;
    if (x.first != y.first) return a.peaks[x.first].original_index < a.peaks[y.first].original_index;
    return b.peaks[x.second].original_index < b.peaks[y.second].original_index;
  });

  used_a.assign(a.peaks.size(), 0);
  used_b.assign(b.peaks.size(), 0);
  double dot = 0.0;
  double matched_norm_a = 0.0;
  for (const auto& candidate : candidates) {
    if (used_a[candidate.first] || used_b[candidate.second]) continue;
    used_a[candidate.first] = used_b[candidate.second] = 1;
    double intensity_a = a.peaks[candidate.first].intensity;
    dot += intensity_a * b.peaks[candidate.second].intensity;
    matched_norm_a += intensity_a * intensity_a;
  }
  double norm_a = method == Method::ReverseCosine ? matched_norm_a : a.norm;
  double denominator = std::sqrt(norm_a * b.norm);
  double value = denominator > 0 ? dot / denominator : 0.0;
  if (!(value >= 0.0)) value = std::isnan(value) ? value : 0.0;  // keep NaN, clip negatives
  if (value > 1.0) value = 1.0;
  return static_cast<float>(value);
}

std::vector<PreparedSpectrum> prepare_all(const MSDataset& dataset, const Options& options) {
  std::vector<PreparedSpectrum> prepared(static_cast<std::size_t>(dataset.size()));
  for (int64_t i = 0; i < dataset.size(); ++i) {
    prepared[static_cast<std::size_t>(i)] = prepare_spectrum(dataset.peaks().spectrum(i), options);
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

std::shared_ptr<arrow::Array> float_array(const std::vector<float>& values) {
  arrow::FloatBuilder builder;
  internal::check(builder.AppendValues(values));
  return internal::unwrap(builder.Finish());
}

}  // namespace

std::string method_name(Method method) {
  return method == Method::Cosine ? "cosine" : "reverse_cosine";
}

Method parse_method(const std::string& name) {
  if (name == "cosine") return Method::Cosine;
  if (name == "reverse_cosine") return Method::ReverseCosine;
  throw std::invalid_argument("method must be 'cosine' or 'reverse_cosine'");
}

std::vector<float> pair_scores(const MSDataset& ds1, const std::vector<int64_t>& index1,
                               const MSDataset& ds2, const std::vector<int64_t>& index2,
                               const Options& options) {
  validate(options);
  if (index1.size() != index2.size()) {
    throw std::invalid_argument("index1 and index2 must have the same length");
  }
  check_positions(index1, ds1.size(), "index1");
  check_positions(index2, ds2.size(), "index2");

  // Prepare only the spectra that are used.
  auto prepare_used = [&](const MSDataset& dataset, const std::vector<int64_t>& index) {
    std::unordered_map<int64_t, PreparedSpectrum> prepared;
    for (int64_t position : index) {
      if (prepared.count(position) == 0) {
        prepared.emplace(position, prepare_spectrum(dataset.peaks().spectrum(position), options));
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
      scores[u] = score(prepared1.at(index1[u]), prepared2.at(index2[u]), tolerance, options.method);
    }
  });
  return scores;
}

std::vector<float> all_pairs_matrix(const MSDataset& ds1, const MSDataset& ds2,
                                    const Options& options) {
  validate(options);
  auto prepared1 = prepare_all(ds1, options);
  auto prepared2 = prepare_all(ds2, options);
  const Tolerance tolerance(options);
  const int64_t n1 = ds1.size(), n2 = ds2.size();
  std::vector<float> matrix(static_cast<std::size_t>(n1 * n2));
  parallel_ranges(n1, thread_count(options, n1 * n2), [&](int64_t begin, int64_t end) {
    for (int64_t i = begin; i < end; ++i) {
      for (int64_t j = 0; j < n2; ++j) {
        matrix[static_cast<std::size_t>(i * n2 + j)] =
            score(prepared1[static_cast<std::size_t>(i)], prepared2[static_cast<std::size_t>(j)],
                  tolerance, options.method);
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
  if (!ds1.has_column(key1)) throw KeyError("ds1 does not contain key column: " + key1);
  if (!ds2.has_column(key2)) throw KeyError("ds2 does not contain key column: " + key2);
  for (const auto& key : {key1, key2}) {
    if (key == "index1" || key == "index2" || key == "cosine_similarity") {
      throw std::invalid_argument(
          "Key column names must not collide with result columns: index1, index2, "
          "cosine_similarity");
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
  auto scores = pair_scores(ds1, index1, ds2, index2, options);

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
  add("cosine_similarity", float_array(scores));
  return arrow::Table::Make(arrow::schema(fields), arrays, static_cast<int64_t>(index1.size()));
}

std::shared_ptr<arrow::Table> library_search(const MSDataset& query, const MSDataset& reference,
                                             double threshold, const Options& options,
                                             const ProgressCallback& progress) {
  validate(options);
  if (!std::isfinite(threshold) || threshold < 0 || threshold > 1) {
    throw std::invalid_argument("threshold must be a finite value between 0 and 1");
  }
  const int64_t n_query = query.size(), n_reference = reference.size();
  const int64_t total = n_query * n_reference;
  auto prepared_query = prepare_all(query, options);
  auto prepared_reference = prepare_all(reference, options);
  const Tolerance tolerance(options);

  // Hits of each query row, filled in parallel, joined in order.
  std::vector<std::vector<std::pair<int64_t, float>>> hits(static_cast<std::size_t>(n_query));
  const unsigned threads = thread_count(options, total);
  const int64_t batch = std::max<int64_t>(1, 4096 / std::max<int64_t>(1, n_reference) + 1) * threads;
  if (progress) progress(0, total);
  for (int64_t start = 0; start < n_query; start += batch) {
    int64_t stop = std::min(n_query, start + batch);
    parallel_ranges(stop - start, threads, [&](int64_t begin, int64_t end) {
      for (int64_t q = start + begin; q < start + end; ++q) {
        auto& row = hits[static_cast<std::size_t>(q)];
        for (int64_t r = 0; r < n_reference; ++r) {
          float s = score(prepared_query[static_cast<std::size_t>(q)],
                          prepared_reference[static_cast<std::size_t>(r)], tolerance,
                          options.method);
          if (s >= threshold) row.emplace_back(r, s);
        }
      }
    });
    if (progress) progress(stop * n_reference, total);
  }

  std::vector<int64_t> index1, index2;
  std::vector<float> scores;
  for (int64_t q = 0; q < n_query; ++q) {
    for (const auto& [r, s] : hits[static_cast<std::size_t>(q)]) {
      index1.push_back(q);
      index2.push_back(r);
      scores.push_back(s);
    }
  }
  return arrow::Table::Make(
      arrow::schema({arrow::field("index1", arrow::int64()), arrow::field("index2", arrow::int64()),
                     arrow::field("cosine_similarity", arrow::float32())}),
      {int64_array(index1), int64_array(index2), float_array(scores)},
      static_cast<int64_t>(scores.size()));
}

}  // namespace msentity::similarity
