#include "msentity/similarity/similarity_dataset.hpp"

#include "msentity/similarity/calculation.hpp"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <map>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <unordered_map>

#include "internal/hdf5.hpp"
#include "internal/parquet.hpp"
#include "internal/table.hpp"
#include "internal/text.hpp"

namespace msentity::similarity {

namespace h5 = msentity::internal::h5;
namespace fs = std::filesystem;
using msentity::internal::check;
using msentity::internal::column_index;
using msentity::internal::unwrap;

namespace {

const char* kFormat = "msentity.similarity";
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

std::vector<int64_t> int_column(const arrow::Table& table, const std::string& name) {
  int index = column_index(table, name);
  if (index < 0) throw std::invalid_argument("missing column: " + name);
  std::vector<int64_t> values;
  for (const auto& value : internal::to_values(*table.column(index))) {
    if (const auto* i = std::get_if<int64_t>(&value)) {
      values.push_back(*i);
    } else if (const auto* d = std::get_if<double>(&value);
               d != nullptr && std::isfinite(*d) && *d == std::floor(*d)) {
      values.push_back(static_cast<int64_t>(*d));
    } else {
      throw std::invalid_argument("'" + name + "' must contain integers");
    }
  }
  return values;
}

// Values of a score column; missing values are NaN.
std::vector<double> score_values(const arrow::Table& table, const std::string& column) {
  std::vector<double> scores;
  for (const auto& value : internal::to_values(*table.column(column_index(table, column)))) {
    double score = kNaN;
    if (const auto* d = std::get_if<double>(&value)) score = *d;
    if (const auto* i = std::get_if<int64_t>(&value)) score = static_cast<double>(*i);
    if (const auto* b = std::get_if<bool>(&value)) score = *b ? 1.0 : 0.0;
    scores.push_back(score);
  }
  return scores;
}

// Non-missing values of a score column.
std::vector<double> present_scores(const arrow::Table& table, const std::string& column) {
  std::vector<double> scores;
  for (double score : score_values(table, column)) {
    if (!std::isnan(score)) scores.push_back(score);
  }
  return scores;
}

std::vector<std::string> score_column_names(const arrow::Table& table) {
  std::vector<std::string> names;
  for (Method method : all_methods()) {
    if (column_index(table, method_name(method)) >= 0) names.push_back(method_name(method));
  }
  return names;
}

nlohmann::ordered_json parse_object(const std::string& text) {
  auto json = nlohmann::ordered_json::parse(text, nullptr, false);
  if (json.is_discarded() || !json.is_object()) {
    throw std::invalid_argument("metadata must be a JSON object");
  }
  return json;
}

// A legacy table (one `similarity` column, no method columns): the column
// renamed to the method of the metadata, which gets parameters.methods.
void upgrade_legacy(std::shared_ptr<arrow::Table>& table, std::string& metadata_json) {
  if (!table || column_index(*table, "similarity") < 0 || !score_column_names(*table).empty()) return;
  auto json = nlohmann::ordered_json::parse(metadata_json, nullptr, false);
  std::string method = "cosine";
  if (json.is_object() && json.contains("parameters") && json["parameters"].is_object()) {
    auto& parameters = json["parameters"];
    if (parameters.contains("method") && parameters["method"].is_string()) {
      try {
        method = method_name(parse_method(parameters["method"].get<std::string>()));
      } catch (const std::invalid_argument&) {
        // keep cosine
      }
    }
    parameters["methods"] = nlohmann::ordered_json::array({method});
    metadata_json = json.dump();
  }
  auto names = internal::column_names(*table);
  for (auto& name : names) {
    if (name == "similarity") name = method;
  }
  table = unwrap(table->RenameColumns(names));
}

std::shared_ptr<arrow::Array> int64_array(const std::vector<int64_t>& values) {
  arrow::Int64Builder builder;
  check(builder.AppendValues(values));
  return unwrap(builder.Finish());
}

// Linear interpolation between order statistics (pandas/NumPy default).
double quantile(const std::vector<double>& sorted, double q) {
  if (sorted.empty()) return kNaN;
  double position = q * static_cast<double>(sorted.size() - 1);
  auto lower = static_cast<std::size_t>(std::floor(position));
  auto upper = static_cast<std::size_t>(std::ceil(position));
  double fraction = position - static_cast<double>(lower);
  return sorted[lower] + (sorted[upper] - sorted[lower]) * fraction;
}

// Ordering for sort_by: numbers by value, text by code point, missing last.
int compare(const Value& a, const Value& b) {
  auto numeric = [](const Value& v, double& out) {
    if (const auto* x = std::get_if<double>(&v)) { out = *x; return true; }
    if (const auto* x = std::get_if<int64_t>(&v)) { out = static_cast<double>(*x); return true; }
    if (const auto* x = std::get_if<bool>(&v)) { out = *x ? 1.0 : 0.0; return true; }
    return false;
  };
  double x = 0, y = 0;
  if (numeric(a, x) && numeric(b, y)) return x < y ? -1 : (x > y ? 1 : 0);
  const auto* sa = std::get_if<std::string>(&a);
  const auto* sb = std::get_if<std::string>(&b);
  if (sa && sb) return sa->compare(*sb) < 0 ? -1 : (*sa == *sb ? 0 : 1);
  throw std::invalid_argument("cannot sort a column mixing numbers and strings");
}

bool missing(const Value& value) {
  if (is_null(value)) return true;
  const auto* d = std::get_if<double>(&value);
  return d != nullptr && std::isnan(*d);
}

std::string csv_field(const std::string& value, char delimiter, bool peak_cell = false) {
  const bool peak_quotes = peak_cell && !value.empty() &&
      (value.find(';') != std::string::npos || internal::is_space(value.front()) || internal::is_space(value.back()));
  if (!peak_quotes && value.find_first_of(std::string{delimiter, '"', '\n', '\r'}) == std::string::npos) {
    return value;
  }
  return "\"" + internal::replace_all(value, "\"", "\"\"") + "\"";
}

std::string delimited_peaks(const Spectrum& spectrum) {
  std::vector<std::string> names;
  std::vector<std::vector<std::string>> annotations;
  if (auto metadata = spectrum.metadata()) {
    for (const auto& name : spectrum.metadata_columns()) {
      if (name == "mz" || name == "intensity") continue;
      auto column = metadata->GetColumnByName(name);
      if (!column) continue;
      std::vector<std::string> values;
      for (const auto& value : internal::to_values(*column)) values.push_back(to_display_string(value));
      if (std::none_of(values.begin(), values.end(), [](const auto& value) { return !value.empty(); })) continue;
      names.push_back(name);
      annotations.push_back(std::move(values));
    }
  }
  std::string peaks;
  if (!names.empty()) {
    peaks = "mz,intensity";
    for (const auto& name : names) peaks += "," + csv_field(name, ',', true);
  }
  for (int64_t peak = 0; peak < spectrum.size(); ++peak) {
    if (!peaks.empty()) peaks += ';';
    peaks += internal::format_g17(spectrum.mz(peak)) + "," + internal::format_g17(spectrum.intensity(peak));
    auto index = static_cast<std::size_t>(peak);
    std::size_t last = annotations.size();
    while (last > 0 && annotations[last - 1][index].empty()) --last;
    for (std::size_t c = 0; c < last; ++c) peaks += "," + csv_field(annotations[c][index], ',', true);
  }
  return peaks;
}

}  // namespace

// ---------------------------------------------------------------------------

SimilarityDataset::SimilarityDataset(std::shared_ptr<arrow::Table> table, std::string metadata_json,
                                     std::optional<DatasetPair> matched_datasets,
                                     std::optional<IndexPair> matched_source_indices)
    : table_(table ? internal::normalize_table(table) : nullptr),
      metadata_json_(std::move(metadata_json)),
      matched_datasets_(std::move(matched_datasets)),
      source_indices_(std::move(matched_source_indices)) {
  upgrade_legacy(table_, metadata_json_);
  validate();
}

void SimilarityDataset::validate() const {
  if (!table_) throw std::invalid_argument("table must not be null");
  parse_object(metadata_json_);
  std::vector<std::string> missing_columns;
  for (const char* column : {"index1", "index2"}) {
    if (column_index(*table_, column) < 0) missing_columns.push_back(column);
  }
  if (!missing_columns.empty()) {
    std::string list;
    for (const auto& m : missing_columns) list += (list.empty() ? "'" : ", '") + m + "'";
    throw std::invalid_argument("Similarity table is missing required columns: [" + list + "]");
  }
  auto scores = score_column_names(*table_);
  if (scores.empty()) {
    throw std::invalid_argument(
        "Similarity table needs a score column: cosine, reverse_cosine, modified_cosine, bonanza, "
        "entropy, or unweighted_entropy");
  }
  for (const auto& column : scores) {
    auto type = table_->column(column_index(*table_, column))->type()->id();
    if (!(arrow::is_integer(type) || arrow::is_floating(type) || type == arrow::Type::BOOL)) {
      throw std::invalid_argument("Similarity scores must be numeric: " + column);
    }
    for (double score : present_scores(*table_, column)) {
      if (!std::isfinite(score) || score < 0 || score > 1) {
        throw std::invalid_argument("Similarity scores must be missing or finite values between 0 and 1: " +
                                    column);
      }
    }
  }
  if (matched_datasets_.has_value() != source_indices_.has_value()) {
    throw std::invalid_argument("matched datasets and source indices must be provided together");
  }
  if (!matched_datasets_) return;
  for (int side = 1; side <= 2; ++side) {
    const MSDataset& dataset = side == 1 ? matched_datasets_->first : matched_datasets_->second;
    const auto& sources = side == 1 ? source_indices_->first : source_indices_->second;
    std::string column = "data_index" + std::to_string(side);
    if (column_index(*table_, column) < 0) {
      throw std::invalid_argument("Embedded matches require the '" + column + "' column");
    }
    if (static_cast<int64_t>(sources.size()) != dataset.size()) {
      throw std::invalid_argument("Embedded source indices must match their datasets");
    }
    for (int64_t index : int_column(*table_, column)) {
      if (index < 0 || index >= static_cast<int64_t>(sources.size())) {
        throw std::invalid_argument("'" + column + "' contains an invalid embedded index");
      }
    }
  }
}

SimilarityDataset SimilarityDataset::with_matched_data(std::shared_ptr<arrow::Table> table,
                                                       std::string metadata_json,
                                                       const MSDataset& first,
                                                       const MSDataset& second) {
  if (!table) throw std::invalid_argument("table must not be null");
  std::vector<MSDataset> datasets;
  std::vector<std::vector<int64_t>> sources;
  for (int side = 1; side <= 2; ++side) {
    const MSDataset& dataset = side == 1 ? first : second;
    auto indices = int_column(*table, "index" + std::to_string(side));
    std::unordered_map<int64_t, int64_t> lookup;
    std::vector<int64_t> unique, local;
    for (int64_t index : indices) {
      if (index < 0 || index >= dataset.size()) {
        throw std::invalid_argument("'index" + std::to_string(side) +
                                    "' contains an invalid source index");
      }
      auto [it, inserted] = lookup.emplace(index, static_cast<int64_t>(unique.size()));
      if (inserted) unique.push_back(index);
      local.push_back(it->second);
    }
    table = internal::put_column(table, "data_index" + std::to_string(side), int64_array(local));
    datasets.push_back(dataset.select(unique).copy());
    sources.push_back(std::move(unique));
  }
  return SimilarityDataset(table, std::move(metadata_json),
                           DatasetPair{std::move(datasets[0]), std::move(datasets[1])},
                           IndexPair{std::move(sources[0]), std::move(sources[1])});
}

int64_t SimilarityDataset::size() const { return table_->num_rows(); }

std::vector<std::string> SimilarityDataset::columns() const {
  return internal::column_names(*table_);
}

std::vector<std::string> SimilarityDataset::score_columns() const {
  return score_column_names(*table_);
}

namespace {

std::string score_column(const SimilarityDataset& result, const std::string& column) {
  auto names = result.score_columns();
  if (column.empty()) return names.front();
  if (std::find(names.begin(), names.end(), column) == names.end()) {
    throw KeyError("not a score column: " + column);
  }
  return column;
}

}  // namespace

std::vector<MatchedMetadataColumn> SimilarityDataset::matched_metadata_columns() const {
  if (!matched_datasets_) return {};
  auto used = columns();
  std::vector<MatchedMetadataColumn> fields;
  const MSDataset* datasets[] = {&matched_datasets_->first, &matched_datasets_->second};
  for (int side : {1, 2}) {
    const std::string prefix = side == 1 ? "Query." : "Reference.";
    for (const auto& name : datasets[side - 1]->columns()) {
      std::string key = prefix + name;
      while (std::find(used.begin(), used.end(), key) != used.end()) key = prefix + key;
      used.push_back(key);
      fields.push_back({side, name, std::move(key)});
    }
  }
  return fields;
}

std::vector<Value> SimilarityDataset::column_values(const std::string& name) const {
  int index = column_index(*table_, name);
  if (index >= 0) return internal::to_values(*table_->column(index));
  for (const auto& field : matched_metadata_columns()) {
    if (field.key != name) continue;
    const auto& dataset = field.side == 1 ? matched_datasets_->first : matched_datasets_->second;
    const auto values = dataset.column_values(field.column);
    std::vector<Value> aligned;
    aligned.reserve(static_cast<std::size_t>(size()));
    for (int64_t position : int_column(*table_, "data_index" + std::to_string(field.side))) {
      aligned.push_back(values.at(static_cast<std::size_t>(position)));
    }
    return aligned;
  }
  throw KeyError("column not found: " + name);
}

void SimilarityDataset::set_metadata_json(std::string json) {
  parse_object(json);
  metadata_json_ = std::move(json);
}

void SimilarityDataset::set_annotations(std::string description, std::vector<std::string> tags,
                                        Attributes attributes) {
  MSDataset annotations;
  annotations.set_tags(std::move(tags));
  annotations.set_attributes(std::move(attributes));
  auto metadata = parse_object(metadata_json_);
  auto fields = nlohmann::ordered_json::object();
  for (const auto& [name, value] : annotations.attributes()) fields[name] = value;
  metadata["annotations"] = {{"description", std::move(description)}, {"tags", annotations.tags()},
                             {"attributes", std::move(fields)}};
  metadata_json_ = metadata.dump();
}

SimilarityDataset::DatasetPair SimilarityDataset::match_records(int64_t row) const {
  if (!matched_datasets_) {
    throw std::invalid_argument("This similarity result does not contain matched data");
  }
  if (row < 0 || row >= size()) throw std::out_of_range("row out of range");
  int64_t first = int_column(*table_, "data_index1")[static_cast<std::size_t>(row)];
  int64_t second = int_column(*table_, "data_index2")[static_cast<std::size_t>(row)];
  return {matched_datasets_->first.select({first}), matched_datasets_->second.select({second})};
}

SimilarityDataset SimilarityDataset::filter(const std::vector<bool>& mask) const {
  if (static_cast<int64_t>(mask.size()) != size()) {
    throw std::invalid_argument("mask must have one entry per row");
  }
  std::vector<int64_t> rows;
  for (std::size_t i = 0; i < mask.size(); ++i) {
    if (mask[i]) rows.push_back(static_cast<int64_t>(i));
  }
  return select_rows(rows);
}

SimilarityDataset SimilarityDataset::select_rows(const std::vector<int64_t>& rows) const {
  auto json = parse_object(metadata_json_);
  json["row_count"] = rows.size();
  json["source_row_count"] = size();
  return SimilarityDataset(internal::take_rows(*table_, rows), json.dump(), matched_datasets_,
                           source_indices_);
}

SimilarityDataset SimilarityDataset::filter_by(const FilterGroup& group) const {
  return select_rows(row_positions(group));
}

std::vector<int64_t> SimilarityDataset::row_positions(const FilterGroup& filters,
                                                       const std::vector<std::string>& columns,
                                                       const std::vector<bool>& ascending,
                                                       const std::vector<SortMode>& modes) const {
  if (!ascending.empty() && ascending.size() != columns.size()) {
    throw std::invalid_argument("ascending must have one entry per column");
  }
  if (!modes.empty() && modes.size() != columns.size()) {
    throw std::invalid_argument("modes must have one entry per column");
  }
  std::vector<int64_t> order;
  if (filters.children.empty() && filters.logic == FilterLogic::And) {
    order.resize(static_cast<std::size_t>(size()));
    std::iota(order.begin(), order.end(), 0);
  } else {
    // FilterEvaluator only needs the columns used by the filter. Give it
    // one empty spectrum per result row and no embedded peak data.
    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;
    std::map<std::pair<std::string, bool>, std::string> projections;
    auto group = filters;
    auto prepare = [&](auto&& self, FilterGroup& current) -> void {
      for (auto& node : current.children) {
        if (auto* nested = std::get_if<FilterGroup>(&node)) {
          self(self, *nested);
          continue;
        }
        auto& condition = std::get<FilterCondition>(node);
        if (condition.op == FilterOperator::SpectrumPeak || condition.op == FilterOperator::NeutralLoss) {
          throw std::invalid_argument("Spectrum peak filters require a spectrum dataset");
        }
        if (condition.op == FilterOperator::Mask || condition.op == FilterOperator::RowIds) continue;
        const bool text = condition.op == FilterOperator::TextEquals ||
                          condition.op == FilterOperator::TextNotEquals || condition.op == FilterOperator::Contains;
        auto [entry, inserted] = projections.try_emplace({condition.column, text}, std::to_string(fields.size()));
        if (inserted) {
          auto values = column_values(condition.column);
          if (text) {
            for (auto& value : values) value = to_display_string(value);
          }
          auto array = internal::values_to_array(values);
          fields.push_back(arrow::field(entry->second, array->type()));
          arrays.push_back(std::move(array));
        }
        if (text) condition.value = to_display_string(condition.value);
        condition.column = entry->second;
      }
    };
    prepare(prepare, group);
    MSDataset projection(arrow::Table::Make(arrow::schema(fields), arrays, size()),
                         PeakSeries({}, std::vector<int64_t>(static_cast<std::size_t>(size()) + 1, 0)));
    const auto mask = FilterEvaluator::evaluate(projection, group);
    for (std::size_t i = 0; i < mask.size(); ++i) {
      if (mask[i]) order.push_back(static_cast<int64_t>(i));
    }
  }
  if (columns.empty()) return order;
  std::vector<std::vector<Value>> keys;
  for (std::size_t c = 0; c < columns.size(); ++c) {
    auto values = column_values(columns[c]);
    auto mode = modes.empty() ? SortMode::Native : modes[c];
    if (mode == SortMode::Text) {
      for (auto& value : values) {
        if (!missing(value)) value = to_display_string(value);
      }
    } else if (mode == SortMode::Numeric) {
      for (auto& value : values) {
        if (const auto* text = std::get_if<std::string>(&value)) {
          auto number = internal::parse_number(*text);
          value = number && std::isfinite(*number) ? Value(*number) : Value(std::monostate{});
        } else if (const auto* number = std::get_if<double>(&value); number && !std::isfinite(*number)) {
          value = std::monostate{};
        }
      }
    }
    keys.push_back(std::move(values));
  }
  std::stable_sort(order.begin(), order.end(), [&](int64_t a, int64_t b) {
    for (std::size_t c = 0; c < keys.size(); ++c) {
      const Value& va = keys[c][static_cast<std::size_t>(a)];
      const Value& vb = keys[c][static_cast<std::size_t>(b)];
      bool ma = missing(va), mb = missing(vb);
      if (ma || mb) {
        if (ma != mb) return mb;  // missing values last
        continue;
      }
      int cmp = compare(va, vb);
      if (cmp != 0) return (ascending.empty() || ascending[c]) ? cmp < 0 : cmp > 0;
    }
    return false;
  });
  return order;
}

SimilarityDataset SimilarityDataset::sort_by(const std::vector<std::string>& columns,
                                             const std::vector<bool>& ascending) const {
  if (columns.empty()) throw std::invalid_argument("columns must not be empty");
  return SimilarityDataset(internal::take_rows(*table_, row_positions({}, columns, ascending)),
                           metadata_json_, matched_datasets_, source_indices_);
}

ScoreSummary SimilarityDataset::describe_scores(const std::string& column) const {
  std::vector<double> scores = present_scores(*table_, score_column(*this, column));
  ScoreSummary summary;
  summary.count = static_cast<int64_t>(scores.size());
  if (scores.empty()) {
    summary.mean = summary.q1 = summary.median = summary.q3 = summary.min = summary.max = kNaN;
    return summary;
  }
  std::sort(scores.begin(), scores.end());
  summary.mean = std::accumulate(scores.begin(), scores.end(), 0.0) / scores.size();
  summary.q1 = quantile(scores, 0.25);
  summary.median = quantile(scores, 0.5);
  summary.q3 = quantile(scores, 0.75);
  summary.min = scores.front();
  summary.max = scores.back();
  return summary;
}

Histogram SimilarityDataset::histogram(int bins, const std::string& column) const {
  if (bins < 1 || bins > 200) throw std::invalid_argument("bins must be an integer between 1 and 200");
  const std::string name = score_column(*this, column);
  Histogram result;
  result.counts.assign(static_cast<std::size_t>(bins), 0);
  for (int i = 0; i <= bins; ++i) result.edges.push_back(static_cast<double>(i) / bins);
  for (double score : present_scores(*table_, name)) {
    auto bin = static_cast<int>(std::floor(score * bins));
    bin = std::clamp(bin, 0, bins - 1);  // 1.0 belongs to the last bin
    ++result.counts[static_cast<std::size_t>(bin)];
  }
  return result;
}

void SimilarityDataset::export_table(const std::string& path) const {
  std::string ext = internal::to_lower_ascii(fs::path(path).extension().string());
  fs::path target(path);
  if (target.has_parent_path()) fs::create_directories(target.parent_path());
  if (ext == ".parquet") {
    auto bytes = internal::table_to_parquet(*table_);
    std::ofstream(path, std::ios::binary).write(reinterpret_cast<const char*>(bytes.data()),
                                                static_cast<std::streamsize>(bytes.size()));
    return;
  }
  if (ext != ".csv" && ext != ".tsv") {
    throw std::invalid_argument("Table export format must be parquet, csv, or tsv");
  }
  export_table(path, DelimitedExportOptions{});
}

void SimilarityDataset::export_table(const std::string& path, const DelimitedExportOptions& options) const {
  fs::path target(path);
  std::string ext = internal::to_lower_ascii(target.extension().string());
  if (ext != ".csv" && ext != ".tsv") throw std::invalid_argument("Selected columns and peaks require CSV or TSV");
  const auto names = options.columns.value_or(columns());
  if (names.empty() && options.peak_sides.empty()) throw std::invalid_argument("Select at least one export column");
  std::vector<std::vector<Value>> values;
  for (const auto& name : names) values.push_back(column_values(name));
  auto headers = names;
  std::vector<std::vector<int64_t>> peak_indices;
  for (int side : options.peak_sides) {
    if ((side != 1 && side != 2) || !matched_datasets_) {
      throw std::invalid_argument("Spectrum peaks require embedded matched data");
    }
    std::string prefix = side == 1 ? "Query." : "Reference.";
    std::string name = prefix + "Peak";
    while (std::find(headers.begin(), headers.end(), name) != headers.end()) name = prefix + name;
    headers.push_back(std::move(name));
    peak_indices.push_back(int_column(*table_, "data_index" + std::to_string(side)));
  }
  if (target.has_parent_path()) fs::create_directories(target.parent_path());
  std::ofstream out(target, std::ios::binary);
  if (!out) throw std::runtime_error("cannot open file for writing: " + path);
  char delimiter = ext == ".csv" ? ',' : '\t';
  auto write_field = [&](std::size_t column, const std::string& value) {
    if (column > 0) out << delimiter;
    out << csv_field(value, delimiter);
  };
  for (std::size_t c = 0; c < headers.size(); ++c) write_field(c, headers[c]);
  out << '\n';
  for (int64_t row = 0; row < size(); ++row) {
    auto index = static_cast<std::size_t>(row);
    for (std::size_t c = 0; c < names.size(); ++c) write_field(c, to_display_string(values[c][index]));
    for (std::size_t c = 0; c < options.peak_sides.size(); ++c) {
      const auto& dataset = options.peak_sides[c] == 1 ? matched_datasets_->first : matched_datasets_->second;
      write_field(names.size() + c, delimited_peaks(dataset.peaks().spectrum(peak_indices[c][index])));
    }
    out << '\n';
  }
  out.close();
  if (!out) throw std::runtime_error("failed to write " + path);
}

// ---- files -----------------------------------------------------------------

SimilarityDataset SimilarityDataset::load(const std::string& path) {
  h5::Handle file = h5::open_file(path, h5::OpenMode::ReadOnly);
  if (h5::read_string_attribute(file.id(), "format").value_or("") != kFormat) {
    throw FormatError("Not an msentity similarity file: " + path);
  }
  int64_t version = h5::read_int_attribute(file.id(), "schema_version").value_or(0);
  if (version < 1 || version > 3) throw FormatError("Unsupported similarity schema version");

  auto payload = h5::read_bytes_dataset(file.id(), "table.parquet");
  auto table = internal::parquet_to_table(payload.data(), payload.size());
  auto metadata = h5::read_string_dataset(file.id(), "metadata.json");
  if (metadata.size() != 1) throw FormatError("metadata.json must be a scalar string");

  std::optional<DatasetPair> datasets;
  std::optional<IndexPair> sources;
  if (version >= 2 && h5::link_exists(file.id(), "matched_data")) {
    h5::Handle embedded = h5::open_group(file.id(), "matched_data");
    std::vector<MSDataset> loaded;
    std::vector<std::vector<int64_t>> indices;
    for (const char* side : {"1", "2"}) {
      h5::Handle group = h5::open_group(embedded.id(), side);
      loaded.push_back(io::msds_from_bytes(h5::read_bytes_dataset(group.id(), "dataset.msds")));
      indices.push_back(h5::read_int64_dataset(group.id(), "source_indices"));
    }
    datasets = DatasetPair{std::move(loaded[0]), std::move(loaded[1])};
    sources = IndexPair{std::move(indices[0]), std::move(indices[1])};
  }
  return SimilarityDataset(table, metadata.front(), std::move(datasets), std::move(sources));
}

void SimilarityDataset::save(const std::string& path) const {
  fs::path target(path);
  if (internal::to_lower_ascii(target.extension().string()) != ".mssim") {
    throw std::invalid_argument("Similarity files must use the .mssim extension");
  }
  if (target.has_parent_path()) fs::create_directories(target.parent_path());

  // Keep only the embedded spectra the table refers to, renumbering
  // data_index1/2 in order of first appearance.
  std::shared_ptr<arrow::Table> table = table_;
  std::vector<MSDataset> datasets;
  std::vector<std::vector<int64_t>> sources;
  if (matched_datasets_) {
    for (int side = 1; side <= 2; ++side) {
      std::string column = "data_index" + std::to_string(side);
      std::unordered_map<int64_t, int64_t> lookup;
      std::vector<int64_t> referenced, renumbered;
      for (int64_t old_index : int_column(*table, column)) {
        auto [it, inserted] = lookup.emplace(old_index, static_cast<int64_t>(referenced.size()));
        if (inserted) referenced.push_back(old_index);
        renumbered.push_back(it->second);
      }
      table = internal::put_column(table, column, int64_array(renumbered));
      const MSDataset& dataset = side == 1 ? matched_datasets_->first : matched_datasets_->second;
      const auto& source = side == 1 ? source_indices_->first : source_indices_->second;
      datasets.push_back(dataset.select(referenced).copy());
      std::vector<int64_t> kept;
      for (int64_t index : referenced) kept.push_back(source[static_cast<std::size_t>(index)]);
      sources.push_back(std::move(kept));
    }
  }

  auto payload = internal::table_to_parquet(*table);
  std::random_device random;
  fs::path temporary = target;
  temporary += ".tmp" + std::to_string(random());
  try {
    {
      h5::Handle file = h5::open_file(temporary.string(), h5::OpenMode::Truncate);
      h5::write_string_attribute(file.id(), "format", kFormat);
      h5::write_int_attribute(file.id(), "schema_version", 3);
      h5::write_string_attribute(file.id(), "data_mode",
                                 matched_datasets_ ? "embedded" : "lightweight");
      h5::write_bytes_dataset(file.id(), "table.parquet", payload.data(), payload.size(), false);
      h5::write_scalar_string_dataset(file.id(), "metadata.json", metadata_json_);
      if (matched_datasets_) {
        h5::Handle embedded = h5::create_group(file.id(), "matched_data");
        for (int side = 0; side < 2; ++side) {
          h5::Handle group = h5::create_group(embedded.id(), std::to_string(side + 1));
          h5::write_int64_dataset(group.id(), "source_indices", sources[static_cast<std::size_t>(side)]);
          auto bytes = io::msds_to_bytes(datasets[static_cast<std::size_t>(side)]);
          h5::write_bytes_dataset(group.id(), "dataset.msds", bytes.data(), bytes.size(), true);
        }
      }
    }
    fs::rename(temporary, target);
  } catch (...) {
    std::error_code ignored;
    fs::remove(temporary, ignored);
    throw;
  }
}

}  // namespace msentity::similarity
