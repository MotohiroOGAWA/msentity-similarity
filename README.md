# msentity-similarity

Spectrum similarity for [msentity](https://github.com/MotohiroOGAWA/msentity)
datasets, and the `.mssim` file format for the results. Like msentity, it is a
C++ library with a C API, used from Python, Java, and C# through the same
shared library.

| | |
|---|---|
| Scores | cosine, reverse cosine, modified cosine, and BONANZA (one column each, several at once) with one-to-one peak matching within an m/z tolerance in Da and/or ppm (intensity exponent), multithreaded; per-side peak intensity coverage (upper bound of the cosine) |
| Modes | chosen pairs, all-pairs matrix, library search with threshold, pairs by metadata key |
| Results | `SimilarityDataset`: Arrow table + calculation metadata + optionally the compared spectra; filter, sort, statistics |
| Format | `.mssim` (HDF5 + Parquet + embedded `.msds`) — [specification](docs/format/mssim.md) |
| Languages | C++, C, Python (`msentity_similarity`), Java (`io.github.motohiroogawa.msentity.similarity`), C# (`MsEntity.Similarity`) |
| Depends on | msentity (fetched by CMake), C++20, CMake ≥ 3.21 |

## Installation

Ready-made packages for Linux x64 and Windows x64 are published with each
[GitHub Release](https://github.com/MotohiroOGAWA/msentity-similarity/releases).
The commands below use msentity-similarity 2.0.0, the latest release when this
README was written. To use another version, replace the msentity-similarity version `2.0.0`
(the newest is marked *Latest* on the [releases page](https://github.com/MotohiroOGAWA/msentity-similarity/releases/latest)).
Each release depends on one msentity release, named in its release notes;
for these examples, msentity 1.2.2.

| Platform | Python wheel tag | Java / C# / C platform |
|---|---|---|
| Linux x64 (glibc ≥ 2.35) | `manylinux_2_35_x86_64` | `linux-x64` |
| Windows x64 | `win_amd64` | `win-x64` |

| Language | Install |
|---|---|
| C | Download `msentity-similarity-c-<version>-<platform>.zip` from the release (it includes the msentity C library) |
| C++ | CMake `FetchContent` with the release tag; the compatible msentity source is fetched automatically |
| Python | `pip` with the release pages of msentity-similarity and msentity as `--find-links` (below) |
| Java | Download the msentity and msentity-similarity JARs from their releases |
| C# | Download both NuGet packages from their releases and use the folder as a package source |

### C

Download the archive for the current platform from the
[release page](https://github.com/MotohiroOGAWA/msentity-similarity/releases/tag/v2.0.0).
It contains the similarity and msentity C headers and native libraries.

### C++

```cmake
include(FetchContent)
FetchContent_Declare(msentity_similarity
  GIT_REPOSITORY https://github.com/MotohiroOGAWA/msentity-similarity.git
  GIT_TAG v2.0.0)  # msentity-similarity release
FetchContent_MakeAvailable(msentity_similarity)
target_link_libraries(app PRIVATE msentity::similarity)
```

### Python (≥ 3.10)

Give pip both release pages: the msentity-similarity wheel and the msentity
wheel it depends on:

```console
pip install msentity-similarity==2.0.0 \
  --find-links https://github.com/MotohiroOGAWA/msentity-similarity/releases/expanded_assets/v2.0.0 \
  --find-links https://github.com/MotohiroOGAWA/msentity/releases/expanded_assets/v1.2.2
```

The second URL is the msentity release that this similarity release depends on.

### Java (≥ 22)

Download both JARs with GitHub CLI:

```console
mkdir -p lib
gh release download v1.2.2 --repo MotohiroOGAWA/msentity \
  --pattern 'msentity-*.jar' --dir lib
gh release download v2.0.0 --repo MotohiroOGAWA/msentity-similarity \
  --pattern 'msentity-similarity-*.jar' --dir lib
java --enable-native-access=ALL-UNNAMED -cp 'app.jar:lib/*' app.Main
```

The JARs contain the native libraries for all supported platforms.

### C# (.NET ≥ 8)

Download both packages into one local source, then install similarity:

```console
mkdir -p packages
gh release download v1.2.2 --repo MotohiroOGAWA/msentity \
  --pattern 'MsEntity.*.nupkg' --dir packages
gh release download v2.0.0 --repo MotohiroOGAWA/msentity-similarity \
  --pattern 'MsEntity.Similarity.*.nupkg' --dir packages
dotnet add package MsEntity.Similarity --version 2.0.0 --source ./packages
```

The NuGet packages contain native libraries under `runtimes/<rid>/native/`.

More details and source-build instructions: [docs/install.md](docs/install.md).

## Example

```python
import msentity as ms
import msentity_similarity as mss

ds = ms.load("example.msp")                       # 2 spectra
hits = mss.library_search(ds, ds, threshold=0.7)   # methods="cosine" by default
hits.to_pandas()
#    index1  index2  cosine      peak_intensity_coverage1  peak_intensity_coverage2  data_index1  data_index2
# 0       0       0    1.000000             1.000000             1.000000            0            0
# 1       0       1    0.759714             0.929270             1.000000            0            1
# 2       1       0    0.759714             1.000000             0.929270            1            0
# 3       1       1    1.000000             1.000000             1.000000            1            1
query, reference = hits.match_records(1)          # the two compared spectra
hits.save("hits.mssim")
```

The same in C++:

```cpp
#include <msentity/msentity.hpp>
#include <msentity/similarity/similarity.hpp>

auto ds = msentity::io::load_ms_dataset("example.msp");
auto hits = msentity::similarity::calculate_library_search(ds, ds, 0.7);
hits.save("hits.mssim");
```

and on the command line:

```console
$ mssim-tool search example.msp example.msp hits.mssim --threshold 0.7
4 pairs written to hits.mssim
```

## Building from source

```console
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure    # C++, C API, and binding tests
```

See [docs/install.md](docs/install.md) for every language.

## Documentation

* [Installation](docs/install.md)
* [How scores are calculated](docs/similarity.md)
* [C++ API and mssim-tool](docs/api.md)
* [Python, Java, C#, and C](docs/bindings.md)
* [MSSIM file format](docs/format/mssim.md)

## Relation to msentity

msentity-similarity depends on msentity, never the other way round. Datasets
cross between the two C libraries as `.msds` bytes, so each can be built,
packaged, and versioned on its own.

## License

MIT — see [LICENSE](LICENSE).
