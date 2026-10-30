# msentity-similarity

Spectrum similarity for [msentity](https://github.com/MotohiroOGAWA/msentity)
datasets, and the `.mssim` file format for the results. Like msentity, it is a
C++ library with a C API, used from Python, Java, and C# through the same
shared library.

| | |
|---|---|
| Scores | binned cosine and reverse cosine (bin width, intensity exponent), multithreaded |
| Modes | chosen pairs, all-pairs matrix, library search with threshold, pairs by metadata key |
| Results | `SimilarityDataset`: Arrow table + calculation metadata + optionally the compared spectra; filter, sort, statistics |
| Format | `.mssim` (HDF5 + Parquet + embedded `.msds`) — [specification](docs/format/mssim.md) |
| Languages | C++, C, Python (`msentity_similarity`), Java (`io.github.motohiroogawa.msentity.similarity`), C# (`MsEntity.Similarity`) |
| Depends on | msentity (fetched by CMake), C++20, CMake ≥ 3.21 |

## Installation

Ready-made packages for Linux x64 and Windows x64 are published
with each [GitHub Release](https://github.com/MotohiroOGAWA/msentity-similarity/releases).
`msentity-similarity` uses one pinned, compatible `msentity` release.

| Platform | Python wheel tag | Java / C# / C platform |
|---|---|---|
| Linux x64 (glibc ≥ 2.35) | `manylinux_2_35_x86_64` | `linux-x64` |
| Windows x64 | `win_amd64` | `win-x64` |

| Language | Install |
|---|---|
| C | Download `msentity-similarity-c-<version>-<platform>.zip` from the latest GitHub Release |
| C++ | Use CMake `FetchContent` with a release tag; the compatible msentity source is fetched automatically |
| Python | `pip install msentity-similarity --find-links https://motohiroogawa.github.io/msentity-similarity/` |
| Java | Download both the msentity and msentity-similarity JARs from their GitHub Releases |
| C# | Download both NuGet packages from their GitHub Releases and use the download directory as a package source |

### C

Download the archive for the current platform from the
[latest release](https://github.com/MotohiroOGAWA/msentity-similarity/releases/latest).
It contains the similarity and msentity C headers and native libraries.

### C++

Replace `<tag>` with the release to use, for example `v0.1.0`:

```cmake
include(FetchContent)
FetchContent_Declare(msentity_similarity
  GIT_REPOSITORY https://github.com/MotohiroOGAWA/msentity-similarity.git
  GIT_TAG <tag>)
FetchContent_MakeAvailable(msentity_similarity)
target_link_libraries(app PRIVATE msentity::similarity)
```

### Python (≥ 3.10)

Without a version, pip installs the latest release listed on the GitHub Pages
wheel index. The compatible msentity wheel is installed automatically:

```console
pip install msentity-similarity \
  --find-links https://motohiroogawa.github.io/msentity-similarity/
```

Pin a version for reproducible environments:

```console
pip install msentity-similarity==<version> \
  --find-links https://motohiroogawa.github.io/msentity-similarity/
```

### Java (≥ 22)

Download the latest msentity and similarity JARs with GitHub CLI:

```console
mkdir -p lib
gh release download --repo MotohiroOGAWA/msentity \
  --pattern 'msentity-*.jar' --dir lib
gh release download --repo MotohiroOGAWA/msentity-similarity \
  --pattern 'msentity-similarity-*.jar' --dir lib
java --enable-native-access=ALL-UNNAMED -cp 'app.jar:lib/*' app.Main
```

Pass a tag immediately after `gh release download`, such as
`gh release download v0.1.0 ...`, to select a fixed version. The JARs contain
the native libraries for all supported platforms.

### C# (.NET ≥ 8)

Download both packages into one local source, then install similarity:

```console
mkdir -p packages
gh release download --repo MotohiroOGAWA/msentity \
  --pattern 'MsEntity.*.nupkg' --dir packages
gh release download --repo MotohiroOGAWA/msentity-similarity \
  --pattern 'MsEntity.Similarity.*.nupkg' --dir packages
dotnet add package MsEntity.Similarity --source ./packages
```

Add `--version <version>` to the final command for a fixed release. The NuGet
packages contain native libraries under `runtimes/<rid>/native/`.

More details and source-build instructions: [docs/install.md](docs/install.md).

## Example

```python
import msentity as ms
import msentity_similarity as mss

ds = ms.load("example.msp")                       # 2 spectra
hits = mss.library_search(ds, ds, threshold=0.7)
hits.to_pandas()
#    index1  index2  cosine_similarity  data_index1  data_index2
# 0       0       0           1.000000            0            0
# 1       0       1           0.759714            0            1
# 2       1       0           0.759714            1            0
# 3       1       1           1.000000            1            1
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
