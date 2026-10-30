# Installation

`msentity-similarity` uses the released `msentity` version `0.1.0`. The same
version is pinned in CMake and in every language binding so that a similarity
release is tested against one known data-model release.

## C

Build the shared C API and include
`msentity/similarity/c/msentity_similarity_c.h`:

```console
cmake -S . -B build -G Ninja -DMSENTITY_SIMILARITY_BUILD_C_API=ON
cmake --build build
```

Link `build/libmsentity_similarity_c.so` (`.dylib` or `.dll` on other
platforms). Datasets created with the msentity C API cross the boundary as
`.msds` bytes.

## C++

Use a released `msentity-similarity` tag. Its CMake configuration fetches the
compatible released `msentity` automatically:

```cmake
include(FetchContent)
FetchContent_Declare(msentity_similarity
  GIT_REPOSITORY https://github.com/MotohiroOGAWA/msentity-similarity.git
  GIT_TAG v0.1.0)
FetchContent_MakeAvailable(msentity_similarity)

target_link_libraries(app PRIVATE msentity::similarity)
```

The dependency brings `msentity::msentity` transitively. A separate msentity
checkout and `FETCHCONTENT_SOURCE_DIR_MSENTITY` are only needed when developing
both repositories together.

## Python (≥ 3.10)

Install the latest release from PyPI:

```console
pip install msentity-similarity
```

The package metadata installs the compatible `msentity==0.1.0` automatically.
Pin the similarity version for a reproducible environment:

```console
pip install msentity-similarity==0.1.0
```

```python
import msentity as ms
import msentity_similarity as mss

query = ms.load("query.msp")
reference = ms.load("reference.msp")
hits = mss.library_search(query, reference, threshold=0.7)
hits.save("hits.mssim")
```

## Java

The Java artifact uses `io.github.motohiroogawa:msentity:0.1.0`. Until the
similarity artifact is published to a Maven repository, build it locally:

```console
(cd bindings/java && mvn install -DskipTests)
```

## C#

The project references `MsEntity` version `0.1.0` when a sibling development
checkout is not present:

```console
dotnet pack bindings/csharp/MsEntity.Similarity -c Release -o packages \
  -p:UseMsEntityPackage=true
```

## Building and testing from source

```console
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

CMake downloads `msentity` tag `v0.1.0`. To test unpublished changes in both
repositories, override it explicitly:

```console
cmake -S . -B build -G Ninja \
  -DFETCHCONTENT_SOURCE_DIR_MSENTITY=../msentity
```

Release tags build self-contained wheels for Linux x64, Windows x64, and macOS
arm64 and publish them to PyPI using the
[release workflow](../.github/workflows/release.yml).
