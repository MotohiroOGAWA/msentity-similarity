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

Install the latest release from the GitHub Pages wheel index:

```console
pip install msentity-similarity \
  --find-links https://motohiroogawa.github.io/msentity-similarity/
```

The package metadata installs the compatible `msentity==0.1.0` automatically.
Pin the similarity version for a reproducible environment:

```console
pip install msentity-similarity==0.1.0 \
  --find-links https://motohiroogawa.github.io/msentity-similarity/
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

Download both JARs from the corresponding GitHub Releases. The similarity JAR
contains its native library for every supported platform; the msentity JAR does
the same for the dataset API:

```console
curl -LO https://github.com/MotohiroOGAWA/msentity/releases/download/v0.1.0/msentity-0.1.0.jar
curl -LO https://github.com/MotohiroOGAWA/msentity-similarity/releases/download/v0.1.0/msentity-similarity-0.1.0.jar
java --enable-native-access=ALL-UNNAMED \
  -cp 'app.jar:msentity-0.1.0.jar:msentity-similarity-0.1.0.jar' app.Main
```

## C#

Download both NuGet packages from their GitHub Releases into one directory,
then add that directory as the package source:

```console
curl -LO https://github.com/MotohiroOGAWA/msentity/releases/download/v0.1.0/MsEntity.0.1.0.nupkg
curl -LO https://github.com/MotohiroOGAWA/msentity-similarity/releases/download/v0.1.0/MsEntity.Similarity.0.1.0.nupkg
mkdir packages
mv MsEntity*.nupkg packages/
dotnet add package MsEntity.Similarity --version 0.1.0 --source ./packages
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

Release tags build self-contained wheels for Linux x64 and Windows x64, attach
them to a GitHub Release, and update the GitHub Pages pip index
using the [release workflow](../.github/workflows/release.yml). The similarity
index also links the compatible `msentity==0.1.0` wheels, so one command installs
both packages.
