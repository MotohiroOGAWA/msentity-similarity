# Installation

Each msentity-similarity release depends on one msentity release: CMake fetches
its tag and every language package requires that version. The examples use
msentity-similarity 1.2.1, the latest release when this page was written.
To use another version, replace the msentity-similarity version `1.2.1`.
The msentity release it depends on is named in its release notes and in
`CMakeLists.txt` (`GIT_TAG`); for these examples, msentity 1.2.0.

## C

Each release attaches `msentity-similarity-c-<version>-<platform>.zip` with
the similarity and msentity headers and libraries. To build the C API
yourself, include `msentity/similarity/c/msentity_similarity_c.h`:

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
  GIT_TAG v1.2.1)  # msentity-similarity release
FetchContent_MakeAvailable(msentity_similarity)

target_link_libraries(app PRIVATE msentity::similarity)
```

The dependency brings `msentity::msentity` transitively. A separate msentity
checkout and `FETCHCONTENT_SOURCE_DIR_MSENTITY` are only needed when developing
both repositories together.

## Python (≥ 3.10)

The wheels are attached to the GitHub Releases. Give pip the release page of
msentity-similarity and that of the msentity release it depends on:

```console
pip install msentity-similarity==1.2.1 \
  --find-links https://github.com/MotohiroOGAWA/msentity-similarity/releases/expanded_assets/v1.2.1 \
  --find-links https://github.com/MotohiroOGAWA/msentity/releases/expanded_assets/v1.2.0
```

(Earlier instructions used a GitHub Pages index; it is no longer updated.) <!-- fixed -->

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
mkdir -p lib
curl -L -o lib/msentity.jar \
  https://github.com/MotohiroOGAWA/msentity/releases/download/v1.2.0/msentity-1.2.0.jar
curl -L -o lib/msentity-similarity.jar \
  https://github.com/MotohiroOGAWA/msentity-similarity/releases/download/v1.2.1/msentity-similarity-1.2.1.jar
java --enable-native-access=ALL-UNNAMED -cp 'app.jar:lib/*' app.Main
```

## C#

Download both NuGet packages from their GitHub Releases into one directory,
then add that directory as the package source:

```console
mkdir -p packages
curl -L -o packages/MsEntity.1.2.0.nupkg \
  https://github.com/MotohiroOGAWA/msentity/releases/download/v1.2.0/MsEntity.1.2.0.nupkg
curl -L -o packages/MsEntity.Similarity.1.2.1.nupkg \
  https://github.com/MotohiroOGAWA/msentity-similarity/releases/download/v1.2.1/MsEntity.Similarity.1.2.1.nupkg
dotnet add package MsEntity.Similarity --version 1.2.1 --source ./packages
```

## Building and testing from source

```console
cmake -S . -B build -G Ninja
cmake --build build
ctest --test-dir build --output-on-failure
```

CMake downloads the msentity tag named in `CMakeLists.txt` (`GIT_TAG`). To test
unpublished changes in both repositories, override it explicitly:

```console
cmake -S . -B build -G Ninja \
  -DFETCHCONTENT_SOURCE_DIR_MSENTITY=../msentity
```

## How releases are made

Releases are built and tested on Ubuntu 22.04 (Linux x64, glibc 2.35) and on
Windows (x64), and every package is attached to the GitHub Release of the tag
`v<version>`. The maintainers' scripts are in
[msentity-spectrum-viewer/scripts/release-all](https://github.com/MotohiroOGAWA/msentity-spectrum-viewer/tree/main/scripts/release-all).
