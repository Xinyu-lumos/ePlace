# Repository Guidelines

## Project Structure & Module Organization
`ePlace` is a CMake-based C++ placement engine. Core modules live in top-level directories: `Parser/` for Bookshelf parsing, `PlaceDB/` and `PlaceCommon/` for shared data structures and globals, `QPlace/`, `EPlace/`, `Optimization/`, `Legalization/`, and `DetailedPlacement/` for placement stages, `FFT/` and `Plot/` for numerical and visualization helpers, and `main/` for the primary CLI entry point. Benchmarks and helper scripts live in `benchmarks/`. Third-party code is vendored under `Library/`; treat it as external unless you are intentionally updating a dependency.

## Build, Test, and Development Commands
Configure once from the repository root:

```bash
cmake -S . -B build
cmake --build build -j
```

Build specific targets when iterating:

```bash
cmake --build build --target ePlace
cmake --build build --target ParserTest
```

Run the main binary from `build/main/ePlace`:

```bash
./build/main/ePlace -aux benchmarks/ispd2005/adaptec1/adaptec1.aux -targetDensity 1.0
```

`ParserTest`, `QPPlacerTest`, `ePlaceTest`, and `ePlaceAbacus` are standalone smoke-test executables generated from the module directories.

## Coding Style & Naming Conventions
Match the existing codebase: 4-space indentation, braces on their own lines for functions and control blocks, `.h` headers with `.cpp` implementations, and PascalCase for types such as `BookshelfParser` or `PlaceDB`. Existing code uses `using namespace std;`, macros in `PlaceCommon/global.h`, and some global state such as `gArg`; avoid unrelated style cleanups in feature patches. Name new test drivers `<module>_test.cpp`.

## Testing Guidelines
There is no checked-in `ctest` suite or formatter target. Validate changes by building the affected executable and running the relevant smoke test directly, for example `./build/Parser/ParserTest -aux <case>.aux`. For algorithm changes, record the benchmark used, whether placement/legalization completed, and any HPWL or runtime deltas.

## Commit & Pull Request Guidelines
Recent history uses short, plain-English subjects such as `updated README` and `db info bug fixed`. Keep commit titles brief, imperative, and scoped to one change. Pull requests should describe the affected stage of the flow, list the build and benchmark commands you ran, and include output artifacts or screenshots when a change affects plotting, placement quality, or legalization behavior.
