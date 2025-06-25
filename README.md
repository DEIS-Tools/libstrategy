# libstrategy

Library for learning strategies.

## Dependencies

Library only requires:
- [ptrie](https://github.com/petergjoel/ptrie) (CMake will fetch automatically if not installed)
- [nlohmann_json](https://github.com/nlohmann/json) (CMake will fetch automatically if not installed)

Full project requires:
- [Boost](https://boost.org): [program_options](https://www.boost.org/library/latest/program_options/), [test](https://www.boost.org/library/latest/test/)

For Ubuntu 24.04 install build tools and library dependencies:
```shell
sudo apt install cmake ninja-build g++ libboost-program-options-dev libboost-tests-dev
```

For macOS install build tools and library dependencies:
```shell
brew install cmake ninja boost
```

## Compile and Install
Inspect common workflow presets:
```shell
cmake --workflow --list-presets
```

Run minimal preset to build the library with `Debug` and `Release` settings into `build-libonly/lib`:
```shell
cmake --workflow libonly
```
Install the `Release` build of `build-libonly` into `$PWD/local` path:
```shell
cmake --install build-libonly --config Release --prefix $PWD/local
```

Configure, build and test for Development with Sanitizers (GCC/Clang/AppleClang):
```shell
cmake --workflow debug-san
```

Configuration presets:
```shell
cmake --list-presets=configure
```

Build presets:
```shell
cmake --list-presets=build
```

Test presets:
```shell
cmake --list-presets=test
```

