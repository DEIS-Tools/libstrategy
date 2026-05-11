# libstrategy

Library for learning strategies.

## Dependencies

Library-only build requires:
- [ptrie](https://github.com/petergjoel/ptrie) (CMake will fetch automatically if not installed)
- [nlohmann_json](https://github.com/nlohmann/json) (CMake will fetch automatically if not installed)

Full project requires:
- [Boost](https://boost.org): [program_options](https://www.boost.org/library/latest/program_options/), [test](https://www.boost.org/library/latest/test/)

For Ubuntu 24.04 install build tools and library dependencies:
```shell
sudo apt install cmake ninja-build g++ 
sudo apt install libboost-program-options-dev libboost-test-dev  # optional
```

For macOS install build tools and library dependencies:
```shell
brew install cmake ninja gcc
brew install boost  # optional
```

## Compile and Install
Run minimal compilation (just the library) with `Debug` and `Release` settings into `build-quick/lib`:
```shell
cmake --workflow --preset quick-release
```
Install the `Release` build of `build-libonly` into `$PWD/local` path:
```shell
cmake --install build-quick --config Release --prefix $PWD/local
```

## Other Presets
Inspect workflow presets:
```shell
cmake --workflow --list-presets
```

For example, configure, build and **test** for Development with **Sanitizers** (GCC/Clang/AppleClang):
```shell
cmake --workflow --preset debug-san
```

Other configuration presets:
```shell
cmake --list-presets=configure
```

Other build presets:
```shell
cmake --list-presets=build
```

Other test presets:
```shell
cmake --list-presets=test
```

