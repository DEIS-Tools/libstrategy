# libstrategy

Library for learning strategies.

## Dependencies

Library only:
- ptrie (will fetch if not installed)
- nlohmann_json (will fetch if not installed)

Full project:
- Boost: program_options tests

Ubuntu 24.04:
```shell
sudo apt install libboost-program-options-dev libboost-tests-dev
```

macOS:
```shell
brew install boost
```

## Compile and Install
```shell
```


## Compile for Development

```shell
cmake -B build
cmake --build build
ctest --test-dir build
```