cmake_minimum_required(VERSION 3.28)
if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.30.0")
	cmake_policy(SET CMP0167 NEW) # use FindBoost from boost installation (as opposed to CMake own)
endif ()
project(libstrategy VERSION 1.0.9 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

option(LIBSTRATEGY_TESTS "libstrategy Unit Tests" ON)
option(LIBSTRATEGY_OnlyLibrary "Build only the library (no binary utilities)" OFF)

include(cmake/warnings.cmake)
include(cmake/sanitizers.cmake)
include(cmake/nlohmann_json.cmake)
include(cmake/ptrie.cmake)

if(NOT LIBSTRATEGY_OnlyLibrary)
	cmake_policy(SET CMP0069 NEW) # INTERPROCEDURAL_OPTIMIZATIONS flags
endif()

add_subdirectory(src)

if (LIBSTRATEGY_TESTS)
    include(cmake/doctest.cmake)
    enable_testing()
    add_subdirectory(test)
endif()
