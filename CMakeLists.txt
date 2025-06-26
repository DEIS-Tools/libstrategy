cmake_minimum_required(VERSION 3.28)
project(libstrategy VERSION 1.0.9 LANGUAGES CXX C)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

include(cmake/warnings.cmake)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

option(LIBSTRATEGY_TESTS "libstrategy Unit Tests" ON)
option(LIBSTRATEGY_OnlyLibrary "Build only the library (no binary utilities)" OFF)

include(cmake/ptrie.cmake)
include(cmake/nlohmann_json.cmake)

if(NOT LIBSTRATEGY_OnlyLibrary)
	cmake_policy(SET CMP0069 NEW) # INTERPROCEDURAL_OPTIMIZATIONS flags
	if (CMAKE_MAJOR_VERSION GREATER_EQUAL 3 AND CMAKE_MINOR_VERSION GREATER_EQUAL 30)
    	cmake_policy(SET CMP0167 NEW) # use FindBoost from boost installation (as opposed to CMake own)
	endif ()
	find_package(Boost 1.83 COMPONENTS program_options REQUIRED)
endif()

add_subdirectory(src)

if (LIBSTRATEGY_TESTS)
    enable_testing()
    add_subdirectory(test)
endif()
