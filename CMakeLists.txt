cmake_minimum_required(VERSION 3.12)

cmake_policy(SET CMP0048 NEW)
cmake_policy(SET CMP0069 NEW)
cmake_policy(SET CMP0074 NEW)

set(CMAKE_CXX_STANDARD 17)
if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
project(libstrategy VERSION 1.0.0 LANGUAGES CXX C)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -Wpedantic -I -fPIC")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG -Wall -Wpedantic -fPIC")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/ptrie)
# Unset Boost_ROOT here if set, not needed for JSON or source of strategy
unset(Boost_ROOT)
add_subdirectory(${CMAKE_SOURCE_DIR}/third_party/json)
add_subdirectory(${CMAKE_SOURCE_DIR}/src/)
