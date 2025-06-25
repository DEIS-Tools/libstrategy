cmake_minimum_required(VERSION 3.28)
project(libstrategy VERSION 1.0.9 LANGUAGES CXX C)

cmake_policy(SET CMP0048 NEW)

include(ExternalProject)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if (NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif (NOT CMAKE_BUILD_TYPE)

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -Wpedantic -I -fPIC")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2 -DNDEBUG -Wall -Wpedantic -fPIC")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

option(LIBSTRATEGY_OnlyLibrary "Build only as library." OFF)
option(STRATEGY_GetDependencies "Fetch external dependencies from web." ON)

if (STRATEGY_GetDependencies)
    ExternalProject_add(ptrie
        URL https://github.com/petergjoel/ptrie/archive/v1.1.1.zip
        URL_HASH SHA512=5435D6F8132D273B50ECCD721790D559F7DC67FFDADA75FEC315F4C53705199F31F1DDC80DDA0CB86F883DB27076FA6D324D28516B442FECCE3FFA51213AEF4E
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release
    )

    include(cmake/nlohmann_json.cmake)

    include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
    link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)
else (STRATEGY_GetDependencies)
	find_package(nlohmann_json 3.7.3 REQUIRED)
	find_package(ptrie 1.1.1 REQUIRED)
endif (STRATEGY_GetDependencies)

# Unset Boost_ROOT here if set, not needed for JSON or source of strategy
if(Boost_ROOT)
	unset(Boost_ROOT)
endif()

if(NOT LIBSTRATEGY_OnlyLibrary)
	cmake_policy(SET CMP0069 NEW)
	cmake_policy(SET CMP0074 NEW)
    cmake_policy(SET CMP0167 NEW) # do not rely on FindBoost/CMake, use FindBoost from boost installation

	find_package(Boost 1.83 COMPONENTS program_options REQUIRED)
	include_directories(${Boost_INCLUDE_DIR})
endif()

add_subdirectory(src)

set(LIBSTRATEGY_BuildTests "Build tests of libstrategy" ON)
if(BUILD_TESTING AND LIBSTRATEGY_BuildTests)
    enable_testing()
    add_subdirectory(test)
endif()
