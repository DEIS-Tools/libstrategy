cmake_minimum_required(VERSION 3.7)

cmake_policy(SET CMP0048 NEW)

include(ExternalProject)

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
set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

option(STRATEGY_GetDependencies "Fetch external dependencies from web." ON)

if (STRATEGY_GetDependencies)
    ExternalProject_add(ptrie
        URL https://github.com/petergjoel/ptrie/archive/v1.1.1.zip
        URL_HASH SHA512=5435D6F8132D273B50ECCD721790D559F7DC67FFDADA75FEC315F4C53705199F31F1DDC80DDA0CB86F883DB27076FA6D324D28516B442FECCE3FFA51213AEF4E
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release
    )

    ExternalProject_add(nlohmann_json
        URL https://github.com/nlohmann/json/archive/v3.7.3.zip
        URL_HASH SHA512=b47a07de9a071cce645a173d084df5dd31f7669154fc00f6c99e0506474d30e8376acaee1d3c79a50def4f25a36042951bfa4fca9a704687e59c368d05053158
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DJSON_BuildTests=OFF -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    )

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

add_subdirectory(${CMAKE_SOURCE_DIR}/src/)

set(LIBSTRATEGY_BuildTests "Build tests of libstrategy" ON)
if(BUILD_TESTING AND LIBSTRATEGY_BuildTests)
    enable_testing()
    add_subdirectory(test)
endif()
