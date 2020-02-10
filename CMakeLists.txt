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
    ExternalProject_add(ptrie-ext
        URL https://github.com/petergjoel/ptrie/archive/v1.1.0.zip
        URL_HASH SHA512=092a8f50ca21d1199b19a10c4e0273c93a717a9f6491998a16bf21d21d37e6537ffd8a06ac41a2b623241da6036546d44b754567441944565e2a16646378cf29
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release
    )

    ExternalProject_add(json-ext
        URL https://github.com/nlohmann/json/archive/v3.7.3.zip
        URL_HASH SHA512=b47a07de9a071cce645a173d084df5dd31f7669154fc00f6c99e0506474d30e8376acaee1d3c79a50def4f25a36042951bfa4fca9a704687e59c368d05053158
	CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION} -DJSON_BuildTests=OFF
    )

    include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
    link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)
endif (STRATEGY_GetDependencies)

# Unset Boost_ROOT here if set, not needed for JSON or source of strategy
if(Boost_ROOT)
	unset(Boost_ROOT)
endif()

add_subdirectory(${CMAKE_SOURCE_DIR}/src/)
