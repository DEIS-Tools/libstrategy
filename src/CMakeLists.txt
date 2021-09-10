cmake_minimum_required(VERSION 3.7)
project(libstrategy C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_INCLUDE_CURRENT_DIR ON)
option(LIBSTRATEGY_OnlyLibrary "Build only as library." OFF)

if(NOT LIBSTRATEGY_OnlyLibrary)
	cmake_policy(SET CMP0069 NEW)
	cmake_policy(SET CMP0074 NEW)

	find_package(Boost 1.66 COMPONENTS program_options REQUIRED)
	include_directories(${Boost_INCLUDE_DIR})
endif()

add_library(strategy SHARED ${HEADER_FILES} libz2s.cpp ZonotopStrategy.cpp SimpleTree.cpp)
if (STRATEGY_GetDependencies)
	add_dependencies(strategy ptrie nlohmann_json)
endif (STRATEGY_GetDependencies)
target_include_directories (strategy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(strategyStatic STATIC ${HEADER_FILES} libz2s.cpp ZonotopStrategy.cpp SimpleTree.cpp)
if (STRATEGY_GetDependencies)
	add_dependencies(strategyStatic ptrie nlohmann_json)
endif (STRATEGY_GetDependencies)
target_include_directories (strategyStatic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(strategyStatic PROPERTIES OUTPUT_NAME strategy)

if(NOT LIBSTRATEGY_OnlyLibrary)
	add_executable(z2s ${HEADER_FILES} main.cpp ZonotopStrategy.cpp SimpleTree.cpp)
	if (STRATEGY_GetDependencies)
		add_dependencies(z2s ptrie nlohmann_json)
	endif (STRATEGY_GetDependencies)
	target_link_libraries(z2s PRIVATE stdc++fs ${Boost_LIBRARIES})
endif()


install(TARGETS strategy
	RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install (FILES  errors.h  SimpleTree.h ZonotopStrategy.h DESTINATION include/libstrategy)


