cmake_minimum_required(VERSION 3.7)
project(z2s C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(Boost 1.66 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

add_library(strategy SHARED ${HEADER_FILES} libz2s.cpp ZonotopStrategy.cpp SimpleTree.cpp)
target_include_directories (strategy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(strategy PRIVATE nlohmann_json::nlohmann_json ptrie)

add_library(strategyStatic STATIC ${HEADER_FILES} libz2s.cpp ZonotopStrategy.cpp SimpleTree.cpp)
target_include_directories (strategyStatic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(strategyStatic PRIVATE nlohmann_json::nlohmann_json ptrie)
set_target_properties(strategyStatic PROPERTIES OUTPUT_NAME strategy)

add_executable(z2s ${HEADER_FILES} main.cpp ZonotopStrategy.cpp SimpleTree.cpp)
target_link_libraries(z2s PRIVATE stdc++fs nlohmann_json::nlohmann_json ptrie ${Boost_LIBRARIES})


install(TARGETS strategy
	RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install (FILES  errors.h  libz2s.h SimpleTree.h ZonotopStrategy.h DESTINATION include/libstrategy)


