cmake_minimum_required(VERSION 3.7)
project(z2s C CXX)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_INCLUDE_CURRENT_DIR ON)

find_package(nlohmann_json 3.2.0 REQUIRED)
find_package(Boost 1.66 COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR} ../include)

add_library(zonotop2strategy SHARED ${HEADER_FILES} libz2s.cpp ZonotopStrategy.cpp SimpleTree.cpp)
target_include_directories (zonotop2strategy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
add_executable(z2s ${HEADER_FILES} main.cpp ZonotopStrategy.cpp SimpleTree.cpp)
target_link_libraries(z2s PRIVATE stdc++fs nlohmann_json::nlohmann_json ${Boost_LIBRARIES} /home/petko/Devel/ptrie/lib/libptrie.a)
target_link_libraries(zonotop2strategy PRIVATE nlohmann_json::nlohmann_json  /home/petko/Devel/ptrie/lib/libptrie.a)




