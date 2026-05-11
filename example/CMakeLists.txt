cmake_minimum_required(VERSION 3.28)
project(example VERSION 1.0.9)

include(nlohmann_json)
include(ptrie)
include(libstrategy)

add_executable(example example.cpp)
target_link_libraries(example PRIVATE libstrategy::strategy)
add_test(NAME example COMMAND example)
