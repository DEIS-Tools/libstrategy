
find_package (Boost COMPONENTS unit_test_framework REQUIRED)
include_directories (${TEST_SOURCE_DIR}/src
                     ${Boost_INCLUDE_DIRS}
                     ${libstrategy_SOURCE_DIR}
                     )
add_definitions (-DBOOST_TEST_DYN_LINK)

add_executable (unordered_load unordered_load.cpp)

target_link_libraries(unordered_load ${Boost_LIBRARIES} strategy)

add_test(NAME unordered_load COMMAND unordered_load)
