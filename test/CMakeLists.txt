
find_package (Boost COMPONENTS unit_test_framework REQUIRED)
include_directories (${TEST_SOURCE_DIR}/src
                     ${Boost_INCLUDE_DIRS}
                     ${libstrategy_SOURCE_DIR}
                     )
add_definitions (-DBOOST_TEST_DYN_LINK)

add_executable (unordered_load unordered_load.cpp)
add_executable (inconsistent_lookup inconsistent_lookup.cpp)

target_link_libraries(unordered_load ${Boost_LIBRARIES} strategy)
target_link_libraries(inconsistent_lookup ${Boost_LIBRARIES} strategy)

add_test(NAME unordered_load COMMAND unordered_load)
add_test(NAME inconsistent_lookup COMMAND inconsistent_lookup)
set_tests_properties(inconsistent_lookup PROPERTIES
    ENVIRONMENT STRATEGY_DIR=${CMAKE_CURRENT_SOURCE_DIR}/strategies)
