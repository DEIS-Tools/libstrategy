include_directories(${CMAKE_SOURCE_DIR}/src)

add_executable (utilities_test utilities_test.cpp)
add_test(NAME utilities_test COMMAND utilities_test)

add_executable (SimpleTree_test SimpleTree_test.cpp)
target_link_libraries(SimpleTree_test PRIVATE SimpleTree doctest::doctest_with_main)
add_test(NAME SimpleTree_test COMMAND SimpleTree_test)

add_executable (unordered_load unordered_load.cpp)
target_link_libraries(unordered_load PRIVATE strategy doctest::doctest_with_main)
add_test(NAME unordered_load COMMAND unordered_load)

add_executable (inconsistent_lookup inconsistent_lookup.cpp)
target_link_libraries(inconsistent_lookup PRIVATE strategy doctest::doctest_with_main)
add_test(NAME inconsistent_lookup COMMAND inconsistent_lookup)
set_tests_properties(inconsistent_lookup PROPERTIES ENVIRONMENT STRATEGY_DIR=${CMAKE_CURRENT_SOURCE_DIR}/strategies)
