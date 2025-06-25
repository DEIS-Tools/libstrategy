add_library(strategy SHARED libz2s.cpp ZonotopStrategy.cpp SimpleTree.cpp)
target_link_libraries(strategy PUBLIC ptrie::ptrie PRIVATE nlohmann_json::nlohmann_json)

add_library(strategyStatic STATIC libz2s.cpp ZonotopStrategy.cpp SimpleTree.cpp)
target_link_libraries(strategyStatic PUBLIC ptrie::ptrie PRIVATE nlohmann_json::nlohmann_json)
set_target_properties(strategyStatic PROPERTIES OUTPUT_NAME strategy)

if(NOT LIBSTRATEGY_OnlyLibrary)
	add_executable(z2s main.cpp ZonotopStrategy.cpp SimpleTree.cpp)
    target_link_libraries(z2s PUBLIC ptrie::ptrie PRIVATE nlohmann_json::nlohmann_json Boost::program_options)
endif()

install(TARGETS strategy
	RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install (FILES  errors.h  SimpleTree.h ZonotopStrategy.h DESTINATION include/libstrategy)


