add_library(SimpleTree OBJECT SimpleTree.cpp)
target_link_libraries(SimpleTree PUBLIC ptrie::ptrie PRIVATE nlohmann_json::nlohmann_json)

add_library(strategy SHARED libz2s.cpp ZonotopStrategy.cpp)
target_link_libraries(strategy PUBLIC SimpleTree)

add_library(strategyStatic STATIC libz2s.cpp ZonotopStrategy.cpp)
target_link_libraries(strategyStatic PUBLIC SimpleTree)
set_target_properties(strategyStatic PROPERTIES OUTPUT_NAME strategy)

if(NOT LIBSTRATEGY_OnlyLibrary)
	add_executable(z2s main.cpp ZonotopStrategy.cpp)
    target_link_libraries(z2s PUBLIC SimpleTree Boost::program_options)
endif()

install(TARGETS strategy strategyStatic)
install(FILES  errors.h  SimpleTree.h ZonotopStrategy.h DESTINATION include/libstrategy)
