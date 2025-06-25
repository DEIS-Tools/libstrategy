# set(CMAKE_INCLUDE_CURRENT_DIR ON)

add_library(strategy SHARED ${HEADER_FILES} libz2s.cpp ZonotopStrategy.cpp SimpleTree.cpp)
target_link_libraries(strategy PRIVATE nlohmann_json::nlohmann_json)
if (STRATEGY_GetDependencies)
	add_dependencies(strategy ptrie nlohmann_json)
endif (STRATEGY_GetDependencies)
#target_include_directories (strategy PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

add_library(strategyStatic STATIC ${HEADER_FILES} libz2s.cpp ZonotopStrategy.cpp SimpleTree.cpp)
target_link_libraries(strategyStatic PRIVATE nlohmann_json::nlohmann_json)
if (STRATEGY_GetDependencies)
	add_dependencies(strategyStatic ptrie nlohmann_json)
endif (STRATEGY_GetDependencies)
#target_include_directories (strategyStatic PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_target_properties(strategyStatic PROPERTIES OUTPUT_NAME strategy)

if(NOT LIBSTRATEGY_OnlyLibrary)
	add_executable(z2s ${HEADER_FILES} main.cpp ZonotopStrategy.cpp SimpleTree.cpp)
    target_link_libraries(z2s PRIVATE nlohmann_json::nlohmann_json ${Boost_LIBRARIES})
	if (STRATEGY_GetDependencies)
		add_dependencies(z2s ptrie nlohmann_json)
	endif (STRATEGY_GetDependencies)
endif()

install(TARGETS strategy
	RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)
install (FILES  errors.h  SimpleTree.h ZonotopStrategy.h DESTINATION include/libstrategy)


