include(GNUInstallDirs)

set(CMAKE_VISIBILITY_INLINES_HIDDEN ON)

add_library(SimpleTree OBJECT SimpleTree.cpp)
target_link_libraries(SimpleTree PUBLIC ptrie::ptrie PRIVATE nlohmann_json::nlohmann_json)

add_library(ZonotopStrategy OBJECT ZonotopStrategy.cpp)
target_link_libraries(ZonotopStrategy PUBLIC ptrie::ptrie PRIVATE nlohmann_json::nlohmann_json)

add_library(libz2s OBJECT libz2s.cpp)
target_link_libraries(libz2s PUBLIC ptrie::ptrie)

add_library(strategy SHARED $<TARGET_OBJECTS:libz2s> $<TARGET_OBJECTS:ZonotopStrategy> $<TARGET_OBJECTS:SimpleTree>)
target_link_libraries(strategy PUBLIC ptrie::ptrie)
target_include_directories(strategy PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
add_library(libstrategy::strategy ALIAS strategy)

add_library(strategyStatic STATIC $<TARGET_OBJECTS:libz2s> $<TARGET_OBJECTS:ZonotopStrategy> $<TARGET_OBJECTS:SimpleTree>)
target_link_libraries(strategyStatic PUBLIC ptrie::ptrie)
target_include_directories(strategyStatic PUBLIC $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(strategyStatic PROPERTIES OUTPUT_NAME strategy-static)
add_library(libstrategy::strategyStatic ALIAS strategyStatic)

if(NOT LIBSTRATEGY_OnlyLibrary)
    add_executable(z2s main.cpp)
    target_link_libraries(z2s PUBLIC ZonotopStrategy SimpleTree Boost::program_options)
endif()

include(CMakePackageConfigHelpers)

set(CONFIG_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/libstrategy)

configure_package_config_file(libstrategyConfig.cmake.in
        ${CMAKE_CURRENT_BINARY_DIR}/libstrategyConfig.cmake
        INSTALL_DESTINATION ${CONFIG_INSTALL_DIR})
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/libstrategyConfigVersion.cmake
        VERSION        ${PROJECT_VERSION}
        COMPATIBILITY  SameMajorVersion)

install(TARGETS strategy strategyStatic
        EXPORT libstrategyTargets)
install(FILES errors.h libz2s.h SimpleTree.h ZonotopStrategy.h
        DESTINATION include/libstrategy)
install(EXPORT libstrategyTargets
        FILE libstrategyTargets.cmake
        NAMESPACE libstrategy::
        DESTINATION ${CONFIG_INSTALL_DIR})
install(FILES
        ${CMAKE_CURRENT_BINARY_DIR}/libstrategyConfig.cmake
        ${CMAKE_CURRENT_BINARY_DIR}/libstrategyConfigVersion.cmake
        DESTINATION ${CONFIG_INSTALL_DIR})
