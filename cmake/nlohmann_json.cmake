find_package(nlohmann_json 3.12.0 QUIET)
if (nlohmann_json_FOUND)
   message(STATUS "Found nlohmann_json: ${nlohmann_json_SOURCE_DIR}")
else (nlohmann_json_FOUND)
   message(STATUS "Failed to find nlohmann_json, going to fetch from source")
   set(JSON_BuildTests                 OFF CACHE BOOL "Build the unit tests when BUILD_TESTING is enabled.")
   set(JSON_CI                         OFF CACHE BOOL "Enable CI build targets.")
   set(JSON_Diagnostics                OFF CACHE BOOL "Use extended diagnostic messages.")
   set(JSON_Diagnostic_Positions       OFF CACHE BOOL "Enable diagnostic positions.")
   set(JSON_GlobalUDLs                 ON  CACHE BOOL "Place user-defined string literals in the global namespace.")
   set(JSON_ImplicitConversions        ON  CACHE BOOL "Enable implicit conversions.")
   set(JSON_DisableEnumSerialization   OFF CACHE BOOL "Disable default integer enum serialization.")
   set(JSON_LegacyDiscardedValueComparison  OFF CACHE BOOL "Enable legacy discarded value comparison.")
   set(JSON_Install                    OFF CACHE BOOL "Install CMake targets during install step.")
   set(JSON_MultipleHeaders            ON  CACHE BOOL "Use non-amalgamated version of the library.")
   set(JSON_SystemInclude              OFF CACHE BOOL "Include as system headers (skip for clang-tidy).")
   set(FETCHCONTENT_QUIET ON)
   set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
   include(FetchContent)
   FetchContent_Declare(nlohmann_json
       GIT_REPOSITORY https://github.com/nlohmann/json
       GIT_TAG v3.12.0
       GIT_SHALLOW TRUE  # download specific revision only (git clone --depth 1)
       GIT_PROGRESS TRUE # show download progress in Ninja
       FIND_PACKAGE_ARGS NAMES nlohmann_json
       USES_TERMINAL_DOWNLOAD TRUE)
   FetchContent_MakeAvailable(nlohmann_json)
   message(STATUS "Got nlohmann_json: ${nlohmann_json_SOURCE_DIR}")
endif (nlohmann_json_FOUND)
