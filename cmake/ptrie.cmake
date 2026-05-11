#set(FETCHCONTENT_QUIET ON)
#set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
include(FetchContent)
FetchContent_Declare(ptrie
        GIT_REPOSITORY https://github.com/DEIS-Tools/ptrie
        GIT_TAG v1.1.2
        GIT_SHALLOW TRUE  # download specific revision only (git clone --depth 1)
        GIT_PROGRESS TRUE # show download progress in Ninja
        USES_TERMINAL_DOWNLOAD TRUE
        EXCLUDE_FROM_ALL # don't build if not used
        FIND_PACKAGE_ARGS 1.1.2)

set(PTRIE_BuildTests OFF CACHE BOOL "Build the unit tests when BUILD_TESTING is enabled.")
set(PTRIE_BuildBenchmark OFF CACHE BOOL "Build the simple benchmark suite")
FetchContent_MakeAvailable(ptrie)

if (ptrie_FOUND) # find_package
   get_target_property(ptrie_INCLUDE_DIRS ptrie::ptrie INTERFACE_INCLUDE_DIRECTORIES)
   message(STATUS "Found ptrie: ${ptrie_INCLUDE_DIRS}")
else (ptrie_FOUND) # fetch_content
   message(STATUS "Got ptrie: ${ptrie_SOURCE_DIR}")
   # Workaround until ptrie exports proper cmake config:
   add_library(ptrie::ptrie INTERFACE IMPORTED GLOBAL)
   target_include_directories(ptrie::ptrie INTERFACE ${ptrie_SOURCE_DIR}/src)
endif (ptrie_FOUND)

if (TARGET ptrie::ptrie)
   message(STATUS "    Available target: ptrie::ptrie")
endif ()
