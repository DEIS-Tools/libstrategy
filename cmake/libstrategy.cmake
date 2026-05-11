include(FetchContent)
FetchContent_Declare(libstrategy
        GIT_REPOSITORY https://github.com/DEIS-Tools/libstrategy
        GIT_TAG fixing_version_check
        GIT_SHALLOW TRUE  # download specific revision only (git clone --depth 1)
        GIT_PROGRESS TRUE # show download progress in Ninja
        USES_TERMINAL_DOWNLOAD TRUE
        EXCLUDE_FROM_ALL # don't build if not used
        FIND_PACKAGE_ARGS 1.1.2)

set(LIBSTRATEGY_TESTS OFF CACHE BOOL "libstrategy Unit Tests")
set(LIBSTRATEGY_OnlyLibrary ON CACHE BOOL "Build only the library (no binary utilities)")
FetchContent_MakeAvailable(libstrategy)

if (libstrategy_FOUND) # find_package
    message(STATUS "Found libstrategy: ${libstrategy_DIR}")
else (libstrategy_FOUND) # fetch_content
    message(STATUS "Fetched libstrategy: ${libstrategy_SOURCE_DIR}")
endif (libstrategy_FOUND)

if (TARGET libstrategy::strategy)
    message(STATUS "    Available target: libstrategy::strategy")
endif ()
if (TARGET libstrategy::strategyStatic)
    message(STATUS "    Available target: libstrategy::strategyStatic")
endif ()
