#set(FETCHCONTENT_QUIET ON)
#set(FETCHCONTENT_UPDATES_DISCONNECTED ON)
include(FetchContent)
FetchContent_Declare(Boost
        # URL https://github.com/boostorg/boost/releases/download/boost-1.91.0-1/boost-1.91.0-1-cmake.7z
        URL https://people.cs.aau.dk/~marius/mirrors/boost/boost-1.91.0-1-cmake.7z
        URL_HASH SHA256=29c7d4f4ac36ad853b6765d03571ea60d90286775df026b4efd9f3281131972b
        #GIT_REPOSITORY https://github.com/boostorg/boost.git
        #GIT_TAG boost-1.91.0-1
        #GIT_SHALLOW TRUE  # download specific revision only (git clone --depth 1)
        #GIT_PROGRESS TRUE # show download progress in Ninja
        USES_TERMINAL_DOWNLOAD TRUE
        EXCLUDE_FROM_ALL # don't build if not used
        FIND_PACKAGE_ARGS 1.91.0 COMPONENTS program_options unit_test_framework)

set(BOOST_ENABLE_MPI OFF CACHE BOOL "Boost.MPI and its dependents (requires MPI, CMake 3.10)")
set(BOOST_ENABLE_PYTHON OFF CACHE BOOL "Boost.Python and its dependents (requires Python, CMake 3.14)")
set(BUILD_TESTING OFF CACHE BOOL "Build the tests.")
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build shared libraries")


FetchContent_MakeAvailable(Boost)

if (Boost_FOUND) # find_package
   message(STATUS "Found Boost: ${Boost_INCLUDE_DIRS}")
else (Boost_program_options_FOUND) # fetch_content
   message(STATUS "Got Boost: ${Boost_SOURCE_DIR}")
endif (Boost_FOUND)

if (TARGET Boost::program_options)
   message(STATUS "    Available target: Boost::program_options")
endif ()
if (TARGET Boost::unit_test_framework)
   message(STATUS "    Available target: Boost::unit_test_framework")
endif ()

