# Consider: -Wconversion -Wsign-conversion -Wshadow -Weffc++
set(LIBSTRATEGY_COMMON_WARNINGS -Wpedantic -Wall -Wextra -Wold-style-cast -Wcast-align -Wformat=2 -Wimplicit-fallthrough -Wnon-virtual-dtor -Woverloaded-virtual -Wnull-dereference -Wmisleading-indentation)
if (CMAKE_CXX_COMPILER_ID MATCHES GNU)
    set(LIBSTRATEGY_WARNINGS ${LIBSTRATEGY_COMMON_WARNINGS} -Wuseless-cast -Wduplicated-cond -Wduplicated-branches -Wlogical-op)
elseif (CMAKE_CXX_COMPILER_ID MATCHES Clang)
    set(LIBSTRATEGY_WARNINGS ${LIBSTRATEGY_COMMON_WARNINGS})
elseif (CMAKE_CXX_COMPILER_ID MATCHES AppleClang)
    set(LIBSTRATEGY_WARNINGS ${LIBSTRATEGY_COMMON_WARNINGS})
elseif (CMAKE_CXX_COMPILER_ID MATCHES MSVC)
    set(LIBSTRATEGY_WARNINGS /W4 /w14619 /w14265 /w14545 /w14546 /w14640 /w14905 /w14906 /w14928)
endif ()
# Use target_compile_options to avoid warnings on third party code
# add_compile_options(${LIBSTRATEGY_WARNINGS})
