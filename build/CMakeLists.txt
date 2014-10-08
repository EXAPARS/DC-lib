cmake_minimum_required (VERSION 2.6)

# Compilers
set (CMAKE_CXX_COMPILER icpc)

# Project
project (DC_lib)
find_package (MPI)
set (LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

# Flags
set (flags "-O2 -mkl:sequential")
if (${hybrid})
	set (flags "${flags} -DHYBRID")
endif (${hybrid})
if (${vec} STREQUAL "sse")
    set (flags "${flags} -msse")
    add_definitions (-DVEC_SIZE=2)
elseif (${vec} STREQUAL "avx")
    set (flags "${flags} -mavx")
    add_definitions (-DVEC_SIZE=4)
elseif (${vec} STREQUAL "mic")
    set (flags "${flags} -mmic")
    add_definitions (-DVEC_SIZE=8)
endif (${vec} STREQUAL "sse")
if (${dump})
    set (flags "${flags} -DCREATE_PERM_AND_TREE")
endif (${dump})
if (${stats})
    set (flags "${flags} -DSTATS")
endif (${stats})
if (${debug})
    set (flags "${flags} -g -O0")
endif (${debug})
if (${verbose})
    set (CMAKE_VERBOSE_MAKEFILE 1)
    set (flags "${flags} -vec-report6")
endif (${verbose})
set (CMAKE_CXX_FLAGS "${flags} -std=c++11")

# Headers
include_directories (../include)
include_directories (../src/headers)
include_directories (${MPI_INCLUDE_PATH})
include_directories (${METIS_INCLUDE_PATH})

# Sources
file (GLOB_RECURSE source_list ../src/* ../include/*)

# Library
set (lib DC_lib)
if (${hybrid})
	set (lib ${lib}_HYBRID)
endif (${hybrid})
if (${vec} STREQUAL "sse")
    set (lib ${lib}_SSE)
elseif (${vec} STREQUAL "avx")
    set (lib ${lib}_AVX)
elseif (${vec} STREQUAL "mic")
    set (lib ${lib}_MIC)
endif (${vec} STREQUAL "sse")
if (${dump})
	set (lib ${lib}_Dump)
endif (${dump})
if (${stats})
	set (lib ${lib}_Stats)
endif (${stats})
if (${debug})
	set (lib ${lib}_Debug)
endif (${debug})
add_library (
	${lib}
    SHARED
	${source_list}
)

# Link
#if (${dump})
#    target_link_libraries (
#        ${lib}
#        ${MPI_LIBRARIES}
#        ${METIS_LIBRARIES}
#    )
#else (${dump})
#    target_link_libraries (
#        ${lib}
#        ${MPI_LIBRARIES}
#    )
#endif (${dump})
