cmake_minimum_required (VERSION 2.6)

# Compilers
set (LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})
set (CMAKE_CXX_COMPILER icpc)

# Project
project (DC-lib)
if (${DISTRI} STREQUAL "XMPI")
    find_package (MPI)
endif (${DISTRI} STREQUAL "XMPI")

# Flags
set (flags "-O2 -mkl:sequential")
add_definitions (-D${DISTRI} -D${SHARED})
if (${SHARED} STREQUAL "OMP")
    set (flags "${flags} -openmp")
endif (${SHARED} STREQUAL "OMP")
if (${HYBRID})
    add_definitions (-DHYBRID)
endif (${HYBRID})
if (${ARCHI} STREQUAL "SSE")
    set (flags "${flags} -msse")
    add_definitions (-DVEC_SIZE=2)
elseif (${ARCHI} STREQUAL "AVX")
    set (flags "${flags} -mavx")
    add_definitions (-DVEC_SIZE=4)
elseif (${ARCHI} STREQUAL "MIC")
    set (flags "${flags} -mmic")
    add_definitions (-DVEC_SIZE=8)
endif (${ARCHI} STREQUAL "SSE")
if (${TREE})
    add_definitions (-DTREE_CREATION)
endif (${TREE})
if (${STATS})
    add_definitions (-DSTATS)
endif (${STATS})
if (${DEBUG})
    set (flags "${flags} -g -O0")
endif (${DEBUG})
if (${VERBOSE})
    set (CMAKE_VERBOSE_MAKEFILE 1)
    set (flags "${flags} -vec-report6")
endif (${VERBOSE})
set (CMAKE_CXX_FLAGS "${flags} -std=c++11")

# Headers
include_directories (../include)
include_directories (../src/headers)
if (${DISTRI} STREQUAL "XMPI")
    include_directories (${MPI_INCLUDE_PATH})
endif (${DISTRI} STREQUAL "XMPI")
if (${TREE})
    include_directories (${METIS_INCLUDE})
endif (${TREE})

# Sources
file (GLOB_RECURSE source_list ../src/* ../include/*)

# Library
set (lib DC_${DISTRI}_${SHARED})
if (${HYBRID})
	set (lib ${lib}_HYBRID_${ARCHI})
elseif (${ARCHI} STREQUAL "MIC")
	set (lib ${lib}_${ARCHI})
endif (${HYBRID})
if (${TREE})
	set (lib ${lib}_TreeCreation)
endif (${TREE})
if (${STATS})
	set (lib ${lib}_Stats)
endif (${STATS})
if (${DEBUG})
	set (lib ${lib}_Debug)
endif (${DEBUG})

add_library (
	${lib}
    SHARED
	${source_list}
)
