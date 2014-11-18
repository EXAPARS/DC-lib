cmake_minimum_required (VERSION 2.6)

# Compilers
set (CMAKE_CXX_COMPILER icpc)

# Project
project (DC)
find_package (MPI)
set (LIBRARY_OUTPUT_PATH ../lib/${CMAKE_BUILD_TYPE})

# Flags
set (flags "-O2 -mkl:sequential")
if (${HYBRID})
	set (flags "${flags} -DHYBRID")
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
    set (flags "${flags} -DTREE_CREATION")
endif (${TREE})
if (${STATS})
    set (flags "${flags} -DSTATS")
endif (${STATS})
if (${DEBUG})
    set (flags "${flags} -g -O0")
endif (${DEBUG})
if (${VERBOSE})
    set (CMAKE_VERBOSE_MAKEFILE 1)
    set (flags "${flags} -vec-report6")
endif (${VERBOSE})
if (${OMP})
    set (flags "${flags} -openmp")
    add_definitions (-DOMP)
endif (${OMP})
set (CMAKE_CXX_FLAGS "${flags} -std=c++11")

# Headers
include_directories (../include)
include_directories (../src/headers)
include_directories (${MPI_INCLUDE_PATH})
include_directories (${METIS_INCLUDE})

# Sources
file (GLOB_RECURSE source_list ../src/* ../include/*)

# Library
set (lib DC)
if (${HYBRID})
	set (lib ${lib}_Hybrid_${ARCHI})
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
if (${CILK})
    set (lib ${lib}_Cilk)
endif (${CILK})
if (${OMP})
    set (lib ${lib}_OpenMP)
endif (${OMP})
add_library (
	${lib}
    SHARED
	${source_list}
)
