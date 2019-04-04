# \BUILD\---------------------------------------------------------------
#
#  CONTENTS      : CMake rules
#
#  DESCRIPTION   : Top-Level build rules for repetition
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# -----------------------------------------------------------------------
# Copyright (c) 2018-2019, Pay Giesselmann, Max Planck Institute for Molecular Genetics
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
# Written by Pay Giesselmann
# -----------------------------------------------------------------------
cmake_minimum_required(VERSION 3.2.2)
project(pyseqan VERSION 0.2.0 LANGUAGES CXX)

# global settings
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set(BUILD_SHARED_LIBS OFF CACHE BOOL "Build Shared Lib" FORCE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON CACHE BOOL "Build position independent code" FORCE)
include(ExternalProject)

# check openmp
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
             
# set up SeqAn environment
set(CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/submodules/seqan/util/cmake")
set(SEQAN_INCLUDE_PATH "${CMAKE_SOURCE_DIR}/submodules/seqan/include")
set(SEQAN_HAS_OPENMP ${OPENMP_FOUND} CACHE BOOL "SeqAn has openmp")
set(SEQAN_ARCH_SSE4 ON CACHE BOOL "Enable SSE instructions")
set(SEQAN_ARCH_AVX2 ON CACHE BOOL "Enable AVX instructions")
set(SEQAN_ENABLE_DEBUG OFF CACHE BOOL "SeqAn enable debug" FORCE)
find_package(SeqAn REQUIRED)

# set up build environment
add_definitions (${SEQAN_DEFINITIONS})
include_directories (${SEQAN_INCLUDE_DIRS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")

# set platform dependencies
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    message(STATUS "Set flags for GNUCXX")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -pedantic -O3 -DNDEBUG -Wfatal-errors")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lpthread -lrt")
endif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    message(STATUS "Set flags for MSVC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /std:c++14 /W2 /wd4996 /Ox -D_CRT_SECURE_NO_WARNINGS")
endif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")

# add sources
add_subdirectory("${CMAKE_SOURCE_DIR}/submodules/pybind11")
add_subdirectory("${CMAKE_SOURCE_DIR}/src")

# config summary
message(STATUS "SeqAn has openmp: " ${SEQAN_HAS_OPENMP})
message(STATUS "ZLIB library: " ${ZLIB_STATIC_LIBRARY})
message(STATUS "Position independent code: " ${POSITION_INDEPENDENT_CODE})

