# \BUILD\---------------------------------------------------------------
#
#  CONTENTS      : CMake rules
#
#  DESCRIPTION   : Build rules for python SeqAn wrapper
#
#  RESTRICTIONS  : none
#
#  REQUIRES      : none
#
# -----------------------------------------------------------------------
#  All rights reserved to Max Planck Institute for Molecular Genetics
#  Berlin, Germany
#  Written by Pay Giesselmann
# -----------------------------------------------------------------------
cmake_minimum_required(VERSION 3.2.2)
project(repetition CXX)

# jana python extension
file(GLOB repetition_src
    "${CMAKE_CURRENT_SOURCE_DIR}/pyalign.cpp"
)
pybind11_add_module(pyseqan "${repetition_src}")
target_link_libraries(pyseqan)
