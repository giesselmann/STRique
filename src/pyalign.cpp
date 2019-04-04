// \MODULE\-------------------------------------------------------------------------
//
//  CONTENTS      : Top level python bindings
//
//  DESCRIPTION   : Python interface for raw signal alignment
//
//  RESTRICTIONS  : none
//
//  REQUIRES      : none
//
// ---------------------------------------------------------------------------------
// Copyright (c) 2018-2019,  Pay Giesselmann, Max Planck Institute for Molecular Genetics
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// Written by Pay Giesselmann
// ---------------------------------------------------------------------------------

//-- standard headers --------------------------------------------------------------

//-- private headers ---------------------------------------------------------------
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "align_raw.h"

// -- forward declarations ---------------------------------------------------------

// -- exported constants, types, classes -------------------------------------------
namespace py = pybind11;

PYBIND11_MODULE(pyseqan, m) {
    py::class_<align_raw<float, float>> RawAlign(m, "align_raw");
    RawAlign
        .def(py::init<>())
        .def_property("gap_open", &align_raw<float, float>::get_gap_open, &align_raw<float, float>::set_gap_open)
        .def_property("gap_extension", &align_raw<float, float>::get_gap_extension, &align_raw<float, float>::set_gap_extension)
        .def_property("gap_open_h", &align_raw<float, float>::get_gap_open_h, &align_raw<float, float>::set_gap_open_h)
        .def_property("gap_extension_h", &align_raw<float, float>::get_gap_extension_h, &align_raw<float, float>::set_gap_extension_h)
        .def_property("gap_open_v", &align_raw<float, float>::get_gap_open_v, &align_raw<float, float>::set_gap_open_v)
        .def_property("gap_extension_v", &align_raw<float, float>::get_gap_extension_v, &align_raw<float, float>::set_gap_extension_v)
        .def_property("dist_offset", &align_raw<float, float>::get_dist_offset, &align_raw<float, float>::set_dist_offset)
        .def_property("dist_min", &align_raw<float, float>::get_dist_min, &align_raw<float, float>::set_dist_min)
        .def("align_overlap", &align_raw<float, float>::semiglobal,
            py::return_value_policy::take_ownership, "Semi-Global signal alignment",
            py::arg("a"), py::arg("b"))
        ;
}