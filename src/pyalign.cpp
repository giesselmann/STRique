// \MODULE\---------------------------------------------------------------
//
//  CONTENTS      : Top level python bindings
//
//  DESCRIPTION   : Python interface for raw signal alignment
//
//  RESTRICTIONS  : none
//
//  REQUIRES      : none
//
// -----------------------------------------------------------------------
//  All rights reserved to Max Planck Institute for Molecular Genetics
//  Berlin, Germany
//  Written by Pay Giesselmann
// -----------------------------------------------------------------------

//-- standard headers ----------------------------------------------------

//-- private headers -----------------------------------------------------
#include "pybind11/pybind11.h"
#include "pybind11/stl.h"
#include "align_raw.h"

// -- forward declarations -----------------------------------------------

// -- exported constants, types, classes ---------------------------------
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