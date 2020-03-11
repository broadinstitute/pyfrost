#include <string>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "Kmer.h"

namespace py = pybind11;


PYBIND11_MODULE(bifrost_python, m) {
    m.doc() = R"doc(
        Python bindings for Bifrost
        ===========================

        This module provides a Python interface to the Bifrost colored compacted de Bruijn graph library, with a
        NetworkX compatible API.
    )doc";

    pyfrost::define_Kmer(m);
}
