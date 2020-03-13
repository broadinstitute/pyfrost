pybind11_add_module(bifrost_python
    Kmer.h
    Unitig.h
    UnitigDataDict.h
    KmerOnUnitig.h
    Pyfrost.h
    NodeView.h
    NodeView.cpp
    BifrostDiGraph.h
    ColoredBifrostGraph.h
    bifrost_python.cpp
    )

target_link_libraries(bifrost_python PUBLIC bifrost_static)
