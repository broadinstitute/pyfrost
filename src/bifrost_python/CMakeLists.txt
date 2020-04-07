pybind11_add_module(bifrost_python
    Kmer.h
    UnitigDataDict.h
    UnitigColors.h
    UnitigMapping.h
    UnitigDataProxy.h
    UnitigDataProxy.cpp
    AdjacencyProxy.h
    Pyfrost.h
    NodeView.h
    BifrostDiGraph.h
    bifrost_python.cpp
    )

target_compile_options(bifrost_python PUBLIC "-fvisibility=hidden -march=native")
target_link_options(bifrost_python PUBLIC "-fvisibility=hidden -march=native")
target_link_libraries(bifrost_python PUBLIC bifrost_static)
