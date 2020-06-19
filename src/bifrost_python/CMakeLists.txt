pybind11_add_module(bifrost_python
    Pyfrost.h
    UnitigDataDict.h
    Kmer.h
    Kmer.cpp
    UnitigColors.h
    UnitigColors.cpp
    UnitigMapping.h
    UnitigDataProxy.h
    UnitigDataProxy.cpp
    AdjacencyProxy.h
    AdjacencyProxy.cpp
    NodeIterator.h
    NodeView.h
    NodeView.cpp
    EdgeView.h
    EdgeView.cpp
    DegreeView.h
    BifrostDiGraph.h
    BifrostDiGraph.cpp
    bifrost_python.cpp
    )

target_compile_features(bifrost_python PUBLIC cxx_std_14)
target_compile_options(bifrost_python PUBLIC "-fvisibility=hidden -march=native")
target_link_options(bifrost_python PUBLIC "-fvisibility=hidden -march=native")
target_link_libraries(bifrost_python PUBLIC bifrost_static)
