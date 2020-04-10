pybind11_add_module(bifrost_python
    Pyfrost.h
    Kmer.h
    UnitigDataDict.h
    UnitigColors.h
    UnitigColors.cpp
    UnitigMapping.h
    UnitigDataProxy.h
    UnitigDataProxy.cpp
    AdjacencyProxy.h
    AdjacencyProxy.cpp
    NodeView.h
    NodeView.cpp
    NodeBunchIter.h
    NodeBunchIter.cpp
    BifrostDiGraph.h
    BifrostDiGraph.cpp
    bifrost_python.cpp
    )

target_compile_features(bifrost_python PUBLIC cxx_std_14)
target_compile_options(bifrost_python PUBLIC "-fvisibility=hidden -march=native")
target_link_options(bifrost_python PUBLIC "-fvisibility=hidden -march=native")
target_link_libraries(bifrost_python PUBLIC bifrost_static)
