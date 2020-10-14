pybind11_add_module(pyfrostcpp
    UnitigDataDict.h
    Kmer.h
    Kmer.cpp
    KmerCounter.h
    KmerCounter.cpp
    UnitigColors.h
    UnitigColors.cpp
    UnitigMapping.h
    NodeDataDict.h
    NodeDataDict.cpp
    AdjacencyOuterDict.h
    AdjacencyOuterDict.cpp
    AdjacencyInnerDict.h
    AdjacencyInnerDict.cpp
    NodeIterator.h
    NodesDict.h
    NodesDict.cpp
    pyfrost.h
    pyfrost.cpp
    )

target_compile_features(pyfrostcpp PRIVATE cxx_std_14)
target_compile_options(pyfrostcpp PUBLIC -fvisibility=hidden)
target_link_libraries(pyfrostcpp PUBLIC bifrost_dynamic robin_hood cereal -fvisibility=hidden)
