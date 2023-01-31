pybind11_add_module(pyfrostcpp
        UnitigDataDict.h
        Kmer.h
        Kmer.cpp
        Minimizers.h
        Minimizers.cpp
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
        JunctionTree.h
        JunctionTree.cpp
        LinkDB.h
        LinkDB.cpp
        MemLinkDB.h
        LinkAnnotator.h
        LinkAnnotator.cpp
        Serialize.h
        Neighbors.h
        pyfrost.h
        pyfrost.cpp
        )

target_compile_features(pyfrostcpp PRIVATE cxx_std_14)
target_compile_options(pyfrostcpp PUBLIC -fvisibility=hidden -fPIC)
target_link_libraries(pyfrostcpp PUBLIC bifrost_static robin_hood cereal optional -fvisibility=hidden)
