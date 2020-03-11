cmake_minimum_required(VERSION 3.12)
project(pyfrost VERSION 0.1
    DESCRIPTION "Pan-genome graph assisted variant caller for microbial genomes"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Setup third party dependencies
add_subdirectory(vendor)

pybind11_add_module(bifrost_python
    src/Kmer.h
    src/KmerOnUnitig.h
    src/NodeView.h
    src/BifrostGraph.h
    src/ColoredBifrostGraph.h
    src/bifrost_python.cpp
    )

target_link_libraries(bifrost_python PUBLIC bifrost_static)

# Setup unit testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/vendor/Catch2/contrib")
    include(CTest)
    include(Catch)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    message("Enabling testing")
    add_subdirectory(tests)
endif()
