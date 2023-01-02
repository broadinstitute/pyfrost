cmake_minimum_required(VERSION 3.10..3.18)
project(pyfrost VERSION 0.1
    DESCRIPTION "Python bindings to Bifrost with a NetworkX compatible API"
    LANGUAGES CXX)

# To enable a larger default k-mer size, replace MAX_KMER_SIZE with a larger
# multiple of 32: actual maximum k-mer size will be MAX_KMER_SIZE-1.
SET(MAX_KMER_SIZE "32" CACHE STRING "MAX_KMER_SIZE")
SET(MAX_GMER_SIZE "${MAX_KMER_SIZE}" CACHE STRING "MAX_GMER_SIZE")

# Setup third party dependencies
add_subdirectory(vendor)

# Setup our python library
add_subdirectory(src)

# Setup unit testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/vendor/Catch2/contrib")
    include(CTest)
    include(Catch)
endif()

if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME AND BUILD_TESTING)
    add_subdirectory(tests)
endif()
