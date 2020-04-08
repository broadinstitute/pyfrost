cmake_minimum_required(VERSION 3.12..3.16)
project(pyfrost VERSION 0.1
    DESCRIPTION "Python bindings to Bifrost with a NetworkX compatible API"
    LANGUAGES CXX)

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
