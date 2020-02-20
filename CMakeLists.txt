cmake_minimum_required(VERSION 3.12)
project(ugraph VERSION 0.1
    DESCRIPTION "Pan-genome graph assisted variant caller for microbial genomes"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Setup third party dependencies
add_subdirectory(vendor)

# Define ugraph (header only) library
add_library(ugraphlib INTERFACE)

target_sources(ugraphlib INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/include/ugraph/bifrost_lemon_graph.h)

target_include_directories(ugraphlib INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)

target_link_libraries(ugraphlib INTERFACE lemon bifrost_static)

# ugraph CLI targets
add_subdirectory(cli)

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
