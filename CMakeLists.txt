cmake_minimum_required(VERSION 3.12)
project(ugraph VERSION 0.1 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

add_library(ugraph INTERFACE)

target_sources(ugraph
    INTERFACE
    )

target_include_directories(ugraph INTERFACE
    INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    )

# Include Bifrost
add_subdirectory(vendor/bifrost)

# Include testing framework
add_subdirectory(vendor/Catch2)

add_executable(tests_main
    tests/test_main.cpp
    tests/test_linked_dbg.cpp
    )

target_link_libraries(tests_main Catch2::Catch2 ugraph)

