cmake_minimum_required(VERSION 3.12)
project(ugraph VERSION 0.1
    DESCRIPTION "Pan-genome assisted variant caller for microbial genomes"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 14)

# Include Bifrost
add_subdirectory(vendor/bifrost)

# Define LEMON graph library
add_library(lemon INTERFACE)
target_include_directories(lemon INTERFACE
    INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/vendor/lemon/include>
    )

# Include testing framework
add_subdirectory(vendor/Catch2)

# Our directory with CLI programs
add_subdirectory(cli)

# Our directory with tests
add_subdirectory(tests)

