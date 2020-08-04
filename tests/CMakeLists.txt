add_executable(tests_main
    test_main.cpp
    test_mccortex.cpp
    test_node_removal.cpp
    )

target_link_libraries(tests_main PRIVATE bifrost_static Catch2::Catch2)

# Copy data to build directory
file(COPY data/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)

# Automatically register tests with CTest
catch_discover_tests(tests_main)
