add_executable(tests_main
    ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_mccortex.cpp
    )

target_link_libraries(tests_main PRIVATE bifrost_static Catch2::Catch2)

# Copy data to build directory
file(COPY data/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)

# Automatically register tests with CTest
catch_discover_tests(tests_main)
