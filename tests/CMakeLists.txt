add_executable(tests_main
    ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_simple_bubble.cpp
    )

target_link_libraries(tests_main PRIVATE ugraphlib Catch2::Catch2)

# Copy data to build directory
message(${CMAKE_CURRENT_BINARY_DIR})
file(COPY data/ DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/data)

# Automatically register tests with CTest
catch_discover_tests(tests_main)
