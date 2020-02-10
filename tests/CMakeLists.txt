add_executable(tests_main
    ${CMAKE_CURRENT_SOURCE_DIR}/test_main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/test_simple_bubble.cpp
    )

target_link_libraries(tests_main PRIVATE bifrost_static Catch2::Catch2)
