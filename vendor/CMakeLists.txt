# Include Bifrost
add_subdirectory(bifrost)
target_compile_options(bifrost_static PUBLIC -fvisibility=hidden)

# Include pybind11
add_subdirectory(pybind11)

# Include testing framework
add_subdirectory(Catch2)

# Include fast hash table
add_subdirectory(robin-hood-hashing)

# Include cereal serialization library
add_subdirectory(cereal)
