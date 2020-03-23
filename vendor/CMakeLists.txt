# Include Bifrost
add_subdirectory(bifrost)
target_compile_options(bifrost_static PUBLIC "-fvisibility=hidden")
target_link_options(bifrost_static PUBLIC "-fvisibility=hidden")

# Include pybind11
add_subdirectory(pybind11)

# Include testing framework
add_subdirectory(Catch2)
