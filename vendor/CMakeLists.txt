# Include Bifrost
add_subdirectory(bifrost)
target_compile_options(bifrost_static PUBLIC -fvisibility=hidden)

# Include pybind11
add_subdirectory(pybind11)

# Include testing framework
add_subdirectory(Catch2)

# Include fast hash table
add_library(robin_hood INTERFACE)
target_include_directories(robin_hood INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/robin-hood-hashing/include)

# Include cereal serialization library
add_library(cereal INTERFACE)
target_include_directories(cereal INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/cereal/include)

add_library(optional INTERFACE)
target_include_directories(optional INTERFACE ${CMAKE_CURRENT_SOURCE_DIR}/optional/include)
