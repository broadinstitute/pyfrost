# Include Bifrost
add_subdirectory(bifrost)

# Define LEMON graph library
add_library(lemon INTERFACE)
target_include_directories(lemon INTERFACE
    INTERFACE
    $<INSTALL_INTERFACE:include>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/lemon/include>
    )

# Include testing framework
add_subdirectory(Catch2)
