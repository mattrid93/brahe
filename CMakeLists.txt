cmake_minimum_required(VERSION 3.24)

project(brahe
    VERSION 0.1.0
    DESCRIPTION "Brahe"
    LANGUAGES CXX
)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Prevent in-source builds
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
    message(FATAL_ERROR "In-source builds are not allowed. Use a build directory.")
endif()

add_subdirectory(src)

option(BRAHE_BUILD_TESTS "Build tests" ON)
if(BRAHE_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()
