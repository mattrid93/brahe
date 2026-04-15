include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.1
)
FetchContent_MakeAvailable(Catch2)

add_executable(brahe_tests test_main.cpp)
target_link_libraries(brahe_tests PRIVATE brahe_lib Catch2::Catch2WithMain)
target_compile_options(brahe_tests PRIVATE -Wall -Wextra -Wpedantic)

include(Catch)
catch_discover_tests(brahe_tests)
