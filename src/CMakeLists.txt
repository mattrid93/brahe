add_library(brahe_lib brahe.cpp)
target_include_directories(brahe_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_compile_features(brahe_lib PUBLIC cxx_std_17)
target_compile_options(brahe_lib PRIVATE -Wall -Wextra -Wpedantic)

add_executable(brahe main.cpp)
target_link_libraries(brahe PRIVATE brahe_lib)
target_compile_options(brahe PRIVATE -Wall -Wextra -Wpedantic)
