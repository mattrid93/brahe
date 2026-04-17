add_library(brahe_lib body_system.cpp two_body.cpp event_detector.cpp patcher.cpp
                      trajectory_builder.cpp)
target_include_directories(brahe_lib PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_compile_features(brahe_lib PUBLIC cxx_std_17)
target_compile_options(brahe_lib PRIVATE -Wall -Wextra -Wpedantic)
set_target_properties(brahe_lib PROPERTIES POSITION_INDEPENDENT_CODE ON)

add_executable(brahe main.cpp)
target_link_libraries(brahe PRIVATE brahe_lib)
target_compile_options(brahe PRIVATE -Wall -Wextra -Wpedantic)
