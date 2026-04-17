find_package(pybind11 CONFIG REQUIRED)

pybind11_add_module(_brahe brahe_py.cpp)
target_link_libraries(_brahe PRIVATE brahe_lib)
target_compile_features(_brahe PRIVATE cxx_std_17)
target_compile_options(_brahe PRIVATE -Wall -Wextra -Wpedantic)

install(TARGETS _brahe LIBRARY DESTINATION brahe)
