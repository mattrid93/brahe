include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG v3.7.1
)
FetchContent_MakeAvailable(Catch2)

add_executable(brahe_tests
    test_types.cpp
    test_body_builder.cpp
    test_ephemeris.cpp
    test_frame_transforms.cpp
    test_conic_invariants.cpp
    test_conic_classification.cpp
    test_conic_roundtrip.cpp
    test_anomaly_conversions.cpp
    test_two_body_propagation.cpp
    test_invalid_input.cpp
    test_event_api.cpp
    test_event_root_functions.cpp
    test_event_refinement.cpp
    test_event_coarse_scan.cpp
    test_event_detector_impact.cpp
    test_event_detector_soi_exit.cpp
    test_event_detector_child_entry.cpp
    test_event_ordering.cpp
    test_event_boundaries.cpp
    test_event_invalid_input.cpp
    test_event_properties.cpp
    test_event_fuzz.cpp
    test_event_moving_bodies.cpp
)
target_link_libraries(brahe_tests PRIVATE brahe_lib Catch2::Catch2WithMain)
target_compile_options(brahe_tests PRIVATE -Wall -Wextra -Wpedantic)

include(Catch)
catch_discover_tests(brahe_tests)
