
add_executable(run_tests test.cpp)

set_target_properties(run_tests PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

target_include_directories(run_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(run_tests figtree)

add_dependencies(run_tests figtree)