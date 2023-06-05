
add_executable(sample sample.cpp)

set_target_properties(sample PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

target_include_directories(sample PRIVATE ${PROJECT_SOURCE_DIR}/include)

target_link_libraries(sample figtree)

add_dependencies(sample figtree)