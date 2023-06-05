
add_library(figtree SHARED
    figtree.cpp
    KCenterClustering.cpp
    ../include/figtree.h
    KCenterClustering.h)
    
set_target_properties(figtree PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

target_link_libraries(figtree ann_figtree_version)

target_include_directories(figtree PRIVATE ${PROJECT_SOURCE_DIR}/include ${ANN_dir}/include)

add_dependencies(figtree ann_figtree_version)

add_subdirectory(mex)