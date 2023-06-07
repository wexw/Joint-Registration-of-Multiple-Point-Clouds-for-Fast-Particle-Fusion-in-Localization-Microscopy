
set(sources
    src/ANN.cpp
    src/bd_fix_rad_search.cpp
    src/bd_pr_search.cpp
    src/bd_search.cpp
    src/bd_tree.cpp
    src/brute.cpp
    src/kd_dump.cpp
    src/kd_fix_rad_search.cpp
    src/kd_pr_search.cpp
    src/kd_search.cpp
    src/kd_split.cpp
    src/kd_tree.cpp
    src/kd_util.cpp
    src/perf.cpp)
    
set(headers
    include/ANN/ANN.h
    include/ANN/ANNperf.h
    include/ANN/ANNx.h
    src/bd_tree.h
    src/kd_fix_rad_search.h
    src/kd_pr_search.h
    src/kd_search.h
    src/kd_split.h
    src/kd_tree.h
    src/kd_util.h
    src/pr_queue.h
    src/pr_queue_k.h)
    
add_library(ann_figtree_version SHARED ${sources} ${headers})

set_target_properties(ann_figtree_version PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}")

target_include_directories(ann_figtree_version PRIVATE ${ANN_dir}/src ${ANN_dir}/include)