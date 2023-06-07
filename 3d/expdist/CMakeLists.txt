
if(CUDA_FOUND)

    set(sources
        expdist.cu
        expdist_functions.cuh
        expdist_ref.cu
        matrix_functions.cuh
        expdist.h
        expdist_ref.h
        ${CUB_ROOT_DIR}/cub/cub.cuh)

    cuda_add_library(expdist SHARED ${sources})

    target_include_directories(expdist PRIVATE ${CUB_ROOT_DIR})

    # set output directory
    set_target_properties(expdist
      PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    )
    
endif()

add_subdirectory(mex)
