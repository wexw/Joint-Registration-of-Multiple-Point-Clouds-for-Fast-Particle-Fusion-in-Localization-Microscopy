
# return if matlab not found
if (NOT Matlab_FOUND)
    message("No  suitable MATLAB and/or MX_LIBRARY version found - skipping MATLAB binding!")
    return()
endif()

# GPU
if (CUDA_FOUND)

    add_library(mex_gausstransform SHARED mex_gausstransform.cpp)

    set_target_properties(mex_gausstransform
      PROPERTIES
        PREFIX ""
        SUFFIX .${Matlab_MEX_EXTENSION}
        RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
    )

    add_dependencies(mex_gausstransform gausstransform)

    target_include_directories(mex_gausstransform PRIVATE ${Matlab_INCLUDE_DIRS} ${Matlab_SOURCE_DIR} ../)

    target_link_libraries(mex_gausstransform gausstransform ${Matlab_LIBRARIES})

endif()


# CPU
add_library(mex_gausstransform_cpu SHARED mex_gausstransform_cpu.cpp)

set_target_properties(mex_gausstransform_cpu
  PROPERTIES
    PREFIX ""
    SUFFIX .${Matlab_MEX_EXTENSION}
    RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
)

target_include_directories(mex_gausstransform_cpu PRIVATE ${Matlab_INCLUDE_DIRS} ${Matlab_SOURCE_DIR} ../)

target_link_libraries(mex_gausstransform_cpu ${Matlab_LIBRARIES})

# export mexFunction
if(WIN32)
  SET(CMAKE_SHARED_LINKER_FLAGS "/export:mexFunction")
endif()

# MATLAB package

if(CUDA_FOUND)
    set( binary $<TARGET_FILE:mex_gausstransform> $<TARGET_FILE:mex_gausstransform_cpu> $<TARGET_FILE:gausstransform> )
else()
    set( binary $<TARGET_FILE:mex_gausstransform_cpu> )
endif()

add_custom_target( MATLAB_gausstransform_PACKAGE ALL
  COMMAND ${CMAKE_COMMAND} -E
    make_directory ${matlab_package_build_directory}
  COMMAND ${CMAKE_COMMAND} -E
    copy_if_different ${binary} ${matlab_package_build_directory}
  COMMAND ${CMAKE_COMMAND} -E
    remove ${binary}
  COMMENT "Creating gausstransform Matlab package"
)

if(CUDA_FOUND)
    add_dependencies(MATLAB_gausstransform_PACKAGE mex_gausstransform mex_gausstransform_cpu)
else()
    add_dependencies(MATLAB_gausstransform_PACKAGE mex_gausstransform_cpu)
endif()