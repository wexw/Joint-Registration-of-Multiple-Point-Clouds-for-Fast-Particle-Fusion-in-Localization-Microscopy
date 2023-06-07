
find_package( Matlab COMPONENTS MX_LIBRARY )

if( NOT Matlab_FOUND )
  message( STATUS "No suitable Matlab and/or MX_Library version found - skipping MEX-File compilation!" )
  return()
endif()

# mexFigtree
add_library(mexFigtree SHARED 
    mexFigtree.cpp
    ../figtree.cpp
    ../KCenterClustering.cpp
    ../../include/figtree.h
    ../KCenterClustering.h)

add_dependencies(mexFigtree figtree)

set_target_properties( mexFigtree
	PROPERTIES
		PREFIX ""
		SUFFIX .${Matlab_MEX_EXTENSION}
		RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        OUTPUT_NAME "figtree"
)

target_include_directories(mexFigtree PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${ANN_dir}/include
    ${Matlab_INCLUDE_DIRS})

target_link_libraries(mexFigtree
    ann_figtree_version
    ${Matlab_LIBRARIES})

# mexFigtreeChooseEvaluationMethod
add_library(mexFigtreeChooseEvaluationMethod SHARED 
    mexFigtreeChooseEvaluationMethod.cpp
    ../figtree.cpp
    ../KCenterClustering.cpp
    ../../include/figtree.h
    ../KCenterClustering.h)

add_dependencies(mexFigtreeChooseEvaluationMethod figtree)

set_target_properties( mexFigtreeChooseEvaluationMethod
	PROPERTIES
		PREFIX ""
		SUFFIX .${Matlab_MEX_EXTENSION}
		RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        OUTPUT_NAME "figtreeChooseEvaluationMethod"
)

target_include_directories(mexFigtreeChooseEvaluationMethod PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${ANN_dir}/include
    ${Matlab_INCLUDE_DIRS})

target_link_libraries(mexFigtreeChooseEvaluationMethod
    ann_figtree_version
    ${Matlab_LIBRARIES})

# mexFigtreeKCenterClustering
add_library(mexFigtreeKCenterClustering SHARED 
    mexFigtreeKCenterClustering.cpp
    ../figtree.cpp
    ../KCenterClustering.cpp
    ../../include/figtree.h
    ../KCenterClustering.h)
    
add_dependencies(mexFigtreeKCenterClustering figtree)
    
set_target_properties( mexFigtreeKCenterClustering
	PROPERTIES
		PREFIX ""
		SUFFIX .${Matlab_MEX_EXTENSION}
		RUNTIME_OUTPUT_DIRECTORY "${PROJECT_BINARY_DIR}"
        OUTPUT_NAME "figtreeKCenterClustering"
)

target_include_directories(mexFigtreeKCenterClustering PRIVATE
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
    ${ANN_dir}/include
    ${Matlab_INCLUDE_DIRS})

target_link_libraries(mexFigtreeKCenterClustering
    ann_figtree_version
    ${Matlab_LIBRARIES})

# export mexFunction    
if(WIN32)
  SET(CMAKE_SHARED_LINKER_FLAGS "/export:mexFunction")
endif()