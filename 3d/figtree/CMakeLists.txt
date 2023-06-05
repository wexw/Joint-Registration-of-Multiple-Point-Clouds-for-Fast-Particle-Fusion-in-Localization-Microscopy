
### CMake ###

cmake_minimum_required(VERSION 3.11)

set_property(GLOBAL PROPERTY USE_FOLDERS ON)

project(figtree)

### Preprocessor ###

if(WIN32)
    add_definitions(-DWIN32)
    add_definitions(-DDLL_EXPORTS)
endif()

### Paths ###

set(ANN_dir ${PROJECT_SOURCE_DIR}/external/ann_1.1.1)

### Matlab ###

find_package(Matlab COMPONENTS MX_LIBRARY)

### Subdirectories ###

add_subdirectory(external/ann_1.1.1)
add_subdirectory(src)
add_subdirectory(samples)
add_subdirectory(test)