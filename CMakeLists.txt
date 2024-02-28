cmake_minimum_required(VERSION 3.24.3)
project(StringNetFrequency)

set(CMAKE_CXX_STANDARD 17)

add_executable(main
        src/main.cpp
        src/experiments.cpp src/experiments.hpp
        src/suffix_arrays.cpp src/suffix_arrays.hpp
        )

include(ExternalProject)

set(EXTERNAL_INSTALL_LOCATION ${CMAKE_BINARY_DIR}/external)

ExternalProject_Add(sdsl
        GIT_REPOSITORY https://github.com/simongog/sdsl-lite
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${EXTERNAL_INSTALL_LOCATION}
        )

include_directories(${EXTERNAL_INSTALL_LOCATION}/include)
link_directories(${EXTERNAL_INSTALL_LOCATION}/lib)

add_dependencies(main sdsl)

target_link_libraries(main ${EXTERNAL_INSTALL_LOCATION}/lib/libsdsl.a)
target_link_libraries(main ${EXTERNAL_INSTALL_LOCATION}/lib/libdivsufsort.a)
target_link_libraries(main ${EXTERNAL_INSTALL_LOCATION}/lib/libdivsufsort64.a)
