cmake_minimum_required (VERSION 2.8) 
project (CODDE) 
set (DEOM_VERSION_MAJOR 0) 
set (DEOM_VERSION_MINOR 0) 

set (CMAKE_VERBOSE_MAKEFILE OFF) 

file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/bin) 
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/include) 
file(MAKE_DIRECTORY ${PROJECT_SOURCE_DIR}/lib) 

set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/bin)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib) 
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${PROJECT_SOURCE_DIR}/lib) 

set (CMAKE_CXX_FLAGS "-std=c++11 -g -Wall -Ofast -march=native")

include_directories(${PROJECT_SOURCE_DIR}/include 
                    ${PROJECT_SOURCE_DIR}/src/codde
                   )

link_directories(${PROJECT_SOURCE_DIR}/lib)

add_subdirectory(${PROJECT_SOURCE_DIR}/src/codde)

FILE(GLOB apps_src ${PROJECT_SOURCE_DIR}/src/apps/*.cpp)  

foreach (p ${apps_src})
    get_filename_component(p_bin ${p} NAME_WE)
    add_executable(${p_bin} ${p})
    target_link_libraries(${p_bin} codde json11 armadillo lapack blas gfortran)
endforeach ()

include(ExternalProject)

if (EXISTS ${PROJECT_SOURCE_DIR}/lib/libjson11.a)
else()
  ExternalProject_Add(json11
    GIT_REPOSITORY https://github.com/dropbox/json11.git
    PREFIX ${PROJECT_BINARY_DIR}/thirdparty
    INSTALL_DIR ${PROJECT_SOURCE_DIR}
    INSTALL_COMMAND cp <SOURCE_DIR>/json11.hpp <INSTALL_DIR>/include && cp <BINARY_DIR>/libjson11.a <INSTALL_DIR>/lib
  )
endif()
#if (EXISTS ${PROJECT_SOURCE_DIR}/include/armadillo)
#else()
#  ExternalProject_Add(armadillo
#    URL file:///home/codes/armadillo-9.800.2.tar.xz
#    PREFIX ${PROJECT_BINARY_DIR}/thirdparty
#    INSTALL_DIR ${PROJECT_SOURCE_DIR}
#    CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=<INSTALL_DIR>
#  )
#endif()

