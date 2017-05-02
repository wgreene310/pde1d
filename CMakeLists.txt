cmake_minimum_required (VERSION 2.6)

project(pde1d)
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 19)
 set(VC_VER 14)
else()
 set(VC_VER 12)
endif()
add_compile_options("/wd4996")
endif()

set(EIGEN_ROOT $ENV{HOME}/src/Eigen/eigen-3.3)
include_directories(${EIGEN_ROOT})

# suitesparse to get KLU for sundials ida
set(SUITESPARSE_LIB_DIR $ENV{HOME}/src/sparse_solvers/SuiteSparse/libs/x64)
set(SUITESPARSE_LIBS 
libklu.lib
libamd.lib
libcolamd.lib
libbtf.lib
suitesparseconfig.lib)

# sundials ida
set(SUNDIALS_ROOT C:/temp/sundials_2.7_vs2015_x64)
set(SUNDIALS_LIB_DIR ${SUNDIALS_ROOT}/lib)
set(SUNDIALS_LIBS
sundials_ida.lib
sundials_kinsol.lib
sundials_nvecserial.lib)
include_directories(${SUNDIALS_ROOT}/include)

set(BOOST_ROOT c:/boost)
set(Boost_FIND_VERSION_EXACT ON)
set(Boost_LIBRARY_DIR ${BOOST_ROOT}/boost_1_63_0/lib64-msvc-${VC_VER}.0)
find_package( Boost 1.63.0 EXACT REQUIRED )
# COMPONENTS python
#include_directories(${Boost_INCLUDE_DIRS}/boost/tr1/tr1)
include_directories(${Boost_INCLUDE_DIRS})

set(OCTAVE_REL 4.2.0)
set(OCTAVE_ROOT C:/Octave/Octave-${OCTAVE_REL})
include_directories($(OCTAVE_ROOT)/include/octave-$(REL)/octave)

FILE(GLOB PDE_LIB_SRC pde1dlib/*.cpp)
FILE(GLOB PDE_H_FILES pde1dlib/*.h)
FILE(GLOB PDE_MEX_SRC pde1dmex/*.cpp)
FILE(GLOB PDE_MEX_H_FILES pde1dmex/*.h)
FILE(GLOB PDE_JAC_SRC FDJacobian/*.cpp)
FILE(GLOB PDE_JAC_H_FILES FDJacobian/*.h)

#message("PDE_LIB_SRC=" ${PDE_LIB_SRC})

link_directories(${SUNDIALS_LIB_DIR})

include_directories(./util)
include_directories(./pde1dlib)
include_directories(./FDJacobian)

add_library (pde1d SHARED 
${PDE_LIB_SRC} 
${PDE_H_FILES}
${PDE_MEX_SRC}
${PDE_MEX_H_FILES}
${PDE_JAC_SRC}
${PDE_JAC_H_FILES}
util/util.h
util/util.cpp
)

target_link_libraries(pde1d PRIVATE
${SUNDIALS_LIBS}
)