cmake_minimum_required (VERSION 3.1)

project(pde1d)
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
cmake_policy(SET CMP0008 NEW)

set(PSE "Octave" CACHE STRING "Problem solving environment")
set_property(CACHE PSE PROPERTY STRINGS Octave MATLAB)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 19)
 set(VC_VER 14)
else()
 set(VC_VER 12)
endif()
add_compile_options("/wd4996")
add_definitions(-DNOMINMAX)
endif()

set(EIGEN_ROOT "$ENV{HOME}/src/Eigen/eigen-3.3.4" 
    CACHE FILEPATH "Path to Eigen library")
include_directories(${EIGEN_ROOT})


set(LINEAR_SOLVER "KLU" CACHE STRING "Sparse Linear Solver")
set_property(CACHE LINEAR_SOLVER PROPERTY STRINGS Eigen KLU)

if(${LINEAR_SOLVER} STREQUAL KLU)
# suitesparse to get KLU for sundials ida
set(SUITESPARSE_ROOT "C:/temp/Projects/suitesparse5.1/install")
set(SUITESPARSE_LIB_DIR "${SUITESPARSE_ROOT}/lib64"
    CACHE FILEPATH "Path to Suitesparse libraries"
		)
include_directories(${SUITESPARSE_ROOT}/include/suitesparse)
# 
# With Visual Studio we want different debug/release libraries
# depending on the IDE setting. So they can't be identified explicitly
# at configuration time.
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
set(SUITESPARSE_LIBS 
libklu
libamd
libcolamd
libbtf
suitesparseconfig)
else()
find_library(SS_KLU NAMES klu PATHS ${SUITESPARSE_LIB_DIR})
find_library(SS_AMD NAMES amd PATHS ${SUITESPARSE_LIB_DIR})
find_library(SS_COLAMD NAMES colamd PATHS ${SUITESPARSE_LIB_DIR})
find_library(SS_BTF NAMES btf PATHS ${SUITESPARSE_LIB_DIR})
find_library(SS_CONFIG NAMES suitesparseconfig PATHS ${SUITESPARSE_LIB_DIR})
set(SUITESPARSE_LIBS
${SS_KLU}
${SS_AMD}
${SS_COLAMD}
${SS_BTF}
${SS_CONFIG}
)
endif()
else()
add_definitions(-DUSE_EIGEN_LU)
endif()

# sundials ida
set(SUNDIALS_ROOT "C:/temp/sundials_3.1_vs2015_x64"
    CACHE FILEPATH "Path to Sundials"
)
set(SUNDIALS_LIBS
sundials_ida
sundials_sunlinsolklu.lib
sundials_nvecserial)
set(SUNDIALS_LIB_DIR ${SUNDIALS_ROOT}/lib)
include_directories(${SUNDIALS_ROOT}/include)

find_package( Boost 1.67 REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})

if(${PSE} STREQUAL Octave)
set(OCTAVE_REL 4.4.0)
set(OCTAVE_ROOT C:/Octave/Octave-${OCTAVE_REL})
include_directories(${OCTAVE_ROOT}/include/octave-${OCTAVE_REL}/octave)
set(MEX_LIB_DIR ${OCTAVE_ROOT}/lib/octave/${OCTAVE_REL})
set(MEX_LIBS ${MEX_LIB_DIR}/liboctave.dll.a
${MEX_LIB_DIR}/liboctinterp.dll.a)
else()
set(MATLAB_ROOT "C:/Program Files/MATLAB/R2014a")
include_directories(${MATLAB_ROOT}/extern/include)
set(MEX_LIB_DIR ${MATLAB_ROOT}/extern/lib/win64/microsoft)
set(MEX_LIBS libmex libmx)
endif()

FILE(GLOB PDE_LIB_SRC pde1dlib/*.cpp)
FILE(GLOB PDE_H_FILES pde1dlib/*.h)
FILE(GLOB PDE_MEX_SRC pde1dmex/*.cpp)
FILE(GLOB PDE_MEX_H_FILES pde1dmex/*.h)
FILE(GLOB PDE_JAC_SRC FDJacobian/*.cpp FDJacobian/FDJacobian.c)
FILE(GLOB PDE_JAC_H_FILES FDJacobian/*.h FDJacobian/FDJacobian.h)

#message("PDE_LIB_SRC=" ${PDE_LIB_SRC})

link_directories(${SUNDIALS_LIB_DIR} ${SUITESPARSE_LIB_DIR} ${MEX_LIB_DIR})

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

target_compile_features(pde1d PUBLIC cxx_std_11)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
target_sources (pde1d PUBLIC mexFunction.def)
endif()

target_link_libraries(pde1d PRIVATE
${SUNDIALS_LIBS}
${MEX_LIBS}
)

if(${LINEAR_SOLVER} STREQUAL KLU)
target_link_libraries(pde1d PRIVATE ${SUITESPARSE_LIBS})
endif()

add_custom_command(TARGET pde1d POST_BUILD        # Adds a post-build event to MyTest
    COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."      
				$<TARGET_FILE:pde1d>
        "${CMAKE_INSTALL_PREFIX}/pde1d.mex") 