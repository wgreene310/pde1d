cmake_minimum_required (VERSION 3.0)

project(pde1d)
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
cmake_policy(SET CMP0008 NEW)

set(PSE "Octave" CACHE STRING "Problem solving environment")
set_property(CACHE PSE PROPERTY STRINGS Octave MATLAB)

#message("CMAKE_CXX_COMPILER_ID=" ${CMAKE_CXX_COMPILER_ID})
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
add_compile_options("/wd4996")
add_compile_options("/wd4477")
add_definitions(-DNOMINMAX)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
add_compile_options("-std=gnu++0x")
#add_compile_options("-H")
endif()

find_path(EIGEN_ROOT signature_of_eigen3_matrix_library 
          DOC "Path to Eigen library")
message("EIGEN_ROOT=" ${EIGEN_ROOT})
include_directories(${EIGEN_ROOT})


set(LINEAR_SOLVER "KLU" CACHE STRING "Sparse Linear Solver")
set_property(CACHE LINEAR_SOLVER PROPERTY STRINGS Eigen KLU)
message("LINEAR_SOLVER=" ${LINEAR_SOLVER})

if(${LINEAR_SOLVER} STREQUAL KLU)
# suitesparse to get KLU for sundials ida
set(SUITESPARSE_ROOT SUITESPARSE_ROOT-NOTFOUND CACHE PATH "Path to SuiteSparse")
if(SUITESPARSE_ROOT)
set(SUITESPARSE_LIB_DIR "${SUITESPARSE_ROOT}/lib64"
    CACHE PATH "Path to Suitesparse libraries"
		)
set(SUITESPARSE_INC_DIR "${SUITESPARSE_ROOT}/include/suitesparse"
    CACHE PATH "Path to Suitesparse includes"
		)
endif()
if(NOT SUITESPARSE_LIB_DIR)
  message(SEND_ERROR "Path to SuiteSparse libraries not found")
endif()
if(NOT SUITESPARSE_INC_DIR)
  message(SEND_ERROR "Path to SuiteSparse include files not found")
endif()
message("SUITESPARSE_INC_DIR=" ${SUITESPARSE_INC_DIR})
include_directories(${SUITESPARSE_INC_DIR})
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
set(SUNDIALS_ROOT SUNDIALS_ROOT-NOTFOUND CACHE PATH "Path to Sundials")
if(NOT SUNDIALS_ROOT)
  message(SEND_ERROR "Path to Sundials not found")
endif()
set(SUNDIALS_LIBS
sundials_ida
sundials_sunlinsolklu
sundials_nvecserial)
set(SUNDIALS_LIB_DIR ${SUNDIALS_ROOT}/lib)
include_directories(${SUNDIALS_ROOT}/include)
message("SUNDIALS_ROOT=" ${SUNDIALS_ROOT})

find_package( Boost 1.55 REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})

if(${PSE} STREQUAL Octave)
set(OCTAVE_REL 4.4.0 CACHE STRING "Octave Version")
set(OCTAVE_ROOT OCTAVE_ROOT-NOTFOUND CACHE PATH "Path to Octave")
if(NOT OCTAVE_ROOT)
  message(SEND_ERROR "Path to Octave not found")
endif()
include_directories(${OCTAVE_ROOT}/include/octave-${OCTAVE_REL}/octave)
set(MEX_LIB_DIR ${OCTAVE_ROOT}/lib/octave/${OCTAVE_REL})
if(${CMAKE_SYSTEM_NAME} STREQUAL Windows)
  set(MEX_LIBS ${MEX_LIB_DIR}/liboctave.dll.a
               ${MEX_LIB_DIR}/liboctinterp.dll.a)
else()
  set(MEX_LIBS octave octinterp)
endif()
else()
set(MATLAB_ROOT MATLAB_ROOT-NOTFOUND CACHE PATH "Path to MATLAB")
if(NOT MATLAB_ROOT)
  message(SEND_ERROR "Path to MATLAB not found")
endif()
include_directories(${MATLAB_ROOT}/extern/include)
set(MEX_LIB_DIR ${MATLAB_ROOT}/extern/lib/win64/microsoft)
set(MEX_LIBS libmex libmx)
endif()
message("OCTAVE_ROOT=" ${OCTAVE_ROOT})
message("OCTAVE_REL=" ${OCTAVE_REL})

FILE(GLOB PDE_LIB_SRC pde1dlib/*.cpp)
FILE(GLOB PDE_H_FILES pde1dlib/*.h)
FILE(GLOB PDE_MEX_SRC pde1dmex/*.cpp)
FILE(GLOB PDE_MEX_H_FILES pde1dmex/*.h)
FILE(GLOB PDE_JAC_SRC FDJacobian/*.cpp FDJacobian/FDJacobian.c)
FILE(GLOB PDE_JAC_H_FILES FDJacobian/*.h FDJacobian/FDJacobian.h)

message("MEX_LIB_DIR=" ${MEX_LIB_DIR})

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

#target_compile_features(pde1d PUBLIC cxx_std_11)

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

string(TOLOWER ${PSE} PSEL)
if(${PSEL} STREQUAL matlab)
if(${CMAKE_SYSTEM_NAME} STREQUAL Windows)
set(mexext mexw64)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL Linux)
set(mexext mexa64)
elseif(${CMAKE_SYSTEM_NAME} STREQUAL Darwin)
set(mexext mexmaci64)
endif()
else()
  set(mexext mex)
endif()

add_custom_command(TARGET pde1d POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different # which executes "cmake - E copy_if_different..."      
				$<TARGET_FILE:pde1d>
        "${CMAKE_INSTALL_PREFIX}/pde1d.${mexext}") 