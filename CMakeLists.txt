cmake_minimum_required (VERSION 3.0)

project(pde1d)
SET( CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS TRUE )
cmake_policy(SET CMP0008 NEW)
cmake_policy(SET CMP0074 NEW)

set(PSE "Octave" CACHE STRING "Problem solving environment")
set_property(CACHE PSE PROPERTY STRINGS Octave MATLAB)

#message(STATUS "CMAKE_CXX_COMPILER_ID=" ${CMAKE_CXX_COMPILER_ID})
if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
add_compile_options("/wd4996")
add_compile_options("/wd4477")
add_definitions(-DNOMINMAX)
elseif(${CMAKE_CXX_COMPILER_ID} STREQUAL GNU)
add_compile_options("-std=gnu++0x")
#add_compile_options("-H")
set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# Want to find Eigen and check version even if it hasn't been formally "installed".
# In Eigen 3.3 this requires find_package in module mode for the un-installed case and
# config mode for the installed case.
if(EIGEN_ROOT)
  set(CMAKE_MODULE_PATH ${EIGEN_ROOT} ; ${EIGEN_ROOT}/cmake)
	set(CMAKE_PREFIX_PATH ${EIGEN_ROOT}/share) # config mode
	set(EIGEN3_INCLUDE_DIR  ${EIGEN_ROOT})
	set(Eigen3_FIND_VERSION 3.3)
endif() # EIGEN_ROOT
find_package(Eigen3 ${Eigen3_FIND_VERSION} REQUIRED)
if(NOT EIGEN3_FOUND)
  message(SEND_ERROR "Eigen matrix library not found.")
endif()
# Eigen 3.3 FindEigen3.cmake not issuing error for wrong version
if(EIGEN3_VERSION AND EIGEN3_VERSION VERSION_LESS Eigen3_FIND_VERSION)
  message(SEND_ERROR "Eigen version must be at least " ${Eigen3_FIND_VERSION})
endif()
message(STATUS "EIGEN3_INCLUDE_DIR=" ${EIGEN3_INCLUDE_DIR})
include_directories(${EIGEN3_INCLUDE_DIR})

set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

set(LINEAR_SOLVER "Eigen" CACHE STRING "Sparse Linear Solver")
set_property(CACHE LINEAR_SOLVER PROPERTY STRINGS Eigen KLU)
message(STATUS "LINEAR_SOLVER=" ${LINEAR_SOLVER})

string(TOLOWER ${LINEAR_SOLVER} LSL)
if(${LSL} STREQUAL klu)
	# suitesparse to get KLU for sundials ida
	set(SUITESPARSE_COMP_NAMES suitesparseconfig; KLU; AMD; COLAMD; BTF)
	#set(SuiteSparse_VERBOSE ON)
  set(SUITESPARSE_ROOT SUITESPARSE_ROOT-NOTFOUND CACHE PATH "Path to SuiteSparse solver")
	set(SuiteSparse_DIR ${SUITESPARSE_ROOT})
	find_package(SuiteSparse REQUIRED 
	COMPONENTS ${SUITESPARSE_COMP_NAMES})
	foreach(cn IN LISTS SUITESPARSE_COMP_NAMES)
	  string(TOUPPER ${cn} cnUC)
		include_directories(${SuiteSparse_${cnUC}_INCLUDE_DIR})
		list(APPEND SUITESPARSE_LIBS_RELEASE ${SuiteSparse_${cnUC}_LIBRARY_RELEASE})
		list(APPEND SUITESPARSE_LIBS_DEBUG ${SuiteSparse_${cnUC}_LIBRARY_DEBUG})
	endforeach()
	set(USE_KLU TRUE)
elseif(${LSL} STREQUAL eigen)
  add_definitions(-DUSE_EIGEN_LU)
else()
  message(SEND_ERROR "Unsupported linear solver: " ${LINEAR_SOLVER})
endif() # LINEAR_SOLVER

#
# find sundials includes and libs
#
set(SUNDIALS_LIBS_NAMES sundials_ida sundials_nvecserial sundials_sunmatrixsparse)
if(USE_KLU)
  set(SUNDIALS_LIBS_NAMES ${SUNDIALS_LIBS_NAMES} sundials_sunlinsolklu)
endif()
find_package(SUNDIALS 3.1 REQUIRED COMPONENTS ${SUNDIALS_LIBS_NAMES})
include_directories(${SUNDIALS_INCLUDE_DIR}/../)

find_package( Boost 1.55 REQUIRED )
include_directories(${Boost_INCLUDE_DIRS})

string(TOLOWER ${PSE} PSEL)
if(${PSEL} STREQUAL octave)
	find_package(Octave 5.5 REQUIRED COMPONENTS Development)
	include_directories(${Octave_INCLUDE_DIRS})
	set(MEX_LIBS ${Octave_OCTAVE_LIBRARY} ${Octave_INTERP_LIBRARY})
	message(STATUS "Path to Octave include files: " ${Octave_INCLUDE_DIRS})
elseif(${PSEL} STREQUAL matlab)
	find_package( Matlab REQUIRED COMPONENTS MX_LIBRARY)
	set(MATLAB_ROOT MATLAB_ROOT-NOTFOUND CACHE PATH "Path to MATLAB")
	if(NOT MATLAB_ROOT)
		message(SEND_ERROR "Path to MATLAB not found")
	endif()
	include_directories(${Matlab_INCLUDE_DIRS})
	set(MEX_LIBS 
	${Matlab_MEX_LIBRARY} 
	${Matlab_MX_LIBRARY}
	)
else()
  message(SEND_ERROR "Invalid PSE: " ${PSE})
endif()

FILE(GLOB PDE_LIB_SRC pde1dlib/*.cpp)
FILE(GLOB PDE_H_FILES pde1dlib/*.h)
FILE(GLOB PDE_MEX_SRC pde1dmex/*.cpp)
FILE(GLOB PDE_MEX_H_FILES pde1dmex/*.h)
FILE(GLOB PDE_JAC_SRC FDJacobian/*.cpp FDJacobian/FDJacobian.c)
FILE(GLOB PDE_JAC_H_FILES FDJacobian/*.h FDJacobian/FDJacobian.h)

include_directories(./util)
include_directories(./pde1dlib)
include_directories(./FDJacobian)

add_library(pde1dLib 
${PDE_LIB_SRC} 
${PDE_H_FILES}
${PDE_JAC_SRC}
${PDE_JAC_H_FILES}
util/util.h
util/util.cpp
)

add_library (pde1d SHARED 
${PDE_MEX_SRC}
${PDE_MEX_H_FILES}
)

#target_compile_features(pde1d PUBLIC cxx_std_11)

if(${CMAKE_CXX_COMPILER_ID} STREQUAL MSVC)
target_sources (pde1d PUBLIC ${CMAKE_SOURCE_DIR}/mexFunction.def)
endif()

target_link_libraries(pde1d PRIVATE
  pde1dLib
  debug "${SUNDIALS_LIBS_DEBUG}"
  optimized "${SUNDIALS_LIBS_RELEASE}"
  ${MEX_LIBS}
)

if(USE_KLU)
  target_link_libraries(pde1d PRIVATE 
	debug "${SUITESPARSE_LIBS_DEBUG}"
	optimized "${SUITESPARSE_LIBS_RELEASE}"
	)
endif()

string(TOLOWER ${PSE} PSEL)
if(${PSEL} STREQUAL matlab)
  set(mexext ${Matlab_MEX_EXTENSION})
else()
  set(mexext mex)
endif()

add_custom_command(TARGET pde1d POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different # which executes "cmake - E copy_if_different..."      
				$<TARGET_FILE:pde1d>
        "${CMAKE_INSTALL_PREFIX}/pde1d.${mexext}") 
				
add_subdirectory (tests)