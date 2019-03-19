
add_executable (testPde1d testPde1d.cpp 
PDE1dTestDefn.h
PDE1dTestDefn.cpp
ExampleHeatCond.h)

#message("SUNDIALS_LIBS_RELEASE=" "${SUNDIALS_LIBS_RELEASE}")
target_link_libraries(testPde1d PRIVATE
  pde1dLib
  debug "${SUNDIALS_LIBS_DEBUG}"
  optimized "${SUNDIALS_LIBS_RELEASE}"
)

if(USE_KLU)
target_link_libraries(testPde1d PRIVATE
	debug "${SUITESPARSE_LIBS_DEBUG}"
	optimized "${SUITESPARSE_LIBS_RELEASE}"
)
endif()