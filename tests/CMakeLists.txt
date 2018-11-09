
add_executable (testPde1d testPde1d.cpp 
PDE1dTestDefn.h
PDE1dTestDefn.cpp
ExampleHeatCond.h)

target_link_libraries(testPde1d PRIVATE
${SUNDIALS_LIBS}
pde1dLib
)

if(${LINEAR_SOLVER} STREQUAL KLU)
target_link_libraries(testPde1d PRIVATE ${SUITESPARSE_LIBS})
endif()