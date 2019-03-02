## Start PDE1D build configuration
# Eigen-related directories
EIGEN_ROOT?=/usr/include/eigen3/
EIGEN_INC?=-I${EIGEN_ROOT}

# Sundials-related directories
SUNDIALS_ROOT?=/usr/
SUNDIALS_LDIR=${SUNDIALS_ROOT}/lib/
SUNDIALS_INC=-I${SUNDIALS_ROOT}/include/

# Suitesparse-related directories
SUITESPARSE_ROOT?=/usr/include/suitesparse
SUITESPARSE_LDIR?=/usr/lib/x86_64-linux-gnu
SUITESPARSE_INC?=-I${SUITESPARSE_ROOT}

# TODO: Add Boost if necessary
# Search paths for include files
INC=${EIGEN_INC} ${SUNDIALS_INC} ${SUITESPARSE_INC} \
    -Ipde1dlib -IFDJacobian -Iutil

# Flags for compilation step. NOTE: Build SUNDIALS with int32_t
# index type, which will then match with the Eigen index type
# being set here.

CPPFLAGS=-g -Wno-deprecated-declarations -DEIGEN_DEFAULT_DENSE_INDEX_TYPE=int32_t

# Flags for link step
LDFLAGS?=

# Libraries for link step
LDLIBS=-L${SUNDIALS_LDIR} \
			 -lsundials_ida \
			 -lsundials_nvecserial \
			 -lsundials_sunlinsolklu \
			 -L${SUITESPARSE_LDIR} \
			 -lklu

# Set USE_OCTAVE=true to compile against GNU Octave
ifneq (${USE_OCTAVE},) # If USE_OCTAVE is Nonempty
# Octave-related directories
OCTAVE_ROOT?=/usr/include/octave-4.0.0/octave/
OCTAVE_INC?=-I${OCTAVE_ROOT}
OCTAVE_LDIR?=/usr/lib/x86_64-linux-gnu/

INC+=${OCTAVE_INC}
LDLIBS+=-L${OCTAVE_LDIR} -loctave -loctinterp
else # End USE_OCTAVE Nonempty
CPPFLAGS+= -O2 -std=gnu++0x
LDFLAGS+= -shared
endif # End USE_OCTAVE Empty

CPPFLAGS+=${INC}
## End PDE1D Build Configuration

##
# Usage target
##
usage:
	@echo "make TARGET [OPTIONS]"
	@echo ""
	@echo " Valid Targets:"
	@echo "    - pde1d.mex: Build main PDE solver Mex-file"
	@echo "    - clean: Remove object files."
	@echo " Options/Variables"
	@echo "    - USE_OCTAVE: Set to true to build against Octave."
	@echo "    - OCTAVE_ROOT: Adjust to Octave root directory when USE_OCTAVE=true"
	@echo "    - OCTAVE_INC: Command line needed to compile with Octave include files"
	@echo "    - OCTAVE_LDIR: Path to Octave library files"
	@echo "    - SUNDIALS_ROOT: Adjust to root path for Sundials includes/libraries."
	@echo "    - SUNDIALS_INC: Command line needed to compile with Sundials include files"
	@echo "    - SUNDIALS_LDIR: Path to Sundials library files"
	@echo "    - SUITESPARSE_ROOT: Adjust to root path for Suitesparse includes/libraries."
	@echo "    - SUITESPARSE_INC: Command line needed to compile with Suitesparse include files"
	@echo "    - SUITESPARSE_LDIR: Path to Suitesparse library files"
	@echo "    - EIGEN_ROOT: Adjust to root path for Eigen includes."
	@echo "    - EIGEN_INC: Command line needed to compile with Eigen include files"

# Files under pde1dlib we expect to have built
PDELIBFILES=PDE1dImpl \
						PDE1dDefn \
						PDEInitConditions \
						GausLegendreIntRule \
						SunVector \
						PDEMeshMapper \
						PDEModel \
						PDEElement \
						PDEEvents \
						PDESolution \
						ShapeFunctionManager \
						ShapeFunctionHierarchical \
						ShapeFunction

# All built objects
OBJS=$(addprefix FDJacobian/,FDJacobian.o FiniteDiffJacobian.o) \
		 $(addprefix pde1dmex/,pde1dmex.o PDE1dMexInt.o MexInterface.o) \
		 $(addprefix util/,util.o) \
		 $(addprefix pde1dlib/,$(addsuffix .o,${PDELIBFILES}))

# Scrape active subdirectories from OBJS
SUBDIRS=$(subst /,,$(sort $(dir ${OBJS})))

objects: ${OBJS}

# Define build instructions for object files and final mex product
# separately. The instructions differ when building against Octave.
ifneq (${USE_OCTAVE},)
%.o: %.c
	mkoctfile --mex ${CPPFLAGS} -c -o $@ $<

%.o: %.cpp
	mkoctfile --mex ${CPPFLAGS} -c -o $@ $<

pde1d.mex: ${OBJS}
	mkoctfile --mex ${LDFLAGS} -o $@ $^ ${LDLIBS}

else # end USE_OCTAVE Nonempty
%.o: %.c
	$(CXX) ${CPPFLAGS} -c -o $@ $<

%.o: %.cpp
	$(CXX) ${CPPFLAGS} -c -o $@ $<

pde1d.mex: ${OBJS}
	$(CXX) ${LDFLAGS} -o $@ $^ ${LDLIBS}

endif # end USE_OCTAVE Empty
# End build instructions for object files

clean:
	$(RM) $(foreach dir,${SUBDIRS},${dir}/*.o)

clean-all: clean
	$(RM) pde1d.mex
