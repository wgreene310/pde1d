# UNAME_S contains platform information
UNAME_S:=$(shell uname -s)

# Eigen-related directories
EIGEN_ROOT?=/usr/include/eigen3/
EIGEN_INC:=-I${EIGEN_ROOT}/include/

# Sundials-related directories
SUNDIALS_ROOT?=/usr/
SUNDIALS_LDIR=${SUNDIALS_ROOT}/lib/
SUNDIALS_INC=-I${SUNDIALS_ROOT}/include/

# TODO: Add suitesparse if necessary
# TODO: Add Boost if necessary
# Search paths for include files
INC:=$(EIGEN_INC) $(SUNDIALS_INC) \
     -Ipde1dlib -IFDJacobian -Iutil

# Flags to set for compilation step
CXXFLAGS:=-g -O2 -std=gnu++0x -Wno-deprecated-declarations ${INC}

# Libraries for link step
LIBS:=-l:${SUNDIALS_LDIR}/libsundials_ida.a \
			-l:${SUNDIALS_LDIR}/libsundials_nvecserial.a


# Set USE_OCTAVE=true to override build script
ifneq (${USE_OCTAVE},)
# Octave-related directories
OCTAVE_ROOT?=/usr/include/octave-4.0.0/octave/
OCTAVE_INC:=-I${OCTAVE_ROOT}
OCTAVE_LDIR:=/usr/lib/x86_64-linux-gnu/

INC+=${OCTAVE_INC}
LIBS+=-L${OCTAVE_LDIR} -loctave -loctinterp

override CC:=mkoctfile -mex
override CXX:=mkoctfile -mex
override CXXFLAGS:=${INC}
LDFLAGS?=
else
LDFLAGS+=-shared
endif

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
	@echo "    - SUNDIALS_ROOT: Adjust to root path for Sundials includes/libraries."
	@echo "    - EIGEN_ROOT: Adjust to root path for Eigen includes."

# Files under pde1dlib we expect to have built
PDELIBFILES:=	PDE1dImpl \
							PDE1dDefn \
							PDEInitConditions \
							GausLegendreIntRule \
							SunVector \
							PDEMeshMapper \
							PDEModel \
							PDEElement \
							PDESolution \
							ShapeFunctionManager \
							ShapeFunctionHierarchical \
							ShapeFunction

# All built objects
OBJS:=$(addprefix FDJacobian/,FDJacobian.o FiniteDiffJacobian.o) \
			$(addprefix pde1dmex/,pde1dmex.o PDE1dMexInt.o MexInterface.o) \
			$(addprefix util/,util.o) \
			$(addprefix pde1dlib/,$(addsuffix .o,${PDELIBFILES}))

# Scrape active subdirectories from OBJS
SUBDIRS:=$(subst /,,$(sort $(dir ${OBJS})))

# Build instructions for object files
define BUILD_DIR_template =
$(1)/%.o: $(1)/%.c
	$${CC} $${CXXFLAGS} -o $$@ $$<
endef

# end build template
# Now, instantiate that template for all directories
$(foreach dir,${SUBDIRS},$(eval $(call BUILD_DIR_template,${dir})))
# End build instructions for object files

objects: ${OBJS}

pde1d.mex:
	$(CXX) ${CXXFLAGS} ${LDFLAGS} -o $@ $^ ${LIBS}

clean:
	$(RM) $(foreach dir,${SUBDIRS},${dir}/*.o)
