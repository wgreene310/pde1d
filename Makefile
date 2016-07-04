CC=gcc
CXX=g++

EIGEN=$(HOME)/src/Eigen/eigen-3.2.2
OCTAVE_ROOT=/opt/octave-4.0.0
OCTAVE_INC=$(OCTAVE_ROOT)/include/octave-4.0.2/octave
SUNDIALS_ROOT=/usr/local/include/sundials
SUNDIALS_LDIR=$(SUNDIALS_ROOT)/lib/Release
SUNDIALS_INC=-I$(SUNDIALS_ROOT)/include -I./ 

INC=-I$(EIGEN) -I$(OCTAVE_INC) $(SUNDIALS_INC) -I../pde1dlib

CXXFLAGS= -g -O2 -std=gnu++0x -Wno-deprecated-declarations $(INC)


VPATH=../pde1dlib ../pde1dmex 

LIBS=-L$(OCTAVE_ROOT)/lib/octave/4.0.2 -loctave -loctinterp \
      -l:libsundials_ida.a -l:libsundials_nvecserial.a

OBJS=pde1dmex.o PDE1dMexInt.o PDE1dImpl.o PDE1dDefn.o \
     GausLegendreIntRule.o SunVector.o

pde1d.mex: $(OBJS)
	$(CXX) -shared -o pde1d.mex $(OBJS) $(LIBS)
	
clean:
	rm -f *.o

	
