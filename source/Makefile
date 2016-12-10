###################################################################
# This Makefile was created using the CreateProject.sh script.
# CreateProject.sh is part of Bayesian Analysis Toolkit (BAT).
# BAT can be downloaded from http://www.mppmu.mpg.de/bat
###################################################################
#
# Run 'make' to compile the program and 'make clean' to remove
# all compiled parts and 'clean' the directory.
#
# You might need to adjust the CFLAGS and LIBS based on
# the BAT installation on your system. Consult the gmake manual
# for details.
#
###################################################################

# Root variables
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := -lMinuit $(shell root-config --libs)

ROOTLIBS     += -lRooFitCore -lRooFit -lRooStats -lFoam -lMathMore

# compiler and flags
CXX          = g++
CXXFLAGS     =  -g -Wall -fPIC -Wno-deprecated -O2
LD           = /usr/llvm-gcc-4.2/libexec/gcc/i686-apple-darwin11/4.2.1/ld
LDFLAGS      = -g -O2 -L/Users/francesco/Cuba-3.0/lib
SOFLAGS      = -shared

# standard commands
RM           = rm -f
MV           = mv
ECHO         = echo
CINT         = $(ROOTSYS)/bin/rootcint

# add ROOT flags
CXXFLAGS    += $(ROOTCFLAGS) -I/Users/francesco/Cuba-3.0/include

# ----------------------------------------------------------------------
# The following definitions depend on the setup of the system where
# the project is being compiled. If BAT is installed in the standard
# system search path or the installation directory is defined in the
# BATINSTALLDIR environmental variable then the lines below are correct
# and the compilation will work
CXXFLAGS    += -I./include -I$(BATINSTALLDIR)/include
LIBS        += -L$(BATINSTALLDIR)/lib -lBATmodels -lBAT $(ROOTLIBS)

# List of all classes (models) used in the program
# Add classes to the end. Baskslash indicates continuation
# on the next line
#CXXSRCS      =	BinomialModel.cxx
CXXSRCS      =	$(shell ls src/*.cc)

# ----------------------------------------------------------------------
# don't change lines below unless you know what you're doing
#

#CXXOBJS      = $(patsubst %.cxx,%.o,$(CXXSRCS))
CXXOBJS      = $(patsubst %.cc,%.o,$(CXXSRCS))
EXEOBJS      =

MYPROGS      = sampleanalysis_v2.1
#MYPROGS      = GatorBATCalib
#MYPROGS      = GatorFinalCalib
#MYPROGS      = TestMetropolis
#MYPROGS      = GatorProfilePlot
#MYPROGS      = GatorProfileTest
#MYPROGS      = GatorTreeLikelihood
#MYPROGS      = GatorTreePLstat
#MYPROGS      = GatorStatAnalisys


GARBAGE      = $(CXXOBJS) $(EXEOBJS) *.o *~ link.d $(MYPROGS)

# targets
all : project

link.d : $(patsubst %.cc,%.h,$(CXXSRCS))
	$(CXX) -MM $(CXXFLAGS) $(CXXSRCS) > link.d;

-include link.d

%.o : %.cc
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean :
	$(RM) $(GARBAGE)

project : $(MYPROGS).cpp $(CXXOBJS)
	$(CXX) $(CXXFLAGS) -c $<
	$(CXX) $(LDFLAGS) $(MYPROGS).o $(CXXOBJS) $(LIBS) -o $(MYPROGS)

print :
	echo compiler  : $(CXX)
	echo c++ srcs  : $(CXXSRCS)
	echo c++ objs  : $(CXXOBJS)
	echo c++ flags : $(CXXFLAGS)
	echo libs      : $(LIBS)
	echo so flags  : $(SOFLAGS)

	echo rootlibs  : $(ROOTLIBS)
