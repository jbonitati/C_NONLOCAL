#Author: Joey Bonitati
# This makefile finds all of the ".cpp" files recursively in its directory
# and compiles them into the executable specified by the appname variable

appname := run

CXX := g++
LDLIBS := -lgsl -lgslcblas -larmadillo -lboost_system -lboost_timer 
DFLAGS := -lboost_unit_test_framework
LDFlAGS := 
CXXFLAGS :=  -Wall -g -fopenmp

srcfiles := OpticalPotential.cpp system.cpp channel.cpp particle.cpp lagbasis.cpp potential.cpp wavefunction.cpp
mainfiles = $(srcfiles)
mainfiles += multichannel.cpp
testfiles = $(srcfiles)
testfiles += systemtest.cpp
#$(shell find . -name "*.cpp")
mainobjects := $(patsubst %.cpp, %.o, $(mainfiles))
testobjects := $(patsubst %.cpp, %.o, $(testfiles))

all: $(appname)

$(appname): $(mainobjects)
	$(CXX) $(CXXFLAGS) -o $(appname) -O3 $(mainobjects) $(LDFLAGS) $(LDLIBS) 

test: $(testobjects)
	$(CXX) $(CXXFLAGS) -o test -O3 $(testobjects) $(LDFLAGS) $(LDLIBS) $(DFLAGS)

clean:
	$(RM) *.o
	
