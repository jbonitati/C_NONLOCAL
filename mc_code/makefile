#Author: Joey Bonitati
# This makefile finds all of the ".cpp" files recursively in its directory
# and compiles them into the executable specified by the appname variable

appname := run

CXX := g++
LDLIBS := -lgsl -lgslcblas -larmadillo -lboost_system -lboost_timer
LDFlAGS := 
CXXFLAGS :=  -Wall -g -fopenmp

srcfiles := $(shell find . -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) -o $(appname) -O3 $(objects) $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o
