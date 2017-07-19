appname := run

CXX := g++
LDLIBS := -lgsl -lgslcblas -larmadillo -lboost_system -lboost_timer
CXXFLAGS :=  -Wall -g -fopenmp

srcfiles := $(shell find . -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) -o $(appname) -O3 $(objects) $(LDFLAGS) $(LDLIBS)

clean:
	$(RM) *.o
