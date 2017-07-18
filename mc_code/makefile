appname := run

CXX := g++
LDLIBS := -lgsl -lgslcblas -larmadillo -lboost_system -lboost_timer
CXXFLAGS :=  -Wall -g -fopenmp

srcfiles := $(shell find . -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) -O3 $(objects) $(LDLIBS)

clean:
	$(RM) *.o
