appname := run

CXX := g++
LDLIBS := -lgsl -lgslcblas -larmadillo -lboost_system
CXXFLAGS :=  -Wall -g

srcfiles := $(shell find . -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS)

clean:
	$(RM) *.o *~
