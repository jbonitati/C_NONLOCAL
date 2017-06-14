appname := run

CXX := g++
LDLIBS := -lgsl -lgslcblas -larmadillo
CXXFLAGS := -Wall -g

srcfiles := $(shell find . -name "*.cpp")
objects  := $(patsubst %.cpp, %.o, $(srcfiles))

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS)

clean:
	$(RM) $(objects)
