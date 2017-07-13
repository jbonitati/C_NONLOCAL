appname := run

CXX := g++
LDLIBS := -lgsl -lgslcblas -larmadillo -lboost_system
CXXFLAGS :=  -g

srcfiles := new_wavefuntion.cpp
objects  := new_wavefunction.o

all: $(appname)

$(appname): $(objects)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(appname) $(objects) $(LDLIBS)

clean:
	$(RM) *.o
