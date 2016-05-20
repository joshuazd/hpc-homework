SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
LIB=$(wildcard lib/*.cpp)
LOBJ=$(LIB:.cpp=.o)
CXX=mpicxx
CXXFLAGS=-std=c++0x -Wall -O3 -fopenmp

LDFLAGS=-limf

image_denoise : $(OBJ) $(LOBJ)
	$(CXX) $(LDFLAGS) -o $@ -fopenmp $^

clean:
	rm *.o image_denoise
