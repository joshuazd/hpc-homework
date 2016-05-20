SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)
LIB=$(wildcard lib/*.cpp)
LOBJ=$(LIB:.cpp=.o)
CXX=g++
CXXFLAGS=-std=c++0x -Wall -O3 

LDFLAGS=-limf

image_denoise : $(OBJ) $(LOBJ)
	$(CXX) $(LDFLAGS) -o $@ $^

clean:
	rm *.o image_denoise
