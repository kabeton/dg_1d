CXXFLAGS=-g -std=c++17 -Wall

all: solver

solver: second.o dgm.o
	g++ -o dg second.o dgm.o

second.o: second.cpp dgm.h quad.h
	g++ -c second.cpp

dgm.o: dgm.cpp dgm.h quad.h
	g++ -c dgm.cpp

clean:
	rm *.o