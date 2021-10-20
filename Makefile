all: solver

solver: second.o dgm.o
	g++ -g -Wall -o dg second.o dgm.o

second.o: second.cpp dgm.h quad.h
	g++ -g -c second.cpp

dgm.o: dgm.cpp dgm.h quad.h
	g++ -g -c dgm.cpp

clean:
	rm *.o