all: solver

solver: second.o dgm.o
	g++ -o dg second.o dgm.o

second.o: second.cpp dgm.h
	g++ -c second.cpp

dgm.o: dgm.cpp dgm.h
	g++ -c dgm.cpp

clean:
	rm *.o