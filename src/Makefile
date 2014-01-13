# Makefile for respect

#G++FLAG = -g
G++FLAG = -O3

respect : respect.o eps_primitives.o
	g++ -o respect respect.o eps_primitives.o -I/usr/local/include -L/usr/local/lib -lgsl -lgslcblas -lm -lpthread

respect.o : respect.cpp
	g++ $(G++FLAG) -c respect.cpp -I/opt/local/include

eps_primitives.o: eps_primitives.cpp
	g++ $(G++FLAG) -c eps_primitives.cpp -I/opt/local/include	
clean:
	rm *.o
