all: compila clean

compila: esercizio3.3.o lib.o
	g++ esercizio3.3.o lib.o -o esercizio3.3

esercizio3.3.o: esercizio3.3.cpp
	g++ esercizio3.3.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
