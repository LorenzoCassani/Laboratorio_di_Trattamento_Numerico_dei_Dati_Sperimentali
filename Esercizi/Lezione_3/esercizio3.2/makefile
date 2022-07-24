all: compila clean

compila: esercizio3.2.o lib.o
	g++ esercizio3.2.o lib.o -o esercizio3.2

esercizio3.2.o: esercizio3.2.cpp
	g++ esercizio3.2.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
