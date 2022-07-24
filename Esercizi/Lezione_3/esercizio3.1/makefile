all: compila clean

compila: esercizio3.1.o lib.o
	g++ esercizio3.1.o lib.o -o esercizio3.1

esercizio3.1.o: esercizio3.1.cpp
	g++ esercizio3.1.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
