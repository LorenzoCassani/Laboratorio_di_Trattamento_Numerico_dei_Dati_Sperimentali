all: compila clean

compila: esercizio3.0.o lib.o
	g++ esercizio3.0.o lib.o -o esercizio3.0

esercizio3.0.o: esercizio3.0.cpp
	g++ esercizio3.0.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
