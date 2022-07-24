all: compila clean

compila: esercizio1.1.o lib.o
	g++ esercizio1.1.o lib.o -o esercizio1.1

esercizio1.1.o: esercizio1.1.cpp
	g++ esercizio1.1.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
