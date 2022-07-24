all: compila clean

compila: esercizio1.3.o lib.o
	g++ esercizio1.3.o lib.o -o esercizio1.3

esercizio1.3.o: esercizio1.3.cpp
	g++ esercizio1.3.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
