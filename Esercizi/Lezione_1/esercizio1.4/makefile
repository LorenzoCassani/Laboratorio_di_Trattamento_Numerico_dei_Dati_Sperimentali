all: compila clean

compila: esercizio1.4.o lib.o
	g++ esercizio1.4.o lib.o -o esercizio1.4

esercizio1.4.o: esercizio1.4.cpp
	g++ esercizio1.4.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
