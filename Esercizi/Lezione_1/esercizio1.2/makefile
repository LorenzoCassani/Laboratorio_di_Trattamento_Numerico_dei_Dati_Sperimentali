all: compila clean

compila: esercizio1.2.o lib.o
	g++ esercizio1.2.o lib.o -o esercizio1.2

esercizio1.2.o: esercizio1.2.cpp
	g++ esercizio1.2.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
