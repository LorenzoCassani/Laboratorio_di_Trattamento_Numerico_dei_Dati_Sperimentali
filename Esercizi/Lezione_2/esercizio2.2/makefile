all: compila clean

compila: esercizio2.2.o lib.o
	g++ esercizio2.2.o lib.o -o esercizio2.2

esercizio2.2.o: esercizio2.2.cpp
	g++ esercizio2.2.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
