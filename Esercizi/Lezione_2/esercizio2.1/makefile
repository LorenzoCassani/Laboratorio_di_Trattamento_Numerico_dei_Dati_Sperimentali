all: compila clean

compila: esercizio2.1.o lib.o
	g++ esercizio2.1.o lib.o -o esercizio2.1

esercizio2.1.o: esercizio2.1.cpp
	g++ esercizio2.1.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
