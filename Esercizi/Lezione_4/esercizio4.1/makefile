all: compila clean

compila: esercizio4.1.o lib.o
	g++ esercizio4.1.o lib.o -o esercizio4.1

esercizio4.1.o: esercizio4.1.cpp
	g++ esercizio4.1.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
