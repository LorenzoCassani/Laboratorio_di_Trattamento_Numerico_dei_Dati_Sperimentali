all: compila clean

compila: esercizio4.2.o lib.o
	g++ esercizio4.2.o lib.o -o esercizio4.2

esercizio4.2.o: esercizio4.2.cpp
	g++ esercizio4.2.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
