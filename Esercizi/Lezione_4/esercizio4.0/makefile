all: compila clean

compila: esercizio4.0.o lib.o
	g++ esercizio4.0.o lib.o -o esercizio4.0

esercizio4.0.o: esercizio4.0.cpp
	g++ esercizio4.0.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
