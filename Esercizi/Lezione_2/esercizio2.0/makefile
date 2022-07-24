all: compila clean

compila: esercizio2.0.o lib.o
	g++ esercizio2.0.o lib.o -o esercizio2.0

esercizio2.0.o: esercizio2.0.cpp
	g++ esercizio2.0.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
