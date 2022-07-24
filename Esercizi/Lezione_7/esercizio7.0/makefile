all: compila clean

compila: esercizio7.0.o lib.o
	g++ esercizio7.0.o lib.o -o esercizio7.0

esercizio7.0.o: esercizio7.0.cpp
	g++ esercizio7.0.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
