all: compila clean

compila: esercizio7.1.o lib.o
	g++ esercizio7.1.o lib.o -o esercizio7.1

esercizio7.1.o: esercizio7.1.cpp
	g++ esercizio7.1.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
