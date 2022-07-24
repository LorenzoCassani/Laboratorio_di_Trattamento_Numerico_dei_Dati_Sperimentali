all: compila clean

compila: esercizio7.2.o lib.o
	g++ esercizio7.2.o lib.o -o esercizio7.2

esercizio7.2.o: esercizio7.2.cpp
	g++ esercizio7.2.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
