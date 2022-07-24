all: compila clean

compila: esercizio6.1.o lib.o
	g++ esercizio6.1.o lib.o -o esercizio6.1

esercizio6.1.o: esercizio6.1.cpp
	g++ esercizio6.1.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
