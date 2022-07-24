all: compila clean

compila: esercizio6.2.o lib.o
	g++ esercizio6.2.o lib.o -o esercizio6.2

esercizio6.2.o: esercizio6.2.cpp
	g++ esercizio6.2.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
