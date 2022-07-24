all: compila clean

compila: esercizio6.3.o lib.o
	g++ esercizio6.3.o lib.o -o esercizio6.3

esercizio6.3.o: esercizio6.3.cpp
	g++ esercizio6.3.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
