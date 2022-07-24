all: compila clean

compila: esercizio6.0.o lib.o
	g++ esercizio6.0.o lib.o -o esercizio6.0

esercizio6.0.o: esercizio6.0.cpp
	g++ esercizio6.0.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
