all: compila clean

compila: esercizio5.3.o lib.o
	g++ esercizio5.3.o lib.o -o esercizio5.3

esercizio5.3.o: esercizio5.3.cpp
	g++ esercizio5.3.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
