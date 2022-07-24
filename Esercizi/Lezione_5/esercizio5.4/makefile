all: compila clean

compila: esercizio5.4.o lib.o
	g++ esercizio5.4.o lib.o -o esercizio5.4

esercizio5.4.o: esercizio5.4.cpp
	g++ esercizio5.4.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
