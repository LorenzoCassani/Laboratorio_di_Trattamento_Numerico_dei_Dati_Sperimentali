all: compila clean

compila: esercizio5.2.o lib.o
	g++ esercizio5.2.o lib.o -o esercizio5.2

esercizio5.2.o: esercizio5.2.cpp
	g++ esercizio5.2.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
