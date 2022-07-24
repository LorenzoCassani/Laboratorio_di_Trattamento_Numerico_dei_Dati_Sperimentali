all: compila clean

compila: esercizio5.1.o lib.o
	g++ esercizio5.1.o lib.o -o esercizio5.1

esercizio5.1.o: esercizio5.1.cpp
	g++ esercizio5.1.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
