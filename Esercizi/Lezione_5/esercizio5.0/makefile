all: compila clean

compila: esercizio5.0.o lib.o
	g++ esercizio5.0.o lib.o -o esercizio5.0

esercizio5.0.o: esercizio5.0.cpp
	g++ esercizio5.0.cpp -c

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
