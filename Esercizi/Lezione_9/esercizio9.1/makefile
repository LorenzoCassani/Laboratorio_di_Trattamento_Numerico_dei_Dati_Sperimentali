LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio9.1.o lib.o
	g++ esercizio9.1.o lib.o -o esercizio9.1 ${INCS} ${LIBS}

esercizio9.1.o: esercizio9.1.cpp
	g++ esercizio9.1.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
