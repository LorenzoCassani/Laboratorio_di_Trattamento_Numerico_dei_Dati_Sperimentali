LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio9.2.o lib.o
	g++ esercizio9.2.o lib.o -o esercizio9.2 ${INCS} ${LIBS}

esercizio9.2.o: esercizio9.2.cpp
	g++ esercizio9.2.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
