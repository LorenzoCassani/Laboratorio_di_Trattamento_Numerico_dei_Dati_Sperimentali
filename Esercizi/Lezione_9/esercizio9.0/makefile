LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio9.0.o lib.o
	g++ esercizio9.0.o lib.o -o esercizio9.0 ${INCS} ${LIBS}

esercizio9.0.o: esercizio9.0.cpp
	g++ esercizio9.0.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
