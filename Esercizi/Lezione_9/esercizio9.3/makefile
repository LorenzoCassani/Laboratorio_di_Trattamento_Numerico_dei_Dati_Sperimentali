LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio9.3.o lib.o
	g++ esercizio9.3.o lib.o -o esercizio9.3 ${INCS} ${LIBS}

esercizio9.3.o: esercizio9.3.cpp
	g++ esercizio9.3.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
