LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio10.1.o lib.o
	g++ esercizio10.1.o lib.o -o esercizio10.1 ${INCS} ${LIBS}

esercizio10.1.o: esercizio10.1.cpp
	g++ esercizio10.1.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
