LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio10.2.o lib.o
	g++ esercizio10.2.o lib.o -o esercizio10.2 ${INCS} ${LIBS}

esercizio10.2.o: esercizio10.2.cpp
	g++ esercizio10.2.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
