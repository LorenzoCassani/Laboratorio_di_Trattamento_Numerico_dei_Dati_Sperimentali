LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio10.4.o lib.o
	g++ esercizio10.4.o lib.o -o esercizio10.4 ${INCS} ${LIBS}

esercizio10.4.o: esercizio10.4.cpp
	g++ esercizio10.4.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
