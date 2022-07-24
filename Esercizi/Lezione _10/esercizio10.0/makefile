LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio10.0.o lib.o
	g++ esercizio10.0.o lib.o -o esercizio10.0 ${INCS} ${LIBS}

esercizio10.0.o: esercizio10.0.cpp
	g++ esercizio10.0.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
