LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio8.0.o lib.o
	g++ esercizio8.0.o lib.o -o esercizio8.0 ${INCS} ${LIBS}

esercizio8.0.o: esercizio8.0.cpp
	g++ esercizio8.0.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
