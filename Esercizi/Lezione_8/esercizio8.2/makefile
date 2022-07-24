LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio8.2.o lib.o
	g++ esercizio8.2.o lib.o -o esercizio8.2 ${INCS} ${LIBS}

esercizio8.2.o: esercizio8.2.cpp
	g++ esercizio8.2.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
