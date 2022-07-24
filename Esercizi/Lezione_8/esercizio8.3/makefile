LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: esercizio8.3.o lib.o
	g++ esercizio8.3.o lib.o -o esercizio8.3 ${INCS} ${LIBS}

esercizio8.3.o: esercizio8.3.cpp
	g++ esercizio8.3.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
