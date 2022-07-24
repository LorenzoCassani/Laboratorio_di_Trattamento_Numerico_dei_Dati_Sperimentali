LIBS:=`root-config --libs`
INCS:=`root-config --cflags`

all: compila clean

compila: compito.o lib.o
	g++ compito.o lib.o -o compito ${INCS} ${LIBS}

compito.o: lib.h compito.cpp
	g++ compito.cpp -c ${INCS}

lib.o: lib.h lib.cpp 
	g++ lib.cpp -c

clean:
	rm *.o 
