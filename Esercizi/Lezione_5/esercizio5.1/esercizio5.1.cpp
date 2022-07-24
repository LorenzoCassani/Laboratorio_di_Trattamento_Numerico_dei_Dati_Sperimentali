#include "lib.h"

using namespace std;

int main(int argc, char **argv){

	if(argc!=4){
        	cerr << "Usage: " << argv[0] << " <a> <b> <c>" << endl;
        	exit(-1);
	}

	double a=atof(argv[1]);
	double b=atof(argv[2]);
	double c=atof(argv[3]);
	Parabola p(a,b,c); //Utilizzo il costruttore  della Classe Parabola a partire da a,b,c

    	double x_v=-p.GetB()/(2*p.GetA()); //Calcolo la x del vertice utilizzando i data membri della parabola

	cout << endl << "Equazione parabola: y = " << p.GetA() << "x^2 + " << p.GetB() << "x + " << p.GetC() << endl;
	cout << endl << "x_v = " << x_v << endl;
    	cout << endl << "y_v = " << p.Eval(x_v) << endl << endl; //Valuto la parabola nel suo vertice

    	return 0;
}


