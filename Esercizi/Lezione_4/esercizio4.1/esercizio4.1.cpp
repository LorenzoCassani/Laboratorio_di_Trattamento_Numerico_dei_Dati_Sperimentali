#include "lib.h"

using namespace std;

int main(int argc, char** argv){

    	if(argc!=5){
        	cerr << "Usage: " << argv[0] << " <x> <y> <z> <n> " << endl;
        	exit(-1);
	}

	double x=atof(argv[1]);
	double y=atof(argv[2]);
	double z=atof(argv[3]);
	Posizione r(x,y,z);

	int n=atof(argv[4]); //Ordine del multipolo

	if(n%2!=0){
		cerr << endl << "L'ordine n del multipolo deve essere pari!" << endl << endl;
		exit(-1);
	}

	const double r_0 = 1.E-10; //Raggio del multipolo

	const double e = 1.60217653E-19; //Carica elettrica

	CampoVettoriale E(r);

	for(int i=0; i<n; i++){

		PuntoMateriale particella(0.,pow(-1,i)*e,r_0*cos(2*M_PI*i/n),r_0*sin(2*M_PI*i/n),0.);
		E.Somma(particella.CampoElettrico(r));
	}

	cout << endl << "Campo elettrico generato dal multipolo:" << endl;

	cout << endl << "E(" << x << "," << y << "," << z << ") = ("
	<< E.GetFX() << ","
	<< E.GetFY() << ","
	<< E.GetFZ() << ")"
	<< " V/m" << endl << endl;

	cout << "|E(" << x << "," << y << "," << z << ")| = "
	<< E.Modulo() << " V/m" << endl;	

	cout << endl;

	return 0;
}


