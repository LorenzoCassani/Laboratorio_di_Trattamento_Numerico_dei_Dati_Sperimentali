#include "lib.h"

using namespace std;

int main(int argc, char** argv){

    	if(argc!=4){
        	cerr << "Usage: " << argv[0] << " <x> <y> <z>" << endl;
        	exit(-1);
	}

	double x=atof(argv[1]);
	double y=atof(argv[2]);
	double z=atof(argv[3]);
	Posizione r(x,y,z);

	const double e =1.60217653E-19;
	const double me=9.1093826E-31;
    	const double mp=1.6726217E-27;
	const double d =1.E-10; //Angstrom

	PuntoMateriale elettrone(me,-e,0.,0.,d/2.);
	PuntoMateriale protone	(mp, e,0.,0.,-d/2.);

	cout << endl;

	CampoVettoriale E(r);

	E.Somma(elettrone.CampoElettrico(r));
	E.Somma(protone.CampoElettrico(r));

	cout << "Campo elettrico generato dal dipolo:" << endl;

	cout << endl << "E(" << x << "," << y << "," << z << ") = ("
	<< E.GetFX() << ","
	<< E.GetFY() << ","
	<< E.GetFZ() << ")"
	<< " V/m" << endl << endl;

	cout << "|E(" << x << "," << y << "," << z << ")| = "
	<< E.Modulo() << " V/m" << endl;

	cout << endl << endl;

	CampoVettoriale g(r);

	g.Somma(elettrone.CampoGravitazionale(r));
	g.Somma(protone.CampoGravitazionale(r));

	cout << "Campo gravitazionale generato dal dipolo:" << endl;

	cout << endl << "g(" << x << "," << y << "," << z << ") = ("
	<< g.GetFX() << ","
	<< g.GetFY() << ","
	<< g.GetFZ() << ")"
	<< " m/s^2" << endl << endl;

	cout << "|g(" << x << "," << y << "," << z << ")| = "
	<< g.Modulo() << " m/s^2" << endl;

	cout << endl;

	return 0;
}


