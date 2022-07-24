#include "lib.h"

#define MASSATERRA 5.972E24
#define RAGGIOTERRA 6.371E6

using namespace std;

int main(){

	double h = 2.5E+5;

	Posizione GOCE(0.,0.,RAGGIOTERRA+h); //Posizione relativa del satellite GOCE rispetto al centro della Terra

	PuntoMateriale Terra(MASSATERRA,0.,0.,0.,0.); //Tratto la Terra come un punto materiale posto nell'origine del nostro s.d.r.

	CampoVettoriale g(GOCE);

	g.Somma(Terra.CampoGravitazionale(GOCE));

	cout << endl << "Accelerazione di gravità sul satellite GOCE:" << endl;
	cout << endl << setprecision(14) << "g = " << -g.Modulo() << " m/s^2" << endl;


	//Aggiungo la catena montuosa

	cout << endl << "Aggiungo la catena montuosa..." << endl;

	CampoVettoriale deltag(GOCE); //Utilizzato per calcolare la variazione di g

	double rho_m = 3000.;		     //Densità roccia
	double r_m = 1000.;		     //Raggio montagne (schematizzate come sfere)
	double V_m = pow(r_m,3);		     //Volume montagne (schematizzate come sfere)
	double m_m = (rho_m*(4/3*M_PI*V_m)); //Massa montagne

	PuntoMateriale** montagna = new PuntoMateriale*[100];
	for(int i=0; i<100; i++){
		double phi = i*2.*r_m/RAGGIOTERRA;
		montagna[i] = new PuntoMateriale(m_m,0.,RAGGIOTERRA*cos(phi),RAGGIOTERRA*sin(phi),0.);
		g.Somma(montagna[i]->CampoGravitazionale(GOCE));
		deltag.Somma(montagna[i]->CampoGravitazionale(GOCE)); //Calcolo la variazione di g dovuta all'aggiunta della catena montuosa
	}

	cout << endl << "Accelerazione di gravità sul satellite GOCE:" << endl;
	cout << endl << setprecision(14) << "g = " << -g.Modulo() << " m/s^2" << endl;

	double variazione = (deltag.Modulo()/g.Modulo());

	if(variazione < 1.E-13)
		cout << endl << "La variazione di g dovuta alla catena montuosa è troppo piccola per essere misurata da GOCE!" << endl;

	else
		cout << endl << "Variazione di g rilevata: " << scientific << setprecision(6) << variazione << endl;

	cout << endl;

	delete montagna;

	return 0;
}


