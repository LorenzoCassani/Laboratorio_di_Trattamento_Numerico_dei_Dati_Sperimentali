#include "lib.h"

#define MASSATERRA 5.972E24
#define RAGGIOTERRA 6.371E6
#define MASSASOLE 1.989E30
#define RAGGIOSOLE 6.957E8
#define MASSALUNA 7.348E22
#define RAGGIOLUNA 1.737E6

using namespace std;

int main(){

	CorpoCeleste *Terra = new CorpoCeleste("Terra",MASSATERRA,RAGGIOTERRA);
	CorpoCeleste *Sole = new CorpoCeleste("Sole",MASSASOLE,RAGGIOSOLE);
	CorpoCeleste *Luna = new CorpoCeleste("Luna",MASSALUNA,RAGGIOLUNA);

	double distanza_Sole_Luna = 1.496E11;
	double distanza_Terra_Luna = 3.844E8;
	
	//Assegno le posizione ai rispettivi corpi celesti
	Terra->SetPosizione(0,0,0); 		//Origine del nostro s.d.r. cartesiano
	Sole->SetPosizione(distanza_Sole_Luna,0,0);	//Posizione relativa del Sole
	Luna->SetPosizione(0,distanza_Terra_Luna,0);	//Posizione relativa della Luna

	cout << "Potenziale gravitazionale generato dal Sole sulla Luna: " 
	<< Sole->PotenzialeGravitazionale(Luna->GetPosizione()) << " J/kg " << endl;
	cout << "Potenziale gravitazionale generato dalla Terra sulla Luna: " 
	<< Terra->PotenzialeGravitazionale(Luna->GetPosizione()) << " J/kg " << endl;

	delete Terra,Sole,Luna;

	return 0;
}


