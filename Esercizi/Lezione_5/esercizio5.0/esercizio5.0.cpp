#include "lib.h"

using namespace std;

int main(){

	//Elettrone e CorpoCeleste sono anche Particella,
	//quindi le posso assegnare a dei puntatori a Particella
	Particella *a = new Particella(1.,2.);
	Particella *b = new Elettrone();
	Particella *c = new CorpoCeleste("Terra",5.9742E24,6.372797E6);

	a->Print(); //Metodo Print di Particella
	b->Print(); //Metodo Print di Elettrone
	c->Print(); //Metodo Print di CorpoCeleste

	delete a;
	delete b;
	delete c;

	return 0;
}


