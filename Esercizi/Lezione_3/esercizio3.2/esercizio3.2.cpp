#include "lib.h"

using namespace std;

//Questo programma vuole come input da riga di comando le coordinate x, y, z di un punto e ne restituisce le coordinate sferiche e cilindriche

int main(){

	Particella *a = new Particella(1.,1.6E-19);
	Elettrone *e = new Elettrone();
	CorpoCeleste *c = new CorpoCeleste("Terra",5.9742E24,6.372797E6);

	//Metodi della classe base
	cout << a->GetMassa() << "," << a->GetCarica() << endl;
	//Metodi della classe derivata
	a->Print();
	//Metodi della classe base
	cout << e->GetMassa() << "," << e->GetCarica() << endl;
	//Metodi della classe derivata
	e->Print();
	//Metodi della classe base
	cout << c->GetMassa() << "," << c->GetCarica() << endl;
	//Metodi della classe derivata
	cout << c->GetNome() << endl;
	c->Print();

	Particella b(*a); //Posso costruire una Particella da una Particella
	Particella d(*e); //o da un Elettrone (che è una Particella),
	//Elettrone f(d);   //ma non il contrario!

	//CorpoCeleste g; //Non funziona
	CorpoCeleste g("Test",0.,0.); //Questo invece funziona
	g=(*c);
	g.Print();

	return 0;
}


