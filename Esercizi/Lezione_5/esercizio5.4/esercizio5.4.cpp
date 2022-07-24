#include "lib.h"

using namespace std;

int main(){

	double a, b, precision;

	cout << endl << "Inserisci i due estremi dell'intervallo:" << endl;
	cin >> a;
	cin >> b;

	if(a>b) swap(a,b);

	cout << endl << "L'intervallo di ricerca è [" << a << "," << b << "]" << endl;
	cout << endl << "Inserisci la precisione:" << endl;
	cin >> precision;

	FunzioneBase *f = new Parabola(3.,5.,-2.);

	Bisezione x(a,b,precision,f);

	if(x.Trovato()==false){
		cout << endl << "Non ho trovato zeri nell'intervallo!" << endl << endl;
		exit(-1);
	}

	else{

		int cifre_significative = -log10(x.GetPrecisione());

		cout << fixed;
		cout << endl << "Ho trovato uno zero in x = " << setprecision(cifre_significative) << x.CercaZeri(a,b) << endl;
		cout << "Incertezza: " << x.Incertezza() << endl;
	}

	delete f;

	return 0;
}


