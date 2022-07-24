#include "lib.h"

using namespace std;

int main(){

	double a, b, v;

	for(int i=1; i<=20;i++){

		a=i*M_PI;
		b=i*M_PI+M_PI/2;

		cout << endl << "L'intervallo di ricerca è [" << a << "," << b << "]" << endl;

		FunzioneBase *f = new ENRA();

		Bisezione x(a,b,10E-6,f);

		double zero = x.CercaZeri(a,b);

		int cifre_significative = -log10(x.GetPrecisione());

		cout << fixed;
		cout << "Ho trovato uno zero in x = " << setprecision(cifre_significative) << zero << endl;

		delete f;
	}

        return 0;
}


