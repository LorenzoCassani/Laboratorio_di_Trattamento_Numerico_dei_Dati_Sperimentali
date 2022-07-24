#include "lib.h"

using namespace std;

int main(){

	double a,b;

	cout << endl << "Il programma scambia il contenuto di due variabili double a e b" << endl;
	cout << endl << "Inserire il contenuto della variabile a: ";
	cin >> a;
	cout << "Inserire il contenuto della variabile b: ";
	cin >> b;

	//scambiaByValue(a,b); //NON FUNZIONA
	scambiaByRef(a,b);
	//scambiaByPointer(&a,&b);

	cout << endl << "Scambio..." << endl;
	cout << endl << "a = " << a << endl;
	cout << "b = " << b << endl;

	return 0;
}
