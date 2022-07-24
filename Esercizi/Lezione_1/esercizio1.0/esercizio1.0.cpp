#include <iostream>

using namespace std;

void scambiaByValue(double,double); //NON FUNZIONA

void scambiaByRef(double&,double&);

void scambiaByPointer(double*,double*);

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

void scambiaByValue(double a, double b){ //NON FUNZIONA

	double appo=a;
	a=b;
	b=appo;
}

void scambiaByRef(double &a, double &b){

	double appo=a;
	a=b;
	b=appo;
}

void scambiaByPointer(double *a, double *b){

	double appo=*a;
	*a=*b;
	*b=appo;
}
