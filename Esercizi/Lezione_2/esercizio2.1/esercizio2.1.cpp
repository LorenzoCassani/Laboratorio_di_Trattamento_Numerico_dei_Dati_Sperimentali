#include "lib.h"

using namespace std;

int main(){

	double a, b; 

	cout << "Inserire due numeri: " << endl;
	cin >> a >> b;

	for (int i=0; i<10; i++){
		scambiaByValue(a,b);
	}

	for (int i=0; i<10; i++){
		scambiaByRef(a,b);
	}

	for (int i=0; i<10; i++){
		scambiaByPointer(&a,&b);
	}

	cout << endl << "Numero complessivo di scambi effettuati: " << numeroScambi() << endl;

	return 0;
}



