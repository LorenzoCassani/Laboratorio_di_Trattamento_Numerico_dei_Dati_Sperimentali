#include "lib.h"

using namespace std;

int main(){

	Vettore *a = new Vettore(10);
	for(unsigned int i=0; i<a->GetN(); i++) a->SetComponent(i,10.*i);
	cout << endl << "Vettore a iniziale: " << endl;
	for(unsigned int i=0; i<a->GetN(); i++) cout << i << ") " << a->GetComponent(i) << endl;
	if(a->GetN()>0){
		Vettore b=(*a); //Vettore b copia di a
		cout << endl << "Vettore b iniziale: " << endl;
		for(unsigned int i=0; i<b.GetN(); i++) cout << i << ") " << b.GetComponent(i) << endl;
		//Modifichiamo il vettore b e verifichiamo che le modifiche siano corrette
		for(unsigned int i=0; i<b.GetN(); i++) b.SetComponent(i,9.*i);
		cout << endl << "Vettore b modificato: " << endl;
		for(unsigned int i=0; i<b.GetN(); i++) cout << i << ") " << b.GetComponent(i) << endl;
	}
	//Stampiamo di nuovo il vettore a: il risultato è quello che ci aspettiamo?
	cout << endl << "Vettore a: " << endl;
	for(unsigned int i=0; i<a->GetN(); i++) cout << i << ") " << a->GetComponent(i) << endl;
	cout << endl << "Delete a..." << endl;
	delete a;
	cout << endl << "Deleted a" << endl;

        return 0;
}


