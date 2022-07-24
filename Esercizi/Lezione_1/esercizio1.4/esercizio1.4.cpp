#include "lib.h"

using namespace std;

int main(){

	ifstream in;
	
	int scelta;

	cout << endl << "Specificare che tipo di operazione si intende compiere [1=num/2=vettR2/3=vettRn]" << endl;
	cin >> scelta;

	switch(scelta){

		case 1:
			
			double a;
			cout << endl << "Inserire il numero di cui si vuole calcolare il modulo: ";
			cin >> a;
			cout << endl << "Modulo: " << modulo(a) << endl;
			cout << endl;

			break;

		case 2:

			double v1,v2;
			cout << endl << "Inserire le due componenti del vettore di cui si vuole calcolare il modulo: " << endl;
			cin >> v1 >> v2;
			cout << endl << "Modulo: " << modulo(v1,v2) << endl;
			cout << endl;

			break;

		case 3:

			int n;
			double *arr;

			in.open("dati.dat");

			if(in.fail()){
  				cerr << endl << "ERRORE: il file dati.dat non è presente nella cartella!" << endl;
  				return -1;
			}

			in >> n;
			cout << endl << "Numero componenti: " << n << endl;

			arr=new double[n];

			for(int i=0; i<n; i++)
				in >> arr[i];

			in.close();

			cout << endl << "Modulo: " << modulo(arr,n) << endl;
			cout << endl;

			delete []arr;

			break;

		default:
 
     			cout << endl << "La direttiva inserita non è valida: si prega di inserire una delle tre direttive" << endl;
			cout << endl;
	}

	return 0;
}




