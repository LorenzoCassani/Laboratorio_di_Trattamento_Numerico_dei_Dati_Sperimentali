#include "lib.h"

using namespace std;

int main(){

	unsigned int N;
	unsigned int conta=0;
	double *array;

	ifstream in;
	ofstream out;

	in.open("dati.dat");

	if(in.fail()){
  		cerr << endl << "ERRORE: il file dati.dat non è presente nella cartella!" << endl;
  		return -1;
	}

	in >> N;
	cout << endl << "Numero dati: " << N << endl;

	array=new double[N];

	for(int i=0; i<N; i++)
		in >> array[i];

	in.close();

	//Scambio ciascun elemento di indice pari con il successivo (che è di indice dispari)
	while(conta<N){
		scambiaByRef(array[conta],array[conta+1]);
		conta=conta+2;
	}

	out.open("output.txt");

	for(int i=0; i<N; i++)
		out << array[i] << endl;

	out.close();

	cout << endl << "Programma eseguito: è possibile verificarne il risultato sul file output.txt" << endl;
	cout << endl;

	delete []array;

	return 0;
}
