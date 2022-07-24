#include "lib.h"

using namespace std;

int main(int argc, char** argv){

	//Controllo degli argomenti e stampa di un aiuto se sbagliati
	if(argc<3){
		cerr << "Usage: " << argv[0] << " <N_points> <filename>" << endl;
		return -1;
	}
	
	//Costruzione della struttura Vettore
	unsigned int N = atoi(argv[1]);
	struct vettore v = read(N,argv[2]);
	
	//Stampo il vettore, prima e dopo il riordinamento
	cout << "Before sorting..." << endl;
	print(v,"Before.dat");
	selection_sort(v);
        cout << "After sorting..."<<endl;
	print(v,"After.dat");

	//Finito di lavorare, cancello la memoria occupata ed esco dal programma
	delete[] v.v;

	return 0;
}



