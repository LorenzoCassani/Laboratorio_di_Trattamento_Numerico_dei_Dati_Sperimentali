#include "lib.h"

using namespace std;

int main(int argc, char** argv){

	if(argc<3){
		cerr << "Usage: " << argv[0] << " <N_points> <filename>" << endl;
		return -1;
	}

	unsigned int N = atoi(argv[1]);
	VettoreDati v(N,argv[2]);

	cout << endl << "Statistica vettore caricato dal file: " << argv[2] << endl;

	cout << endl << "Media = " 		 << v.Media() << endl; 
	cout << 	"Varianza = " 		 << v.Var() << endl;
	cout << 	"Deviazione standard = " << v.StdDev() << endl;
	cout << 	"Mediana = " 		 << v.Mediana() << endl << endl;

	return 0;
}



