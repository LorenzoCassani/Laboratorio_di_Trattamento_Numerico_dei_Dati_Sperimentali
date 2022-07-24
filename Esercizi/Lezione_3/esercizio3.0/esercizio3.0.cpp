#include "lib.h"

using namespace std;

//Questo programma vuole come input da riga di comando le coordinate x, y, z di un punto e ne restituisce le coordinate sferiche e cilindriche

int main(int argc, char** argv){

	//Controlla gli argomenti
	if(argc!=4){
		cerr << "Usage: " << argv[0] << " <x> <y> <z>" << endl;
		return -1;
	}

	//Decodifica gli argomenti numerici (La funzione atof converte da char* a double)
	double x = atof(argv[1]);
	double y = atof(argv[2]);
	double z = atof(argv[3]);

	//Crea un oggetto posizione ed accede ai vari metodi
	Posizione P(x,y,z);
	cout << "Coordinate cartesiane: "
	     << P.GetX() << ", " << P.GetY() << ", " << P.GetZ() << endl;
	cout << "Coordinate cilindriche: "
	     << P.Rho() << ", " << P.Phi() << ", " << P.GetZ() << endl;
	cout << "Coordinate sferiche: "
	     << P.R() << ", " << P.Phi() << ", " << P.Theta() << endl;

	return 0;
}


