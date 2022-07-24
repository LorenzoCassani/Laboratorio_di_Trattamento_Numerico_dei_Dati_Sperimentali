#include "lib.h"

using namespace std;

int main(){

	double x,y1,y2;
	double data[4] = {4, 7, 13, 16};
	double data_hugenumber[4] = {4+1E9, 7+1E9, 13+1E9, 16+1E9};
	int dim=4;

	cout << endl;

	cout << "0.3 - 0.2 - 0.1 = " << scientific << 0.3 - 0.2 - 0.1 << endl;
	cout << "0.3 - (0.2 + 0.1) = " << scientific << 0.3 - (0.2 + 0.1) << endl;
	cout << "0.4 - 0.3 - 0.1 = " << scientific << 0.4 - 0.3 - 0.1 << endl;
	cout << "0.4 - (0.3 + 0.1) = " << scientific << 0.4 - (0.3 + 0.1) << endl;
	cout << "(1E15 + 1) - 1E15 = " << scientific << (1E15 + 1) - 1E15 << endl;
	cout << "(1E16 + 1) - 1E16 = " << scientific << (1E16 + 1) - 1E16 << endl;

	cout << endl;

	cout << "---------------------------------------" << endl;
	cout << setw(5) << "x" << setw(17) << "fconerrore(x)" << setw(17) << "fmigliore(x)" << endl;
	cout << "---------------------------------------" << endl;

	for(int i=0; i<19; i++){

        	x = pow(10,i);
        	y1 = sqrt(x+1) - sqrt(x);
		y2 = 1/(sqrt(x+1) + sqrt(x)); //RAZIONALIZZO: aggiro il problema della sottrazione per ottenere una precisione maggiore

		cout << scientific << setprecision(0) << setw(5) << x << setprecision(6) << setw(17) << y1 << setw(17) << y2 << endl;
        }

	cout << endl;

	cout << "Media data: " << fixed << setprecision(0) << media(data,dim) << endl;
	cout << "Varianza data: " << fixed << setprecision(1) << var(data,dim) << endl;
	cout << "Media data_hugenumber: " << scientific << setprecision(6) << media(data_hugenumber,dim) << endl;
	cout << "Varianza data_hugenumber: " << scientific << setprecision(6) << var(data_hugenumber,dim) << endl;

	return 0;
}
