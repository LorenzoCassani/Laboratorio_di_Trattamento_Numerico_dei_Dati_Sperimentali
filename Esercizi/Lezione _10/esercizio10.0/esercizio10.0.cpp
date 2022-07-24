#include "lib.h"

using namespace std;

int main(){

	double x1,y1,z1; //Componenti primo vettore
	double x2,y2,z2; //Componenti secondo vettore

	cout << endl << "Inserire le componenti del primo vettore (v):" << endl;
	cin >> x1 >> y1 >> z1;

	cout << endl << "inserire le componenti del secondo vettore (w)" << endl;
	cin >> x2 >> y2 >> z2;

	VettoreLineare v(3);
	v.SetComponent(0,x1);
	v.SetComponent(1,y1);
	v.SetComponent(2,z1);

	VettoreLineare w(3);
	w.SetComponent(0,x2);
	w.SetComponent(1,y2);
	w.SetComponent(2,z2);

	cout << endl << "v = ("
	<< v.GetComponent(0) << ","
	<< v.GetComponent(1) << ","
	<< v.GetComponent(2) << ")" << endl;

	cout << endl << "w = ("
	<< w.GetComponent(0) << ","
	<< w.GetComponent(1) << ","
	<< w.GetComponent(2) << ")" << endl << endl;


	cout << endl << endl << "r = v + w" << endl;

	VettoreLineare r(3);
	r = v + w;
	cout << endl << "r = ("
	<< r.GetComponent(0) << ","
	<< r.GetComponent(1) << ","
	<< r.GetComponent(2) << ")" << endl;


	VettoreLineare versr(3);
	versr = r.Versore();
	cout << endl << "Versore di r = ("
	<< versr.GetComponent(0) << ","
	<< versr.GetComponent(1) << ","
	<< versr.GetComponent(2) << ")" << endl << endl;

	
	cout << endl << "cos(r,v) = " << r.Cos(v) << endl;
	cout << "cos(r,w) = " << r.Cos(w) << endl << endl;

	return 0;
}
