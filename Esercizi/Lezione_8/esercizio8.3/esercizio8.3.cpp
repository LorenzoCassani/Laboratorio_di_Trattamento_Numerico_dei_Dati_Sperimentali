#include "lib.h"

#include "TH1F.h"
#include "TApplication.h"
#include "TCanvas.h"

using namespace std;

int main(int argc, char **argv){

	if(argc!=3){
		cerr << "Usage: " << argv[0] << " <n> <ns>" << endl;
		return -1;
	}

	TApplication app("app",&argc,argv);

	int n = atoi(argv[1]);
	int ns = atoi(argv[2]);

	if(ns > n/2){
		cerr << endl << "ns deve essere minore o uguale a n/2!" << endl << endl;
		return -1;
	}

	double* v = new double [n];
	double* vs = new double [n-(n%ns)];

	Random rand(3);

	rand.SetA(1664525);    //Valore parametro a
	rand.SetC(1013904223); //Valore parametro c
	rand.SetM(pow(2,31));  //Valore parametro m

	TH1F h1("numeri generati","Distribuzione dei numeri generati",100,0,1);
	TH1F h2("somme","Distribuzione delle somme",100,0,ns);

	//Creazione della serie di numeri
	double r;
	for(int i=0; i<n; i++){
		r = rand.Rand01();
		h1.Fill(r);
		v[i] = r;
	}

	//Calcolo della media e varianza della serie di numeri
	cout << endl << "Media serie numeri: " << media(v,n) << endl;
	cout << "Varianza serie numeri: " << var(v,n) << endl << endl;

	//Creazione della serie delle somme
	for(int j=0, i=0; j<(n-(n%ns))/ns; j++){
		vs[j] = 0;
		for(unsigned int k=0; k<ns; k++, i++) vs[j] = vs[j]+v[i];
		h2.Fill(vs[j]);
	}

	//Calcolo della media e varianza della serie delle somme
	cout << endl << "Media serie somme: " << media(vs,(n-(n%ns))/ns) << endl;
	cout << "Varianza serie somme: " << var(vs,(n-(n%ns))/ns) << endl << endl;

	delete[] v;
	delete[] vs;

	TCanvas *c1 = new TCanvas("c1","c1",1200,600);

	c1->Divide(2,1);

	c1->cd(1);
	h1.Draw();

	c1->cd(2);
	h2.Draw();

	app.Run();

	return 0;
}
