#include "lib.h"

#include "TH1F.h"
#include "TApplication.h"
#include "TCanvas.h"

using namespace std;

int main(int argc, char **argv){

	TApplication app("app",&argc,argv);

	Random rand(3); //Istanzio un (mio) oggetto Random

	rand.SetA(1664525);    //Valore parametro a
	rand.SetC(1013904223); //Valore parametro c
	rand.SetM(pow(2,31));  //Valore parametro m

	double a;

	cout << endl << "Inserire la larghezza (in sigma) in cui visualizzare l'istogramma: ";
	cin >> a;

	TH1F h1("rate = 0.1","Distribuzione esponenziale (Metodo della trasformata)",250,0,0.1*a);
	TH1F h2("mu = 0, sigma = 1","Distribuzione Gaussiana (Metodo Box-Muller)",100,-a,a);
	
	double n;
	for(int i=0; i<100000; i++){
		n = rand.RandEsp(0.1);
		h1.Fill(n);
		n = rand.RandGauss(0,1);
		h2.Fill(n);
	}

	TCanvas *c1 = new TCanvas("c1","c1");
	TCanvas *c2 = new TCanvas("c2","c2");

	c1->cd();
	h1.Draw();
	c2->cd();
	h2.Draw();

	app.Run();

	return 0;
}
