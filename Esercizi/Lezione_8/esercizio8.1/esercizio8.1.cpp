#include "lib.h"

#include "TH1F.h"
#include "TApplication.h"
#include "TCanvas.h"

using namespace std;

int main(int argc, char **argv){

	TApplication app("app",&argc,argv);

	double a;

	cout << endl << "Inserire la larghezza (in sigma) in cui visualizzare l'istogramma: ";
	cin >> a;

	Gaussiana g1(1,1); //Gaussiana: mu = sigma = 1
	Gaussiana g2(0,2); //Gaussiana: mu = 0, sigma = 2

	TH1F h1("mu = sigma = 1","Distribuzione Gaussiana (Metodo Accept-Reject)",100,1-a,1+a);
	TH1F h2("mu = 0, sigma = 2","Distribuzione Gaussiana (Metodo Accept-Reject)",100,-2*a,2*a);

	AcceptReject r1(&g1,1-a,1+ a,1);
	AcceptReject r2(&g2,-2*a,2*a,1);
	
	double n;
	for(int i=0; i<100000; i++){
		n = r1.RandGauss();
		h1.Fill(n);
		n = r2.RandGauss();
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
