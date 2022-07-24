#include "lib.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TGraph.h"

using namespace std;

int main(int argc, char **argv){

	TApplication app("app",&argc,argv);

	Seno *sin = new Seno; //Funzione integranda
	double max=1;	      //Massimo della funzione sin(x)

	Integral integrale(0,M_PI,sin); //Integrale da 0 a Pi di sin(x)

	int nin=100;	//n iniziale
	int nfin=10000; //n finale
	int step=100;	//Valore degli step

	int n = nfin/step;

	double x[n];	 //Array per i valori delle ascisse
	double media[n]; //Array per i valori dell'integrale calcolati col metodo della media
	double hom[n];   //Array per i valori dell'integrale calcolati col metodo Hit-or-miss
	
	
	int conta=0;
	for(int i=nin; i<=nfin; i+=step){
		x[conta]=i;
		media[conta] = integrale.Media(i);
		hom[conta] = integrale.HitOrMiss(i,max);
		conta++;
	}


	TCanvas *c1 = new TCanvas("c1","c1",1000,100,1000,700);
	c1->Divide(1,2);

	TGraph *gr1 = new TGraph(n,x,media); //Grafico metodo della media
	TGraph *gr2 = new TGraph(n,x,hom);   //Grafico metodo Hit-or-miss

	TPad *p1 = new TPad(); //Reticolo metodo della Media
	TPad *p2 = new TPad(); //Reticolo metodo Hit-or-miss

	c1->cd(1);
	p1->SetGrid();
	p1->Draw();
	p1->cd();
	gr1->SetTitle("Andamento del valore dell'integrale (metodo della media)");
	gr1->GetXaxis()->SetTitle("Numero di punti");
	gr1->GetYaxis()->SetTitle("Valore dell'integrale di sin(x) tra [0,#pi]");
	gr1->Draw("AP*");

	c1->cd(2);
	p2->SetGrid();
	p2->Draw();
	p2->cd();
	gr2->SetTitle("Andamento del valore dell'integrale (metodo Hit-or-miss)");
	gr2->GetXaxis()->SetTitle("Numero di punti");
	gr2->GetYaxis()->SetTitle("Valore dell'integrale di sin(x) tra [0,#pi]");
	gr2->Draw("AP*");

	app.Run();

	return 0;
}
