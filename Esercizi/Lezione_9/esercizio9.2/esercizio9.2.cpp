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

	double x[n];	    //Array per i valori delle ascisse
	double errmedia[n]; //Array per gli errori sulla misura dell'integrale calcolati col metodo della media
	double errhom[n];   //Array per gli errori sulla misura dell'integrale calcolati col metodo Hit-or-miss
	
	
	int conta=0;
	for(int i=nin; i<=nfin; i+=step){
		integrale.Media(i);
		integrale.HitOrMiss(i,max);
		x[conta]=i;
		errmedia[conta] = integrale.GetErrorMedia();
		errhom[conta] = integrale.GetErrorHitOrMiss();
		conta++;
	}


	TCanvas *c1 = new TCanvas("c1","c1",1000,100,1000,700);
	c1->Divide(1,2);

	TGraph *gr1 = new TGraph(n,x,errmedia); //Grafico metodo della media
	TGraph *gr2 = new TGraph(n,x,errhom);   //Grafico metodo Hit-or-miss

	TPad *p1 = new TPad(); //Reticolo metodo della Media
	TPad *p2 = new TPad(); //Reticolo metodo Hit-or-miss

	c1->cd(1);
	p1->SetGrid();
	p1->Draw();
	p1->cd();
	gr1->SetTitle("Andamento dell'errore sul calcolo dell'integrale (metodo della media)");
	gr1->GetXaxis()->SetTitle("Numero di punti");
	gr1->GetYaxis()->SetTitle("Errore sul calcolo dell'integrale di sin(x) tra [0,#pi]");
	gr1->Draw("AP*");

	c1->cd(2);
	p2->SetGrid();
	p2->Draw();
	p2->cd();
	gr2->SetTitle("Andamento dell'errore sul calcolo dell'integrale (metodo Hit-or-miss)");
	gr2->GetXaxis()->SetTitle("Numero di punti");
	gr2->GetYaxis()->SetTitle("Errore sul calcolo dell'integrale di sin(x) tra [0,#pi]");
	gr2->Draw("AP*");


	//Calcolo alpha

	integrale.Media(nin);
	integrale.HitOrMiss(nin,max);

	double errmediain = integrale.GetErrorMedia();
	double errhomin   = integrale.GetErrorHitOrMiss();

	integrale.Media(nfin);
	integrale.HitOrMiss(nfin,max);

	double errmediafin = integrale.GetErrorMedia();
	double errhomfin   = integrale.GetErrorHitOrMiss();

	double alphamedia = log(errmediain/errmediafin)/log((double)(nin)/nfin);
	double alphahom   = log(errhomin/errhomfin)/log((double)(nin)/nfin);

	cout << endl << "Metodo della media: l'andamento dell'errore è del tipo K*N^α con α = " << setprecision(1) << alphamedia << endl;
	cout << "Metodo Hit-or-miss: l'andamento dell'errore è del tipo K*N^α con α = " << setprecision(1) << alphahom << endl;

	cout << endl << "Numero di punti necessari per ottenere una precisione di 0.001: " << endl;
	cout << "Metodo della Media: " << scientific << integrale.NMedia(0.001) << endl;
	cout << "Metodo Hit-or-miss: " << scientific << integrale.NHitOrMiss(0.001) << endl;

	app.Run();

	return 0;
}
