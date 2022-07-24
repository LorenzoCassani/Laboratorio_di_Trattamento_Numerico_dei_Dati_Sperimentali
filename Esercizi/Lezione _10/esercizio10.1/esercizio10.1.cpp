#include "lib.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

using namespace std;

int main(int argc, char **argv){

	if(argc!=2){
		cerr << "Usage: " << argv[0] << " <passo_di_integrazione>" << endl;
		exit(-1);
	}
	double h=atof(argv[1]);

	TApplication app("app",0,0);

	double t=0.;	 //tempo iniziale
	double tmax=70.; //tempo finale
	double x0=0.;	 //posizione iniziale
	double v0=1.;	 //velocità iniziale
	double omega=1.; //Pulsazione

	int nstep = int(tmax/h+0.5); //Calcolo il numero di passi

	VettoreLineare x(2);
	x.SetComponent(0,x0); //Setta la posizione iniziale in posizione 0
	x.SetComponent(1,v0); //Setta la velocità iniziale in posizione 1

	Eulero myEuler; //Creo il solutore di equazioni differenziali (metodo di Eulero)
	OscillatoreArmonico *osc = new OscillatoreArmonico(omega); //Creo l'oscillatore armonico

	TGraph *gr1 = new TGraph(); //Grafico posizione
	TGraph *gr2 = new TGraph(); //Grafico errore

	TPad *p1 = new TPad();
	TPad *p2 = new TPad();

	for(int step=0; step<nstep; step++){
		//Metto il punto nel grafico
		gr1->SetPoint(step,t,x.GetComponent(0));
		gr2->SetPoint(step,t,x.GetComponent(0)-sin(t)); //x=sin(t) è la soluzione esatta
		//Calcolo il prossimo punto
		x = myEuler.Passo(t,x,h,osc);
		t = t+h;
	}

	TCanvas *c1 = new TCanvas("c1","c1",1000,100,1000,700);
	c1->Divide(1,2);

	c1->cd(1);
	p1->SetGrid();
	p1->Draw();
	p1->cd();
	gr1->GetXaxis()->SetTitle("t [s]");
	gr1->GetYaxis()->SetTitle("x [m]");
	gr1->SetTitle("Oscillatore armonico (metodo di Eulero)");
	gr1->Draw("AL");

	c1->cd(2);
	p2->SetGrid();
	p2->Draw();
	p2->cd();
	gr2->GetXaxis()->SetTitle("t [s]");
	gr2->GetYaxis()->SetTitle("errore [m]");
	gr2->SetTitle("Errore (metodo di Eulero)");
	gr2->Draw("AL");

	app.Run();

	return 0;
}
