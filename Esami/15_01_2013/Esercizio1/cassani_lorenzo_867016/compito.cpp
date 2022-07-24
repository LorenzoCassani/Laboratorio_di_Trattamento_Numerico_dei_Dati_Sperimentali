#include "lib.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"
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

	double t=0.;	  //tempo iniziale
	double tmax=56.;  //tempo finale
	double x0=1.;	  //posizione iniziale
	double v0=0.;	  //velocità iniziale
	double omega=1.2; //Pulsazione

	int nstep = int(tmax/h+0.5); //Calcolo il numero di passi

	VettoreLineare x(2);
	x.SetComponent(0,x0); //Setta la posizione iniziale in posizione 0
	x.SetComponent(1,v0); //Setta la velocità iniziale in posizione 1

	RungeKutta myRungeKutta; //Creo il solutore di equazioni differenziali (metodo di Runge-Kutta)
	OscillatoreArmonico *osc = new OscillatoreArmonico(omega); //Creo l'oscillatore armonico

	TPad *p1 = new TPad();

	//Domanda 1:

	for(int step=0; step<nstep; step++){
		x = myRungeKutta.Passo(t,x,h,osc);
		t = t+h;
	}

	double pos = x.GetComponent(0);

	cout << endl << "Posizione x all'istante t=56 s = " << pos << " m" << endl;

	//Domanda 2:

	double sigma1=0.005; //Precisione 1 sulla misura di x0

	Random generatore(1); //Creo il generatore di numeri casuali

	generatore.SetA(1664525);    //Valore parametro a
	generatore.SetC(1013904223); //Valore parametro c
	generatore.SetM(pow(2,31));  //Valore parametro m

	TH1F *h1 = new TH1F("#sigma_{x_{0}}","Distribuzione dei valori della posizione dopo 56 s",100,pos-3*sigma1,pos+3*sigma1);

	for(int i=0; i<1000; i++){

		t=0.;
		x.SetComponent(0,generatore.RandGauss(x0,sigma1));
		x.SetComponent(1,v0);

		for(int step=0; step<nstep; step++){
			x = myRungeKutta.Passo(t,x,h,osc);
			t = t+h;
		}

		h1->Fill(x.GetComponent(0));
	}

	//Domanda 3:

	double sigma2=0.010; //Precisione 2 sulla misura di x0
	double sigma3=0.015; //Precisione 3 sulla misura di x0
	double sigma4=0.020; //Precisione 4 sulla misura di x0

	TH1F *h2 = new TH1F("#sigma_{x_{0}}=0.010","Distribuzione dei valori della posizione dopo 56 s",100,pos-3*sigma2,pos+3*sigma2);
	TH1F *h3 = new TH1F("#sigma_{x_{0}}=0.015","Distribuzione dei valori della posizione dopo 56 s",100,pos-3*sigma3,pos+3*sigma3);
	TH1F *h4 = new TH1F("#sigma_{x_{0}}=0.020","Distribuzione dei valori della posizione dopo 56 s",100,pos-3*sigma4,pos+3*sigma4);

	TGraph *gr1 = new TGraph();

	for(int i=0; i<1000; i++){

		t=0.;
		x.SetComponent(0,generatore.RandGauss(x0,sigma2));
		x.SetComponent(1,v0);

		for(int step=0; step<nstep; step++){
			x = myRungeKutta.Passo(t,x,h,osc);
			t = t+h;
		}

		h2->Fill(x.GetComponent(0));
	}

	for(int i=0; i<1000; i++){

		t=0.;
		x.SetComponent(0,generatore.RandGauss(x0,sigma3));
		x.SetComponent(1,v0);

		for(int step=0; step<nstep; step++){
			x = myRungeKutta.Passo(t,x,h,osc);
			t = t+h;
		}

		h3->Fill(x.GetComponent(0));
	}

	for(int i=0; i<1000; i++){

		t=0.;
		x.SetComponent(0,generatore.RandGauss(x0,sigma4));
		x.SetComponent(1,v0);

		for(int step=0; step<nstep; step++){
			x = myRungeKutta.Passo(t,x,h,osc);
			t = t+h;
		}

		h4->Fill(x.GetComponent(0));
	}

	gr1->SetPoint(0,sigma1,h1->GetRMS());
	gr1->SetPoint(1,sigma2,h2->GetRMS());
	gr1->SetPoint(2,sigma3,h3->GetRMS());
	gr1->SetPoint(3,sigma4,h4->GetRMS());

	TCanvas *c1 = new TCanvas("c1","c1",1000,100,2000,700);
	c1->Divide(1,2);

	c1->cd(1);
	h1->Draw();

	c1->cd(2);
	p1->SetGrid();
	p1->Draw();
	p1->cd();
	gr1->GetXaxis()->SetTitle("t [s]");
	gr1->GetYaxis()->SetTitle("errore [m]");
	gr1->SetTitle("Errore (metodo di Runge-Kutta)");
	gr1->Draw("AL*");

	app.Run();

	return 0;
}
