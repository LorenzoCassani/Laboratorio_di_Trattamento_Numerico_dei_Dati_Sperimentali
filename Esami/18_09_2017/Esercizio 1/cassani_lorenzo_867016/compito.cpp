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
	double tmax=1.;   //tempo finale
	double x0=0.;	  //posizione iniziale
	double v0=0.;	  //velocità iniziale

	double eta=0.83;
	double rho=2700;
	double rho0=1250;
	double g=9.81;
	double R=0.02;	

	int nstep = int(tmax/h+0.5); //Calcolo il numero di passi

	VettoreLineare x(2);
	x.SetComponent(0,x0); //Setta la posizione iniziale in posizione 0
	x.SetComponent(1,v0); //Setta la velocità iniziale in posizione 1

	RungeKutta myRungeKutta; //Creo il solutore di equazioni differenziali (metodo di Runge-Kutta)
	Sfera *sfe = new Sfera(eta,rho,rho0,R); //Creo la sfera

	TGraph *gr1 = new TGraph(); //Grafico posizione

	TPad *p1 = new TPad();
	TPad *p2 = new TPad();


	//Domanda 1:

	for(int step=0; step<nstep; step++){
		gr1->SetPoint(step,t,x.GetComponent(0));
		x = myRungeKutta.Passo(t,x,h,sfe);
		t = t+h;
	}

	double pos1 = x.GetComponent(0);

	cout << endl << "Posizione x all'istante t=1 s = " << pos1 << " m" << endl;


	//Domanda 2:

	h = h/2.; //Dimezzo il passo
	nstep = int(tmax/h+0.5); //Calcolo il numero di passi

	t=0.;
	x.SetComponent(0,x0);
	x.SetComponent(1,v0);

	for(int step=0; step<nstep; step++){
		x = myRungeKutta.Passo(t,x,h,sfe);
		t = t+h;
	}

	double pos2 = x.GetComponent(0);

	h=atof(argv[1]); //Riporto il passo a quello iniziale

	double k = (16*abs(pos1-pos2))/(15*pow(h,4));

	cout << endl << "Errore sulla posizione x all'istante t=1 s = " << k*pow(h,4) << " m" << endl;

	


	//Domanda 3:

	double sigma1=1;  //Precisione 1 sulla misura di rho
	double sigma2=2;  //Precisione 2 sulla misura di rho
	double sigma3=5;  //Precisione 3 sulla misura di rho
	double sigma4=7;  //Precisione 4 sulla misura di rho
	double sigma5=10; //Precisione 5 sulla misura di rho
	double sigma6=20; //Precisione 6 sulla misura di rho

	Random generatore(1); //Creo il generatore di numeri casuali

	generatore.SetA(1664525);    //Valore parametro a
	generatore.SetC(1013904223); //Valore parametro c
	generatore.SetM(pow(2,31));  //Valore parametro m

	TH1F *h1 = new TH1F("#sigma_{#rho}=1" ,"Distribuzione dei valori della posizione dopo 1 s",100,pos1-3*sigma1,pos1+3*sigma1);
	TH1F *h2 = new TH1F("#sigma_{#rho}=2" ,"Distribuzione dei valori della posizione dopo 1 s",100,pos1-3*sigma2,pos1+3*sigma2);
	TH1F *h3 = new TH1F("#sigma_{#rho}=5" ,"Distribuzione dei valori della posizione dopo 1 s",100,pos1-3*sigma3,pos1+3*sigma3);
	TH1F *h4 = new TH1F("#sigma_{#rho}=7" ,"Distribuzione dei valori della posizione dopo 1 s",100,pos1-3*sigma4,pos1+3*sigma4);
	TH1F *h5 = new TH1F("#sigma_{#rho}=10","Distribuzione dei valori della posizione dopo 1 s",100,pos1-3*sigma5,pos1+3*sigma5);
	TH1F *h6 = new TH1F("#sigma_{#rho}=20","Distribuzione dei valori della posizione dopo 1 s",100,pos1-3*sigma6,pos1+3*sigma6);

	for(int i=0; i<1000; i++){

		t=0.;
		x.SetComponent(0,x0);
		x.SetComponent(1,v0);
		sfe->SetRho(generatore.RandGauss(rho,sigma1));

		for(int step=0; step<nstep; step++){
			x = myRungeKutta.Passo(t,x,h,sfe);
			t = t+h;
		}

		h1->Fill(x.GetComponent(0));
	}

	for(int i=0; i<1000; i++){

		t=0.;
		x.SetComponent(0,x0);
		x.SetComponent(1,v0);
		sfe->SetRho(generatore.RandGauss(rho,sigma2));

		for(int step=0; step<nstep; step++){
			x = myRungeKutta.Passo(t,x,h,sfe);
			t = t+h;
		}

		h2->Fill(x.GetComponent(0));
	}

	for(int i=0; i<1000; i++){

		t=0.;
		x.SetComponent(0,x0);
		x.SetComponent(1,v0);
		sfe->SetRho(generatore.RandGauss(rho,sigma3));

		for(int step=0; step<nstep; step++){
			x = myRungeKutta.Passo(t,x,h,sfe);
			t = t+h;
		}

		h3->Fill(x.GetComponent(0));
	}

	for(int i=0; i<1000; i++){

		t=0.;
		x.SetComponent(0,x0);
		x.SetComponent(1,v0);
		sfe->SetRho(generatore.RandGauss(rho,sigma4));

		for(int step=0; step<nstep; step++){
			x = myRungeKutta.Passo(t,x,h,sfe);
			t = t+h;
		}

		h4->Fill(x.GetComponent(0));
	}

	for(int i=0; i<1000; i++){

		t=0.;
		x.SetComponent(0,x0);
		x.SetComponent(1,v0);
		sfe->SetRho(generatore.RandGauss(rho,sigma5));

		for(int step=0; step<nstep; step++){
			x = myRungeKutta.Passo(t,x,h,sfe);
			t = t+h;
		}

		h5->Fill(x.GetComponent(0));
	}

	for(int i=0; i<1000; i++){

		t=0.;
		x.SetComponent(0,x0);
		x.SetComponent(1,v0);
		sfe->SetRho(generatore.RandGauss(rho,sigma6));

		for(int step=0; step<nstep; step++){
			x = myRungeKutta.Passo(t,x,h,sfe);
			t = t+h;
		}

		h6->Fill(x.GetComponent(0));
	}


	TGraph *gr2 = new TGraph(); //Grafico dipendenza RMS

	gr2->SetPoint(0,sigma1,h1->GetRMS());
	gr2->SetPoint(1,sigma2,h2->GetRMS());
	gr2->SetPoint(2,sigma3,h3->GetRMS());
	gr2->SetPoint(3,sigma4,h4->GetRMS());
	gr2->SetPoint(4,sigma5,h5->GetRMS());
	gr2->SetPoint(5,sigma6,h6->GetRMS());


	//Domanda 4:

	cout << endl << "L'errore introdotto dall'incertezza 1 kg/m^3 è: "
	<< setprecision(1) << h1->GetRMS() << " m" << endl;

	cout << "Affinchè l'errore così introdotto superi quello di troncamento il passo dev'essere minore di: "
	<< setprecision(4) << pow((h1->GetRMS()/k),0.25) << " m" << endl;


	TCanvas *c1 = new TCanvas("c1","c1",1000,100,2000,700);
	c1->Divide(1,2);

	c1->cd(1);
	p1->SetGrid();
	p1->Draw();
	p1->cd();
	gr1->GetXaxis()->SetTitle("t [s]");
	gr1->GetYaxis()->SetTitle("x [m]");
	gr1->SetTitle("Moto di una sfera in un fluido viscoso (metodo di Runge-Kutta)");
	gr1->Draw("AL");

	c1->cd(2);
	p2->SetGrid();
	p2->Draw();
	p2->cd();
	gr2->GetXaxis()->SetTitle("#sigma_{#rho}");
	gr2->GetYaxis()->SetTitle("RMS");
	gr2->SetTitle("RMS(#sigma_{#rho}) (metodo di Runge-Kutta)");
	gr2->Draw("AL*");

	app.Run();

	return 0;
}
