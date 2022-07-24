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

	double t=0.;	  //tempo iniziale
	double tmax=30.;  //tempo finale
	double l=1.;	  //Lunghezza del filo del pendolo
	double A=-M_PI/4; //Ampiezza angolare
	double omega0=0.; //Velocità angolare iniziale 

	int nstep = int(tmax/h+0.5); //Calcolo il numero di passi

	//(metodo di Eulero)
	VettoreLineare x1(2);
	x1.SetComponent(0,A);      //Setta l'ampiezza angolare in posizione 0
	x1.SetComponent(1,omega0); //Setta la velocità angolare iniziale in posizione 1

	//(metodo di Runge-Kutta)
	VettoreLineare x2(2);
	x2.SetComponent(0,A);      //Setta l'ampiezza angolare in posizione 0
	x2.SetComponent(1,omega0); //Setta la velocità angolare iniziale in posizione 1

	Pendolo *pen = new Pendolo(l); //Creo il pendolo
	Eulero myEulero;	       //Creo il solutore di equazioni differenziali (metodo di Eulero)
	RungeKutta myRungeKutta;       //Creo il solutore di equazioni differenziali (metodo di Runge-Kutta)

	TGraph *gr1 = new TGraph(); //Grafico posizione angolare in funzione del tempo (metodo di Eulero)
	TGraph *gr2 = new TGraph(); //Grafico posizione angolare in funzione del tempo (metodo di Runge-Kutta)
	TGraph *gr3 = new TGraph(); //Grafico spazio delle fasi (metodo di Eulero)
	TGraph *gr4 = new TGraph(); //Grafico spazio delle fasi (metodo di Runge-Kutta)
	TGraph *gr5 = new TGraph(); //Grafico periodo in funzione dell'ampiezza angolare (metodo di Eulero)
	TGraph *gr6 = new TGraph(); //Grafico periodo in funzione dell'ampiezza angolare (metodo di Runge-Kutta)

	TPad *p1 = new TPad();
	TPad *p2 = new TPad();
	TPad *p3 = new TPad();
	TPad *p4 = new TPad();
	TPad *p5 = new TPad();
	TPad *p6 = new TPad();

	for(int step=0; step<nstep; step++){
		gr1->SetPoint(step,t,x1.GetComponent(0)); 
		gr2->SetPoint(step,t,x2.GetComponent(0));
		gr3->SetPoint(step,x1.GetComponent(0),x1.GetComponent(1)); 
		gr4->SetPoint(step,x2.GetComponent(0),x2.GetComponent(1));

		x1 = myEulero.Passo(t,x1,h,pen);
		x2 = myRungeKutta.Passo(t,x2,h,pen);
		t = t+h;
	}


	//Calcolo del periodo:

	double T; //Periodo
	double v; //velocità

	int N=30;	 //Numero di volte per cui calcolare il periodo
	double step=A/N; //Calcolo lo step


	//(metodo di Eulero)

	for(int i=0; i<N; i++){

		//Condizioni iniziali
		t=0.;
		x1.SetComponent(0,A);
		x1.SetComponent(1,0.);

		while(x1.GetComponent(1)>=0){ //Si procede finchè la velocità non inverte di segno
			v  = x1.GetComponent(1);
			x1 = myEulero.Passo(t,x1,h,pen);
			t  = t+h;
		}

		T = (t-h)-v*h/(x1.GetComponent(1)-v); //Interpolazione
		T = 2*T; //Il periodo è due volte il semi-periodo appena calcolato

		gr5->SetPoint(i,-A,T);
		A-=step;
	}

	A=-M_PI/4; //Riporto l'ampiezza angolare al suo valore iniziale


	//(metodo di Runge-Kutta)

	for(int i=0; i<N; i++){

		//Condizioni iniziali
		t=0.;
		x2.SetComponent(0,A);
		x2.SetComponent(1,0.);

		while(x2.GetComponent(1)>=0){ //Si procede finchè la velocità non inverte di segno
			v  = x2.GetComponent(1);
			x2 = myRungeKutta.Passo(t,x2,h,pen);
			t  = t+h;
		}

		T = (t-h)-v*h/(x2.GetComponent(1)-v); //Interpolazione
		T = 2*T; //Il periodo è due volte il semi-periodo appena calcolato

		gr6->SetPoint(i,-A,T);
		A-=step;
	}

	A=-M_PI/4; //Riporto l'ampiezza angolare al suo valore iniziale


	TCanvas *c1 = new TCanvas("c1","c1",1000,100,2000,700);
	c1->Divide(3,2);

	c1->cd(1);
	p1->SetGrid();
	p1->Draw();
	p1->cd();
	gr1->GetXaxis()->SetTitle("t [s]");
	gr1->GetYaxis()->SetTitle("#theta [rad]");
	gr1->SetTitle("#theta(t) (metodo di Eulero)");
	gr1->Draw("AL");

	c1->cd(2);
	p3->SetGrid();
	p3->Draw();
	p3->cd();
	gr3->GetXaxis()->SetTitle("#theta [rad]");
	gr3->GetYaxis()->SetTitle("#omega [rad/s]");
	gr3->SetTitle("Spazio delle fasi (metodo di Eulero)");
	gr3->Draw("AL");

	c1->cd(3);
	p5->SetGrid();
	p5->Draw();
	p5->cd();
	gr5->GetXaxis()->SetTitle("A [rad]");
	gr5->GetYaxis()->SetTitle("T [s]");
	gr5->SetTitle("Periodo in funzione dell'ampiezza (metodo di Eulero)");
	gr5->Draw("AL*");

	c1->cd(4);
	p2->SetGrid();
	p2->Draw();
	p2->cd();
	gr2->GetXaxis()->SetTitle("t [s]");
	gr2->GetYaxis()->SetTitle("#theta [rad]");
	gr2->SetTitle("#theta(t) (metodo di Runge-Kutta)");
	gr2->Draw("AL");

	c1->cd(5);
	p4->SetGrid();
	p4->Draw();
	p4->cd();
	gr4->GetXaxis()->SetTitle("#theta [rad]");
	gr4->GetYaxis()->SetTitle("#omega [rad/s]");
	gr4->SetTitle("Spazio delle fasi (metodo di Runge-Kutta)");
	gr4->Draw("AL");

	c1->cd(6);
	p6->SetGrid();
	p6->Draw();
	p6->cd();
	gr6->GetXaxis()->SetTitle("A [rad]");
	gr6->GetYaxis()->SetTitle("T [s]");
	gr6->SetTitle("Periodo in funzione dell'ampiezza (metodo di Runge-Kutta)");
	gr6->Draw("AL*");

	app.Run();

	return 0;
}
