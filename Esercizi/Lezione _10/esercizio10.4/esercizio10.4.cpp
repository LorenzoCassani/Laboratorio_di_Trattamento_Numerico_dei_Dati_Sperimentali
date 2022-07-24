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

	double t=0.;	    //tempo iniziale
	double tmax=25.;    //tempo finale
	double x0=0.;	    //posizione iniziale
	double v0=0.;	    //velocità iniziale
	double omega0=10.;  //Pulsazione iniziale
	double omega=15.;   //Pulsazione
	double alpha=1/30.; //alpha

	int nstep = int(tmax/h+0.5); //Calcolo il numero di passi

	VettoreLineare x(2);
	x.SetComponent(0,x0); //Setta la posizione iniziale in posizione 0
	x.SetComponent(1,v0); //Setta la velocità iniziale in posizione 1

	RungeKutta myRungeKutta; //Creo il solutore di equazioni differenziali (metodo di Runge-Kutta)
	OscillatoreForzato *osc = new OscillatoreForzato(omega0,omega,alpha); //Creo l'oscillatore forzato

	TGraph *gr1 = new TGraph(); //Grafico posizione in funzione del tempo (metodo di Runge-Kutta)
	TGraph *gr2 = new TGraph(); //Grafico curva di risonanza (metodo di Runge-Kutta)

	TPad *p1 = new TPad();
	TPad *p2 = new TPad();

	for(int step=0; step<nstep; step++){
		//Metto il punto nel grafico
		gr1->SetPoint(step,t,x.GetComponent(0));
		//Calcolo il prossimo punto
		x = myRungeKutta.Passo(t,x,h,osc);
		t = t+h;
	}

	//Curva di risonanza:

	double A; //Ampiezza
	omega=omega0-1;	      //Pongo omega = omega0-1
	osc->SetOmega(omega); //Setto il nuovo valore di omega nell'oscillatore forzato

	int conta=0;	       //Contatore per i punti del grafico
	while(omega<omega0+1){ //Traccio la curva di risonanza per valori di omega da omega0-1 ad omega0+1

		//Condizioni iniziali
		t=0;
		x.SetComponent(0,x0);
		x.SetComponent(1,v0);

		//Per raggiungere una condizione di stabilità, integro per un tempo pari ad almeno dieci volte 1/α
		while(t<=10*1/osc->GetAlpha()){
			x = myRungeKutta.Passo(t,x,h,osc);
			t = t+h;
		}

		//Prima del cambio di segno della velocità
		if(x.GetComponent(1)>0){
			while(x.GetComponent(1)>=0){
				x = myRungeKutta.Passo(t,x,h,osc);
				t = t+h;
			}
		}

		//Dopo il cambio di segno della velocità
		else{
			while(x.GetComponent(1)<=0){
				x = myRungeKutta.Passo(t,x,h,osc);
				t = t+h;
			}
		}

		A = abs(x.GetComponent(0));
		gr2->SetPoint(conta,omega,A);
		omega+=0.05; //incremento omega di 0.05 alla volta
		osc->SetOmega(omega);
		conta++;
	}

	TCanvas *c1 = new TCanvas("c1","c1",1000,100,2000,700);
	c1->Divide(1,2);

	c1->cd(1);
	p1->SetGrid();
	p1->Draw();
	p1->cd();
	gr1->GetXaxis()->SetTitle("t [s]");
	gr1->GetYaxis()->SetTitle("x [m]");
	gr1->SetTitle("Oscillatore armonico smorzato con forzante (metodo di Runge-Kutta)");
	gr1->Draw("AL");

	c1->cd(2);
	p2->SetGrid();
	p2->Draw();
	p2->cd();
	gr2->GetXaxis()->SetTitle("#omega [rad/s]");
	gr2->GetYaxis()->SetTitle("A [rad]");
	gr2->SetTitle("Curva di risonanza (metodo di Runge-Kutta)");
	gr2->Draw("AL*");

	app.Run();

	return 0;
}
