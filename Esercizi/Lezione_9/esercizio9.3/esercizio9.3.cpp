#include "lib.h"

#include "TApplication.h"
#include "TCanvas.h"
#include "TH1F.h"

using namespace std;

int main(int argc, char**argv){

	TApplication app("app",&argc,argv);

	int dim=3; //Dimensioni

	double raggio=1.; //Raggio della sfera n-dimensionale

	double esatto;	  //Valore esatto di Sn/Rn
	double calcolato; //Valore calcolato di Sn/Rn
	int N=100;	  //Ripetizioni del calcolo di Sn/Rn

	double somma, sommaq, media, devstd, devstdmedia;

	TH1F h1("10000 punti" ,"Rapporto S_{n}/R_{n} (n=3)",25,0.5,0.55);
	TH1F h2("40000 punti" ,"Rapporto S_{n}/R_{n} (n=3)",25,0.5,0.55);
	TH1F h3("90000 punti" ,"Rapporto S_{n}/R_{n} (n=3)",25,0.5,0.55);
	TH1F h4("160000 punti","Rapporto S_{n}/R_{n} (n=3)",25,0.5,0.55);

	FunzioneBase *f = new Caratteristica(raggio); //Funzione caratteristica con raggio=1
	Integral integrale(-raggio,raggio,f);	      //Integrale
	integrale.SetDim(dim);			      //Assegno la dimensione all'integrale


	cout << endl << "Andamento del rapporto Sn/Rn dove Sn è il volume di una sfera n-dimensionale di raggio 1 e Rn è il volume del cubo n-dimensionale di lato 2 ad essa circoscritto" << endl;


	//PARTE 1:

	cout << endl << "PARTE 1:" << endl;
	cout << endl << "Per n=3, " << N << " ripetizioni:" << endl;

	//10000 punti
	somma=0.;
	sommaq=0.; 

	for(int i=0; i<N; i++){

		calcolato = (integrale.HitOrMiss(10000))/pow(2*raggio,dim); // = Sn/Rn
		somma+=calcolato;
		h1.Fill(calcolato);
		sommaq+=pow(calcolato,2);
	}

	media = somma/N;
	devstd = sqrt((sommaq/N)-pow(media,2));
	devstdmedia = devstd/N;

	cout << endl << "10000  PUNTI: Media = " << setprecision(5) << media
	<< ", DevStd =  " << setprecision(2) << devstd << ", DevStd media = " << setprecision(2) << devstdmedia << endl;

	//40000 punti
	somma=0.;
	sommaq=0.; 

	for(int i=0; i<N; i++){

		calcolato = (integrale.HitOrMiss(40000))/pow(2*raggio,dim); // = Sn/Rn
		somma+=calcolato;
		h2.Fill(calcolato);
		sommaq+=pow(calcolato,2);
	}

	media = somma/N;
	devstd = sqrt((sommaq/N)-pow(media,2));
	devstdmedia = devstd/N;

	cout << endl << "40000  PUNTI: Media = " << setprecision(5) << media
	<< ", DevStd =  " << setprecision(2) << devstd << ", DevStd media = " << setprecision(2) << devstdmedia << endl;

	//90000 punti
	somma=0.;
	sommaq=0.; 

	for(int i=0; i<N; i++){

		calcolato = (integrale.HitOrMiss(90000))/pow(2*raggio,dim); // = Sn/Rn
		somma+=calcolato;
		h3.Fill(calcolato);
		sommaq+=pow(calcolato,2);
	}

	media = somma/N;
	devstd = sqrt((sommaq/N)-pow(media,2));
	devstdmedia = devstd/N;

	cout << endl << "90000  PUNTI: Media = " << setprecision(5) << media
	<< ", DevStd =  " << setprecision(2) << devstd << ", DevStd media = " << setprecision(2) << devstdmedia << endl;

	//160000 punti
	somma=0.;
	sommaq=0.; 

	for(int i=0; i<N; i++){

		calcolato = (integrale.HitOrMiss(160000))/pow(2*raggio,dim); // = Sn/Rn
		somma+=calcolato;
		h4.Fill(calcolato);
		sommaq+=pow(calcolato,2);
	}

	media = somma/N;
	devstd = sqrt((sommaq/N)-pow(media,2));
	devstdmedia = devstd/N;

	cout << endl << "160000  PUNTI: Media = " << setprecision(5) << media
	<< ", DevStd =  " << setprecision(2) << devstd << ", DevStd media = " << setprecision(2) << devstdmedia << endl;


	//PARTE 2:

	cout << endl << "PARTE 2:" << endl;
	cout << endl << "Per n = 2,...,5, " << N << " ripetizioni, 160000 punti:" << endl << endl;

	for(dim=2; dim<=5; dim++){

		integrale.SetDim(dim);

		calcolato = (integrale.HitOrMiss(160000))/pow(2*raggio,dim); // = Sn/Rn
		esatto = (pow(M_PI,dim/2.)/(tgamma(dim*0.5+1)))/pow(2.,dim);

		cout << "Sn/Rn calcolato (n=" << dim << "): " << setprecision(5) << calcolato << endl;
		cout << "Sn/Rn esatto    (n=" << dim << "): " << setprecision(5) << esatto << endl; //tgamma è contenuta in cmath
		cout << "Discrepanza: " << setprecision(2) << (abs(calcolato-esatto)/esatto)*100 << " %" << endl << endl;
	}

	delete f;

	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	c1->Divide(2,2);

	c1->cd(1);
	h1.Draw();

	c1->cd(2);
	h2.Draw();

	c1->cd(3);
	h3.Draw();

	c1->cd(4);
	h4.Draw();

	app.Run();

	return 0;
}
