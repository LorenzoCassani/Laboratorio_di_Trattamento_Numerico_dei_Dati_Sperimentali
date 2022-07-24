#include "lib.h"

#include "TH1F.h"
#include "TApplication.h"
#include "TCanvas.h"

using namespace std;

int main(int argc, char **argv){

	TApplication app("app",&argc,argv); //Istanzio un oggetto TApplication

	Random rand(3); //Istanzio un (mio) oggetto Random

	rand.SetA(1664525);    //Valore parametro a
	rand.SetC(1013904223); //Valore parametro c
	rand.SetM(pow(2,31));  //Valore parametro m

	TH1F h("Intervallo = [0,1]","Generatore Lineare Congruenziale",100,0,1); //Istanzio un oggetto TH1F
	
	for(int i=0; i<100000; i++){
		double n = rand.Rand01(); //Genero un numero casuale
		h.Fill(n);		  //Metto il numero casuale nell'istogramma
	}

	h.Draw(); // disegna l'istogramma

	app.Run();

	return 0;
}
