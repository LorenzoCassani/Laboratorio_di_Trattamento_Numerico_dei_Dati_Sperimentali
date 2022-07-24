#include "lib.h"

#include "TH1F.h"
#include "TH2F.h"
#include "TApplication.h"
#include "TCanvas.h"

#define N 10000

using namespace std;

int main(int argc, char **argv){

	TApplication app("app",&argc,argv);

	EsperimentoPrisma esp;

	VettoreDati theta0(N);
	VettoreDati theta1(N);
	VettoreDati theta2(N);

	VettoreDati dm1(N);
	VettoreDati dm2(N);

	VettoreDati n1(N);
	VettoreDati n2(N);

	VettoreDati A(N);
	VettoreDati B(N);

	for(int i=0; i<N; i++){

		esp.Esegui();
		esp.Analizza();

		theta0.SetComponent(i,esp.Gett0Misurato());
		theta1.SetComponent(i,esp.Gett1Misurato());
		theta2.SetComponent(i,esp.Gett2Misurato());

		dm1.SetComponent(i,esp.Getdm1Misurato());
		dm2.SetComponent(i,esp.Getdm2Misurato());

		n1.SetComponent(i,esp.Getn1Misurato());
		n2.SetComponent(i,esp.Getn2Misurato());

		A.SetComponent(i,esp.GetAMisurato());
		B.SetComponent(i,esp.GetBMisurato());
	}

	//Istogrammi per t0, t1, t2
	TH1F h1("t0hist","Distribuzione di #theta_{0}",100,esp.Gett0Input()-3*esp.Getsigmat(),esp.Gett0Input()+3*esp.Getsigmat());
	TH1F h2("t1hist","Distribuzione di #theta_{1}",100,esp.Gett1Input()-3*esp.Getsigmat(),esp.Gett1Input()+3*esp.Getsigmat());
	TH1F h3("t2hist","Distribuzione di #theta_{2}",100,esp.Gett2Input()-3*esp.Getsigmat(),esp.Gett2Input()+3*esp.Getsigmat());

	//Istogrammi per dm1, dm2
	TH1F h4("n1hist","Distribuzione di #delta_{m,1}",100,esp.Getdm1Input()-3*dm1.StdDev(),esp.Getdm1Input()+3*dm1.StdDev());
	TH1F h5("n2hist","Distribuzione di #delta_{m,2}",100,esp.Getdm2Input()-3*dm2.StdDev(),esp.Getdm2Input()+3*dm2.StdDev());
	TH2F h6("ncorr","Residui 2D",100,esp.Getdm1Input()-3*dm1.StdDev(),esp.Getdm1Input()+3*dm1.StdDev(),100,esp.Getdm2Input()-3*dm2.StdDev(),esp.Getdm2Input()+3*dm2.StdDev());

	//Istogrammi per n1, n2
	TH1F h7("dm1hist","Distribuzione di n_{1}",100,esp.Getn1Input()-3*n1.StdDev(),esp.Getn1Input()+3*n1.StdDev());
	TH1F h8("dm2hist","Distribuzione di n_{2}",100,esp.Getn2Input()-3*n2.StdDev(),esp.Getn2Input()+3*n2.StdDev());
	TH2F h9("dcorr","Residui n_{1,2}",100,esp.Getn1Input()-3*n1.StdDev(),esp.Getn1Input()+3*n1.StdDev(),100,esp.Getn2Input()-3*n1.StdDev(),esp.Getn2Input()+3*n1.StdDev());

	//Istogrammi per A, B
	TH1F h10("Ahist","Distribuzione di A",100,esp.GetAInput()-3*A.StdDev(), esp.GetAInput()+3*A.StdDev());
	TH1F h11("Bhist","Distribuzione di B",100,esp.GetBInput()-3*B.StdDev(), esp.GetBInput()+3*B.StdDev());
	TH2F h12("ABcorr","Residui A e B",100,esp.GetAInput()-3*A.StdDev(),esp.GetAInput()+3*A.StdDev(),100,esp.GetBInput()-3*B.StdDev(),esp.GetBInput()+3*B.StdDev());


	for(int i=0; i<N; i++){

		h1.Fill(theta0.GetComponent(i));
		h2.Fill(theta1.GetComponent(i));
		h3.Fill(theta2.GetComponent(i));

		h4.Fill(dm1.GetComponent(i));
		h5.Fill(dm2.GetComponent(i));
		h6.Fill(dm1.GetComponent(i),dm2.GetComponent(i));

		h7.Fill(n1.GetComponent(i));
		h8.Fill(n2.GetComponent(i));
		h9.Fill(n1.GetComponent(i),n2.GetComponent(i));

		h10.Fill(A.GetComponent(i));
		h11.Fill(B.GetComponent(i));
		h12.Fill(A.GetComponent(i),B.GetComponent(i));
	}

	cout << endl << "Simulazione dell'esperienza dello spettrometro a prisma" << endl;

	cout << endl << "Parte 1:" << endl;

	cout << endl << "Media t0 = " << setprecision(5) << theta0.Media() << " rad" << endl;
	cout << "Media t1 = " << setprecision(5) << theta1.Media() << " rad" << endl;
	cout << "Media t2 = " << setprecision(5) << theta2.Media() << " rad" << endl << endl;

	cout << endl << "Deviazione standard t0 = " << setprecision(1) << theta0.StdDev() << " rad" << endl;
	cout << "Deviazione standard t1 = " << setprecision(1) << theta1.StdDev() << " rad" << endl;
	cout << "Deviazione standard t2 = " << setprecision(1) << theta2.StdDev() << " rad" << endl << endl;


	cout << endl << "Parte 2:" << endl;

	cout << endl << "Media dm1 = " << setprecision(6) << dm1.Media() << " rad" << endl;
	cout << "Media dm2 = " << setprecision(7) << dm2.Media() << " rad" << endl;

	cout << endl << "Deviazione standard dm1 = " << setprecision(2) << dm1.StdDev() << " rad" << endl;
	cout << "Deviazione standard dm2 = " << setprecision(2) << dm2.StdDev() << " rad" << endl;

	cout << endl << "Coefficiente di correlazione: " << setprecision(2) << dm1.CoeffCorr(dm2)*100 << "%" << endl << endl;


	cout << endl << "Media n1 = " << setprecision(6) << n1.Media() << endl;
	cout << "Media n2 = " << setprecision(7) << n2.Media() << endl;

	cout << endl << "Deviazione standard n1 = " << setprecision(2) << n1.StdDev() << endl;
	cout << "Deviazione standard n2 = " << setprecision(2) << n2.StdDev() << endl;

	cout << endl << "Coefficiente di correlazione: " << setprecision(2) << n1.CoeffCorr(n2)*100 << "%" << endl << endl;


	cout << endl << "Media A = " << setprecision(6) << A.Media() << endl;
	cout << "Media B = " << setprecision(5) << B.Media() << " nm^2" << endl;

	cout << endl << "Deviazione standard A = " << setprecision(2) << A.StdDev() << endl;
	cout << "Deviazione standard B = " << setprecision(3) << B.StdDev() << " nm^2" << endl;

	cout << endl << "Coefficiente di correlazione: " << setprecision(2) << A.CoeffCorr(B)*100 << "%" << endl << endl;

	cout << endl << "FINE SIMULAZIONE (Premere Ctrl+c per uscire dal programma)" << endl;


	//Canvas per t0, t1, t2
	TCanvas *c1 = new TCanvas("c1","c1",1200,600);
	c1->Divide(3,1);

	c1->cd(1);
	h1.GetXaxis()->SetTitle("#theta_{0}[rad]");
	h1.Draw();
	
	c1->cd(2);
	h2.GetXaxis()->SetTitle("#theta_{1}[rad]");
	h2.Draw();
	
	c1->cd(3);
	h3.GetXaxis()->SetTitle("#theta_{2}[rad]");
	h3.Draw();

	//Canvas per dm1, dm2
	TCanvas *c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(3,1);

	c2->cd(1);
	h4.GetXaxis()->SetTitle("Errori #delta_{m,1} [rad]");
	h4.Draw();
	
	c2->cd(2);
	h5.GetXaxis()->SetTitle("Errori #delta_{m,2} [rad]");
	h5.Draw();
	
	c2->cd(3);
	h6.GetXaxis()->SetTitle("Errori #delta_{m,1} [rad]");
	h6.GetYaxis()->SetTitle("Errori #delta_{m,2} [rad]");
	h6.Draw();

	//Canvas per n1, n2
	TCanvas *c3 = new TCanvas("c3","c3",1200,600);
	c3->Divide(3,1);

	c3->cd(1);
	h7.GetXaxis()->SetTitle("n(#lambda_{1})");
	h7.Draw();
	
	c3->cd(2);
	h8.GetXaxis()->SetTitle("n(#lambda_{2})");
	h8.Draw();
	
	c3->cd(3);
	h9.GetXaxis()->SetTitle("n(#lambda_{1})");
	h9.GetYaxis()->SetTitle("n(#lambda_{2})");
	h9.Draw();

	//Canvas per A, B
	TCanvas *c4 = new TCanvas("c4","c4",1200,600);
	c4->Divide(3,1);

	c4->cd(1);
	h10.GetXaxis()->SetTitle("A");
	h10.Draw();
	
	c4->cd(2);
	h11.GetXaxis()->SetTitle("B[m^{2}]");
	h11.Draw();
	
	c4->cd(3);
	h12.GetXaxis()->SetTitle("A");
	h12.GetYaxis()->SetTitle("B[m^{2}]");
	h12.Draw();

	app.Run();

	return 0;
}
