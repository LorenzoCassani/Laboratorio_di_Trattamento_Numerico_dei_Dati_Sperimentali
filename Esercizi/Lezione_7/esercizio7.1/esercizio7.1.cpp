#include "lib.h"

using namespace std;

int main(int argc, const char **argv){

	if(argc!=2){
		cerr << "Usage: " << argv[0] << " <nstep>" << endl;
		return -1;
	}

	double a = 0.;
	double b = M_PI;

	int nstep = atoi(argv[1]);

	FunzioneBase *seno = new Seno();
	Integral *integrale = new Integral(a,b,seno);

	cout << endl << "Integrale = " << setprecision(12) << integrale->Simpson(nstep) << endl << endl;

	return 0;
}



