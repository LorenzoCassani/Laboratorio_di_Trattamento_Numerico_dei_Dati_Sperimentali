#include "lib.h"

using namespace std;

int main(int argc, const char **argv){

	if(argc!=2){
		cerr << "Usage: " << argv[0] << " <precision>" << endl;
		return -1;
	}

	double a = 0.;
	double b = M_PI;

	double prec = atof(argv[1]);

	FunzioneBase *seno = new Seno();
	Integral *integrale = new Integral(a,b,prec,seno);

	cout << endl << "Integrale = " << setprecision(12) << integrale->Trapezi() << endl << endl;

	return 0;
}



