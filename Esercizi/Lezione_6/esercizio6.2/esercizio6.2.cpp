#include "lib.h"

using namespace std;

int main(int argc, char** argv){

	if(argc<3){
		cerr << "Usage: " << argv[0] << " <N_points> <filename>" << endl;
		return -1;
	}

	unsigned int N = atoi(argv[1]);
	VettoreDati v(N,argv[2]);

	cout << "Before sorting..." << endl;
	v.Print("Before.dat");
	v.QuickSort();
        cout << "After sorting..."<<endl;
	v.Print("After.dat");

	return 0;
}



