#include "lib.h"

using namespace std;

namespace math {

   	const int one   = 1;
   	const int tw    = 2;
   	const int three = 3;
   	const int four  = 4;
   	const int five  = 5;
   	const int six   = 6;
   	const int seven = 7;
   	const int eight = 8;
   	const int nine  = 9;

   	int add (const int x, const int y) {return x+y;}
   	int subtract (const int x, const int y) {return x-y;}
}

namespace smath {

   	const int one   = 4;
   	const int two   = 9;
   	const int three = 3;
   	const int four  = 1;
   	const int five  = 8;
   	const int six   = 2;
   	const int seven = 5;
   	const int eight = 6;
   	const int nine  = 7;

   	int add (const int x, const int y) {return x-y;}
   	int subtract (const int x, const int y) {return x+y;}
}

using namespace math;
using namespace smath;

int main(){

   	cout << endl << "Matematica corretta: sommo e sottraggo due numeri" << endl;
   	cout << math::one << "+" << math::nine
   	<< "=" << math::add(math::one,math::nine) << endl;
   	cout << math::seven << "-" << math::three
   	<< "=" << math::subtract(math::seven,math::three) << endl;

   	cout << endl << "Matematica strana:" << endl;
	cout << smath::one << "+" << smath::nine
   	<< "=" << smath::add(smath::one,smath::nine) << endl;
   	cout << smath::seven << "-" << smath::three
   	<< "=" << smath::subtract(smath::seven,smath::three) << endl;

	return 0;
}



