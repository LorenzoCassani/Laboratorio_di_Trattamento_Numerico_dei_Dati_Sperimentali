#include "lib.h" 

using namespace std;

/*************************************************METODI CLASSI**************************************************************************************/

//Costruttore di default
Posizione::Posizione(){

	m_R=0.;
	m_Phi=0.;
	m_Theta=0.;
}

//Costruttore a partire da una terna cartesiana
Posizione::Posizione(double x, double y, double z){

	m_R=sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	m_Phi=atan(y/x);
	m_Theta=acos(z/sqrt(pow(x,2)+pow(y,2)+pow(z,2)));
}

//Distruttore (può essere vuoto)
Posizione::~Posizione(){
}

//Coordinate cartesiane
double Posizione::X() const {
	return m_R*sin(m_Theta)*cos(m_Phi);
}
double Posizione::Y() const {
	return m_R*sin(m_Theta)*sin(m_Phi);
}
double Posizione::Z() const {
	return m_R*cos(m_Theta);
}

//Coordinate sferiche
double Posizione::R() const {
	return m_R;
}
double Posizione::Phi() const {
	return m_Phi;
}
double Posizione::Theta() const {
	return m_Theta;
}

//Raggio delle coordinate cilindriche
double Posizione::Rho() const {
	return m_R*sin(m_Theta);
}

//Distanza da un altro punto
double Posizione::Distanza(const Posizione& b) const {
	return sqrt( pow(X()-b.X(),2)
		    +pow(Y()-b.Y(),2)
	            +pow(Z()-b.Z(),2) );
}

/*************************************************FUNZIONI*******************************************************************************************/

//////////////////////////////////////////////////LEZIONE_1///////////////////////////////////////////////////////////////////////////////////////////

int Nscambi = 0; //Variabile globale: utilizzabile da tutte le funzioni di questo file

int numeroScambi(){ //Questa funzione può venire chiamata da una funzione esterna per consultare il numero di scambi effettuati

	return Nscambi;
}

void scambiaByValue(double a, double b){ //NON FUNZIONA

	static int Ns=0; //Variabile statica locale a questa funzione
	double appo=a;
	a=b;
	b=appo;
	Ns++;
	Nscambi++;
	cout << "Numero di scambi per valore effettuati: " << Ns << endl;
}

void scambiaByRef(double &a, double &b){

	static int Ns=0; //Variabile statica locale a questa funzione
	double appo=a;
	a=b;
	b=appo;
	Ns++;
	Nscambi++;
	cout << "Numero di scambi per referenza effettuati: " << Ns << endl;
}

void scambiaByPointer(double *a, double *b){

	static int Ns=0; //Variabile statica locale a questa funzione
	double appo=*a;
	*a=*b;
	*b=appo;
	Ns++;
	Nscambi++;
	cout << "Numero di scambi per puntatore effettuati: " << Ns << endl;
}


double media(double arr[], int dim){

	double somma=0;

	for(int i=0; i<dim; i++){
		somma+=arr[i];
	}
	
	return somma/dim;
}

double var(double arr[], int dim){

	double somma=0;
	double sommaq=0;

	for(int i=0; i<dim; i++){
		somma+=arr[i];
		sommaq+=pow(arr[i],2);
	}

	return dim/(dim-1)*((sommaq/dim)-pow((somma/dim),2));
}


double modulo(double a){

	return abs(a);
}

double modulo(double v1, double v2){
  
	double sommaq=0;
 
	sommaq=pow(v1,2)+pow(v2,2);
  	
	return sqrt(sommaq);
}

double modulo(double arr[], int n){

	double sommaq=0;

	for(int i=0; i<n; i++)
		sommaq+=pow(arr[i],2);

	return sqrt(sommaq);
}

//////////////////////////////////////////////////LEZIONE_2///////////////////////////////////////////////////////////////////////////////////////////

struct Vettore read(unsigned int N, const char* filename){

	struct Vettore vec;	//Dichiarazione di vec come una struttura Vettore
	vec.N = N;		//Lunghezza del vettore
	vec.v = new double [N]; //Allocazione della memoria per i dati
	ifstream in(filename);	//Apertura del file
	if(!in)			//Controllo che l'apertura abbia avuto successo
		cerr << "Cannot open file" << filename << endl;
	else{
		for(unsigned int i=0; i<N; i++){
			in >> vec.v[i];	//Leggo N valori dal file
			if(in.eof()){	//Controllo di non superare la fine del file
				cerr << "Stop reading after " << i << " entries!" << endl;
				break;
			}
		}
	}

	return vec;	//Restituisco la nuova struttura
}

void print(const struct Vettore vec, const char* filename){

	ofstream out(filename);
	int width = int(log10(vec.N)+1);
	for(unsigned int i=0; i<vec.N; i++){
		out << setw(width) << i << ") " << vec.v[i] << endl;
	}
	out.close();
}

void selection_sort(struct Vettore vec){

	for(int j=0; j<vec.N; j++){
		int iMin=j;
		for(int i=j; i<vec.N; i++){
			if(vec.v[i]<vec.v[iMin])
				iMin=i;
		}
		if(iMin!=j)
			swap(vec.v[j],vec.v[iMin]);
	}
}
