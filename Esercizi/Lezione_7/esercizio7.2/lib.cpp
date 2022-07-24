#include "lib.h" 

using namespace std;

/*************************************************METODI CLASSI**************************************************************************************/

//////////////////////////////////////////////////POSIZIONE///////////////////////////////////////////////////////////////////////////////////////////

//Costruttore di default
Posizione::Posizione(){

	m_x=0.;
	m_y=0.;
	m_z=0.;
}

//Costruttore a partire da una terna cartesiana
Posizione::Posizione(double x, double y, double z){

	m_x=x;
	m_y=y;
	m_z=z;
}

//Distruttore (può essere vuoto)
Posizione::~Posizione(){
}

//Coordinate cartesiane
double Posizione::GetX() const {
	return m_x;
}
double Posizione::GetY() const {
	return m_y;
}
double Posizione::GetZ() const {
	return m_z;
}

void Posizione::SetX(double x){
	m_x=x;
}

void Posizione::SetY(double y){
	m_y=y;
}

void Posizione::SetZ(double z){
	m_z=z;
}

//Coordinate sferiche
double Posizione::R() const {
	return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
}
double Posizione::Phi() const {
	return atan2(m_y,m_x);
}
double Posizione::Theta() const {
	return acos(m_z/R());
}

//Raggio delle coordinate cilindriche
double Posizione::Rho() const {
	return sqrt(m_x*m_x+m_y*m_y);
}

//Distanza da un altro punto
double Posizione::Distanza(const Posizione& b) const {
	return sqrt( pow(GetX()-b.GetX(),2)
		    +pow(GetY()-b.GetY(),2)
	            +pow(GetZ()-b.GetZ(),2) );
}

//////////////////////////////////////////////////PARTICELLA//////////////////////////////////////////////////////////////////////////////////////////

//Metodi per la classe base

Particella::Particella(double massa, double carica){

	m_massa=massa;
	m_carica=carica;
}

Particella::~Particella(){
}

void Particella::Print() const {

	cout << "Particella: m=" << m_massa << ", q=" << m_carica << endl;
}

//////////////////////////////////////////////////ELETTRONE///////////////////////////////////////////////////////////////////////////////////////////

//Metodi per l'elettrone, particella di massa e carica predefinita

Elettrone::Elettrone(): Particella(9.1093826E-31,-1.60217653E-19){ //Invoco il costruttore della classe base con i parametri opportuni
}

Elettrone::~Elettrone(){
}

void Elettrone::Print() const {

	cout << "Elettrone: m=" << m_massa << ", q=" << m_carica << endl;
}

//////////////////////////////////////////////////CORPO CELESTE///////////////////////////////////////////////////////////////////////////////////////

//Metodi per il corpo celeste

CorpoCeleste::CorpoCeleste(string nome, double massa, double raggio): Particella(massa,0){ //Costruttore della classe madre con carica nulla

	m_nome=nome;
	m_raggio=raggio;
}

CorpoCeleste::~CorpoCeleste(){
}

void CorpoCeleste::Print() const {

    	cout << m_nome << ": m=" << m_massa << ", R=" << m_raggio << endl;
}

double CorpoCeleste::PotenzialeGravitazionale(Posizione r){

	return (-(G*m_massa)/m_posizione.Distanza(r));
}

//////////////////////////////////////////////////CAMPO VETTORIALE////////////////////////////////////////////////////////////////////////////////////

CampoVettoriale::CampoVettoriale(const Posizione& r): Posizione(r){ //Crea un vettore nullo nella posizione indicata

	m_x=0.;
	m_y=0.;
	m_z=0.;
}

CampoVettoriale::~CampoVettoriale(){
}

double CampoVettoriale::Modulo() const {

  	return sqrt(pow(m_x,2)+pow(m_y,2)+pow(m_z,2));
}

void CampoVettoriale::Somma(const CampoVettoriale& v){

	m_x=m_x+v.GetFX();
	m_y=m_y+v.GetFY();
	m_z=m_z+v.GetFZ();
}

//////////////////////////////////////////////////PUNTO MATERIALE/////////////////////////////////////////////////////////////////////////////////////

PuntoMateriale::PuntoMateriale(double massa, double carica, double x, double y, double z): Particella(massa,carica), Posizione(x,y,z){
}

PuntoMateriale::~PuntoMateriale(){
}

CampoVettoriale PuntoMateriale::CampoElettrico(const Posizione& r) const {

	CampoVettoriale E(r);
	E.SetFX((m_carica/(4*M_PI*epsilon_0))*(((r.GetX()-GetX())/pow(Distanza(r),3))));
	E.SetFY((m_carica/(4*M_PI*epsilon_0))*(((r.GetY()-GetY())/pow(Distanza(r),3))));
	E.SetFZ((m_carica/(4*M_PI*epsilon_0))*(((r.GetZ()-GetZ())/pow(Distanza(r),3))));

	return E;
}

CampoVettoriale PuntoMateriale::CampoGravitazionale(const Posizione& r) const {

	CampoVettoriale g(r);
	g.SetFX(-G*m_massa*((r.GetX()-GetX())/pow(Distanza(r),3)));
	g.SetFY(-G*m_massa*((r.GetY()-GetY())/pow(Distanza(r),3)));
	g.SetFZ(-G*m_massa*((r.GetZ()-GetZ())/pow(Distanza(r),3)));

	return g;
}

//////////////////////////////////////////////////PARABOLA////////////////////////////////////////////////////////////////////////////////////////////

Parabola::Parabola(){

	m_a=0.;
	m_b=0.;
	m_c=0.;

}

Parabola::Parabola(double a, double b, double c){

	m_a=a;
	m_b=b;
	m_c=c;
}

Parabola::~Parabola(){
}

double Parabola::Eval(double x) const {

	return m_a*pow(x,2)+m_b*x+m_c;
}

//////////////////////////////////////////////////SENO////////////////////////////////////////////////////////////////////////////////////////////////

Seno::Seno(){
}

Seno::~Seno(){
}

double Seno::Eval(double x) const {

    return sin(x);
}

//////////////////////////////////////////////////ENRA////////////////////////////////////////////////////////////////////////////////////////////////

ENRA::ENRA(){
}

ENRA::~ENRA(){
}

double ENRA::Eval(double x) const {

	return sin(x)-x*cos(x); 
}

//////////////////////////////////////////////////BISEZIONE///////////////////////////////////////////////////////////////////////////////////////////

Bisezione::Bisezione(){

	m_a=0.;
	m_b=0.;
	m_prec=0.;
	m_f=NULL;
}

Bisezione::Bisezione(double a, double b, double prec, const FunzioneBase* f){

	m_a=a;
	m_b=b;
	m_prec=prec;
	m_f=f;
}

Bisezione::~Bisezione(){
}

double Bisezione::CercaZeri(double xmin, double xmax){

	double c=xmin; //midpoint

	if(sign(m_f->Eval(xmin))==0){

		c=xmin;
		return c;
	}

	else if(sign(m_f->Eval(xmax))==0){

		c=xmax;
		return c;
	}

	else while(xmax-xmin >= m_prec){

		c = xmin + 0.5*(xmax-xmin);

		if(m_f->Eval(c)==0) break;

		else if(sign(m_f->Eval(xmin))*sign(m_f->Eval(c))<0) xmax=c;

		else if(sign(m_f->Eval(xmax))*sign(m_f->Eval(c))<0) xmin=c;
	}

	return c;
}

bool Bisezione::Trovato(){

	if(sign(m_f->Eval(m_a))==sign(m_f->Eval(m_b))){
		cerr << endl << "f(a) ha lo stesso segno di f(b)!" << endl;
		return false;
	}

	else {

		double xmin=m_a, xmax=m_b, m;
		bool risp;

		if(m_f->Eval(xmin)==0) return true;

		else if (m_f->Eval(xmax)==0) return true;

		else for(int i=0, j=1; xmax-xmin >= m_prec; i++){

			m=xmin+((xmax-xmin)/2);

                	if(m_f->Eval(m)==0) return true;

                	else if(sign(m_f->Eval(m))==sign(m_f->Eval(xmin))) xmin=m;

                	else xmax=m;

                	if(i==100){
				cerr << endl << "L'algoritmo ha iterato " <<i*j <<" volte, continuare? (1: Sì   0: No) " << endl;
				cin >> risp;
				if(risp==false) return false;
				else{
					i=0;
					j++;
				}
			}
		}
	}
	return true;
}

double Bisezione::Incertezza(){

	double xmin=m_a, xmax=m_b, m;

	if(m_f->Eval(xmin)==0) return 0;

	else if(m_f->Eval(xmax)==0) return 0;

	else while(xmax-xmin >= m_prec){

		m=xmin+((xmax-xmin)/2);

		if(m_f->Eval(m)==0) return 0;

		else if(sign(m_f->Eval(m))*sign(m_f->Eval(xmin))>0) xmin=m;

		else xmax=m;
	}
	return (xmax-xmin)/2;
}

//////////////////////////////////////////////////VETTORE/////////////////////////////////////////////////////////////////////////////////////////////

Vettore::Vettore(){
}

Vettore::Vettore(unsigned int N){


	m_N=N;
	m_v = new double[N];
	for(int i=0; i<N; ++i) m_v[i]=0; //Inizializzo tutti gli elementi del vettore a 0
}

Vettore::~Vettore(){
}

void Vettore::SetComponent(unsigned int i, double val){

	while(i >= m_N){

		cerr << "La posizione specificata è non è compatibile con la dimensione dell'array!" << endl;
		return;
	}

	m_v[i]=val;
}

double Vettore::GetComponent(unsigned int i) const {

	while(i >= m_N){

		cerr << "La posizione specificata è non è compatibile con la dimensione dell'array!" << endl;
		return -1;
	}

	return m_v[i];
}

Vettore::Vettore(const Vettore& V){

	m_N = V.m_N;
	m_v = new double[m_N];
	for(int i=0; i<m_N; i++) m_v[i]=V.m_v[i];
}

Vettore& Vettore::operator=(const Vettore& V){

	m_N = V.m_N;
	if(m_v) delete[] m_v;
	m_v = new double[m_N];
	for(int i=0; i<m_N; i++) m_v[i]=V.m_v[i];
	return *this;
}

//////////////////////////////////////////////////VETTORE DATI////////////////////////////////////////////////////////////////////////////////////////

VettoreDati::VettoreDati(unsigned int N, const char* filename){

	m_N = N;
	m_v = new double[m_N];
	ifstream in(filename);
	if(!in)
		cerr << "Cannot open file" << filename << endl;
	else{
		for(unsigned int i=0; i<N; i++){
			in >> m_v[i];
			if(in.eof()){
				cerr << "Stop reading after " << i << " entries!" << endl;
				break;
			}
		}
	}
}

void VettoreDati::Print(){

	int width = int(log10(m_N)+1);
	for(unsigned int i=0; i<m_N; i++){
		cout << setw(width) << i << ") " << m_v[i] << endl;
	}
}

void VettoreDati::Print(const char* filename){

	ofstream out(filename);
	int width = int(log10(m_N)+1);
	for(unsigned int i=0; i<m_N; i++){
		out << setw(width) << i << ") " << m_v[i] << endl;
	}
	out.close();
}

void VettoreDati::SelectionSort(){

	for(int j=0; j<m_N; j++){
		int iMin=j;
		for(int i=j; i<m_N; i++){
			if(m_v[i]<m_v[iMin])
				iMin=i;
		}
		if(iMin!=j)
			swap(m_v[j],m_v[iMin]);
	}
}

void VettoreDati::QuickSort(){

	QuickSort(0,m_N-1);
}

void VettoreDati::QuickSort(unsigned int primo, unsigned int ultimo){

	if( primo>ultimo or ultimo>=m_N ) return; //Controllo di validità degli estremi
	if( ultimo-primo <= 1 ){
		//Intervallo di 1 o 2 elementi, viene ordinato
		if( m_v[primo]>m_v[ultimo] ) swap(m_v[primo],m_v[ultimo]);
		return;
	}
	double pivot = m_v[(primo+ultimo)/2]; //Elemento di pivot a metà dell'intervallo
	unsigned int basso=primo, alto=ultimo;
	while(basso<alto){
		while( m_v[basso]<pivot ) basso++; //Cerca un elemento > pivot
		while( m_v[alto]>pivot ) alto--;   //Cerca un elemento < pivot
		if(basso<alto) {swap(m_v[basso],m_v[alto]); basso++;}
	}
	QuickSort(primo,basso-1);
	QuickSort(basso,ultimo);
}

double VettoreDati::Media(){

	double somma=0;

	for(int i=0; i<m_N; i++) somma+=m_v[i];

	return somma/m_N;
}

double VettoreDati::Var(){

	double somma=0;
	double sommaq=0;

	for(int i=0; i<m_N; i++){
		somma+=m_v[i];
		sommaq+=pow(m_v[i],2);
	}

	return m_N/(m_N-1)*((sommaq/m_N)-pow((somma/m_N),2));
}

double VettoreDati::StdDev(){

	double somma=0;
	double sommaq=0;

	for(int i=0; i<m_N; i++){
		somma+=m_v[i];
		sommaq+=pow(m_v[i],2);
	}

	return sqrt(m_N/(m_N-1)*((sommaq/m_N)-pow((somma/m_N),2)));
}

double VettoreDati::Mediana(){

	QuickSort(); //Mediana: il numero a metà della serie ordinata dei valori del vettore

	double mediana=0;

	if(m_N%2==0) mediana = (m_v[m_N/2]+m_v[(m_N/2)+1])/2.;
	if(m_N%2!=0) mediana = m_v[(m_N/2)+1];

	return mediana;
}

double VettoreDati::CoeffCorr(VettoreDati& V){

	if(m_N!=V.GetN()){
		cerr << endl << "N deve essere lo stesso per entrambi i vettori!" << endl;
		return-1;
	}

	else{

		double somma=0;
		double covarianza=0;

		for(int i=0; i<m_N; i++) somma+=((m_v[i]-Media())*(V.GetComponent(i)-V.Media()));

		covarianza=(somma/m_N);

		return (covarianza/(StdDev()* V.StdDev()));
	}
}

//////////////////////////////////////////////////INTEGRAL////////////////////////////////////////////////////////////////////////////////////////////

Integral::Integral(double a, double b, FunzioneBase *f){

	m_integrand = f;
	m_a = min(a,b);
	m_b = max(a,b);
	if(a>b) m_sign = -1;
	else m_sign = 1;
}

Integral::Integral(double a, double b, double prec, FunzioneBase *f){

	m_prec = prec;
	m_integrand = f;
	m_a = min(a,b);
	m_b = max(a,b);
	if(a>b) m_sign = -1;
	else m_sign = 1;
}

double Integral::Midpoint(int nstep){

	m_sum = 0.;
	m_h = (m_b-m_a)/nstep;

	for(int i=0; i<nstep; ++i){
		double x = m_a+(i+0.5)*m_h;
		m_sum += m_integrand->Eval(x);
	}

	m_integral = m_sign*m_sum*m_h;
	return m_integral;
}

double Integral::Simpson(int nstep){

	if(nstep%2!=0){
		cout << "nstep deve essere pari (come prevede il metodo Simpson)" << endl;
		return -1;
	}

	m_sum = 0.;
	m_h = (m_b-m_a)/nstep;

	for(int k=1; k<nstep-2; k=k+2){
		double x = m_a+k*m_h;
		m_sum += 4./3.*m_integrand->Eval(x)+2./3.*m_integrand->Eval(x+m_h);
	}

	m_sum += 4./3.*m_integrand->Eval(m_b-m_h)+1./3.*(m_integrand->Eval(m_a)+m_integrand->Eval(m_b));

	m_integral = m_sign*m_sum*m_h;
        return m_integral;
}

double Integral::Trapezi(){

	double delta=INT_MAX;

	for(unsigned int k=0; delta>m_prec; k++){
		if(k==0){
			m_sum = (m_integrand->Eval(m_a)+m_integrand->Eval(m_b))/2.;
			m_integral = m_sign*m_sum*(m_b-m_a);
		}
		else{
			m_h = (m_b-m_a)/(pow(2,k-1));
			double sum_f, x;
			double x_1 = m_a+m_h/2.;
			for(unsigned int i=1; i<=(pow(2,k-1)); i++){
				if(i==1) sum_f = m_integrand->Eval(x_1);
				else{
					if(i==(pow(2,k-1))) x = m_b-m_h/2.;
					else x = x_1+(i-1)*m_h;
					sum_f = sum_f+m_integrand->Eval(x);
				}
			}
			m_sum = m_sum+sum_f;
			double I_new = m_sign*m_sum*(m_b-m_a)/(pow(2,k));
			delta = abs(I_new-m_integral);
			m_integral = I_new;
		}
	}
	return m_integral;
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

	for(int i=0; i<dim; i++) somma+=arr[i];
	
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

	for(int i=0; i<n; i++) sommaq+=pow(arr[i],2);

	return sqrt(sommaq);
}

//////////////////////////////////////////////////LEZIONE_2///////////////////////////////////////////////////////////////////////////////////////////

struct vettore read(unsigned int N, const char* filename){

	struct vettore vec;	//Dichiarazione di vec come una struttura Vettore
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

void print(const struct vettore vec, const char* filename){

	ofstream out(filename);
	int width = int(log10(vec.N)+1);
	for(unsigned int i=0; i<vec.N; i++){
		out << setw(width) << i << ") " << vec.v[i] << endl;
	}
	out.close();
}

void selection_sort(struct vettore vec){

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

//////////////////////////////////////////////////LEZIONE_5///////////////////////////////////////////////////////////////////////////////////////////

int sign(double x){

	if(x==0.)
		return 0;
	if(x>0)
		return 1;
	if(x<0)
		return -1;
}

