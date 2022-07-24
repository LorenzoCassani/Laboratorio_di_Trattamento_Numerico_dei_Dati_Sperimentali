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

void Particella::Print() const {

	cout << "Particella: m=" << m_massa << ", q=" << m_carica << endl;
}

//////////////////////////////////////////////////ELETTRONE///////////////////////////////////////////////////////////////////////////////////////////

//Metodi per l'elettrone, particella di massa e carica predefinita

Elettrone::Elettrone(): Particella(9.1093826E-31,-1.60217653E-19){ //Invoco il costruttore della classe base con i parametri opportuni
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

double Parabola::Eval(double x) const {

	return m_a*pow(x,2)+m_b*x+m_c;
}

//////////////////////////////////////////////////ENRA////////////////////////////////////////////////////////////////////////////////////////////////

ENRA::ENRA(){
}

double ENRA::Eval(double x) const {

	return sin(x)-x*cos(x); 
}

//////////////////////////////////////////////////SENO////////////////////////////////////////////////////////////////////////////////////////////////

Seno::Seno(){
}

double Seno::Eval(double x) const {

	return sin(x);
}

//////////////////////////////////////////////////GAUSSIANA///////////////////////////////////////////////////////////////////////////////////////////

Gaussiana::Gaussiana(double mu, double sigma){

	m_mu = mu;
	m_sigma = sigma;
}

double Gaussiana::Eval(double x) const {

	return (exp(-(x-m_mu)*(x-m_mu)/(2*pow(m_sigma,2))))/(m_sigma*sqrt(2*M_PI));
}

//////////////////////////////////////////////////CARATTERISTICA//////////////////////////////////////////////////////////////////////////////////////

double Caratteristica::Eval(double x) const {

	if(x<=m_raggio) return 1;
	else return 0;
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
	m_v= new double[N];
	for(int i=0; i<N; ++i) m_v[i]=0; //Inizializzo tutti gli elementi del vettore a 0
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

VettoreDati::VettoreDati(unsigned int N): Vettore(N){
}

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

//////////////////////////////////////////////////VETTORE LINEARE/////////////////////////////////////////////////////////////////////////////////////

VettoreLineare::VettoreLineare(unsigned int N): Vettore(N){
}

VettoreLineare VettoreLineare::operator+(const VettoreLineare& b) const{
	VettoreLineare r(GetN());
	if(GetN()!=b.GetN()){
		cerr << "Errore: somma di vettori di dimensione " << GetN() << " e " << b.GetN() << endl;
		exit(-2);
	}
	for(unsigned int i=0; i<GetN(); i++)
		r.SetComponent(i,GetComponent(i)+b.GetComponent(i));
	return r;
}

double VettoreLineare::operator*(const VettoreLineare& b) const{
	double r=0;
	if(GetN()!=b.GetN()){
		cerr << "Errore: prodotto scalare di vettori di dimensione " << GetN() << " e " << b.GetN() << endl;
		exit(-2);
	}
	for(unsigned int i=0; i<GetN(); i++)
		r+=GetComponent(i)*b.GetComponent(i);
	return r;
}

VettoreLineare VettoreLineare::operator*(const double alpha) const{
	VettoreLineare r(GetN());
	for(unsigned int i=0; i<GetN(); i++)
		r.SetComponent(i,GetComponent(i)*alpha);
	return r;
}

VettoreLineare VettoreLineare::Versore(){
	VettoreLineare vers(GetN());
	for(unsigned int i=0; i<GetN(); i++)
		vers.SetComponent(i,GetComponent(i)/abs(GetComponent(i)));
	return vers;
}

double VettoreLineare::Modulo() const{
	VettoreLineare v(*this);
	return sqrt(v*v);
}

double VettoreLineare::Cos(const VettoreLineare& b) const{
	VettoreLineare a(*this);
	return (a * b)/(a.Modulo()*b.Modulo());
}

//////////////////////////////////////////////////RANDOM//////////////////////////////////////////////////////////////////////////////////////////////

double Random::Rand01(){

	double risultato = ((m_a*(m_seed)+m_c)%m_m);
	m_seed = risultato;
	return risultato/(m_m-1);
}

double Random::RandEsp(double rate){

	double s = Rand01();

	return -rate*log(Rand01());
}

double Random::RandGauss(double mu, double sigma){

	double s = this->Rand01();
	double t = this->Rand01();
	double g = sqrt(-2.*log(s))*cos(2*M_PI*t);

	return mu + sigma*g;
}

//////////////////////////////////////////////////ACCEPT REJECT///////////////////////////////////////////////////////////////////////////////////////

AcceptReject::AcceptReject(FunzioneBase *f, double a, double b, double max){

	srand(time(NULL)); //Inizializzo srand per estrarne seed casuali

	m_a = a;
	m_b = b;
	m_f = f;
	m_max = max;

	s = new Random (rand()%10+1); //s è un random generato uniformemente in [0,1]
	t = new Random (rand()%10+1); //t è un random generato uniformemente in [0,1]

	s->SetA(1664525);
	s->SetC(1013904223);
	s->SetM(pow(2,31));

	t->SetA(1664525);
	t->SetC(1013904223);
	t->SetM(pow(2,31));
}

AcceptReject::~AcceptReject(){ //Per evitare memory leak

	delete s;
	delete t;
}

double AcceptReject::RandGauss(){

	double x, y, Y;

	do{
		x = m_a+(m_b-m_a)*s->Rand01();
		y = m_max*t->Rand01();
		Y = m_f->Eval(x);

	}while(y>=Y);

	return x;
}

//////////////////////////////////////////////////INTEGRAL////////////////////////////////////////////////////////////////////////////////////////////

Integral::Integral(double a, double b, FunzioneBase *f): m_generatore(1){

	m_generatore.SetA(1664525);    //Valore parametro a
	m_generatore.SetC(1013904223); //Valore parametro c
	m_generatore.SetM(pow(2,31));  //Valore parametro m

	m_integrand = f;
	m_a = min(a,b);
	m_b = max(a,b);
	if(a>b) m_sign = -1;
	else m_sign = 1;
}

Integral::Integral(double a, double b, double prec, FunzioneBase *f): m_generatore(1){

	m_generatore.SetA(1664525);    //Valore parametro a
	m_generatore.SetC(1013904223); //Valore parametro c
	m_generatore.SetM(pow(2,31));  //Valore parametro m

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

//MONODIMENSIONALE
double Integral::Media(int n){

	double somma=0;
	double sommaq=0;

	for(int i=0; i<n; i++){
		double x = m_a+(m_b-m_a)*m_generatore.Rand01();
		somma+=m_integrand->Eval(x);
		sommaq+=pow(m_integrand->Eval(x),2);
	}

	m_integral = (m_b-m_a)*1/n*somma;
	m_errorMedia = (m_b-m_a)*sqrt(((sommaq/n)-pow((somma/n),2))/n);
	m_kMedia = (m_b-m_a)*sqrt((sommaq/n)-pow((somma/n),2));
 
	return m_integral;		
}

//MONODIMENSIONALE
double Integral::HitOrMiss(int n, double max){

	double A = (m_b-m_a)*max;
	double phit;

	int conta=0;
	for(int i=0; i<n; i++){
		double x = m_a+(m_b-m_a)*m_generatore.Rand01();
		double y = max*m_generatore.Rand01();
		double Y = m_integrand->Eval(x);
		if(Y<y) i--; //Se y<f(x), accetto il punto
		conta++;
	}

	m_integral = A*n/conta;
	phit = m_integral/A;
	m_errorHitOrMiss = A*sqrt((phit*(1-phit))/n);
	m_kHitOrMiss = sqrt(phit*(1-phit))*A;

	return m_integral;
}

//MULTIDIMENSIONALE
double Integral::HitOrMiss(int n){

	double *v = new double[m_dim];
	double Ns=0;

	for(int i=0; i<n; i++){

		double d=0;

		for(int j=0; j<m_dim; j++){
			v[j] = m_a+(m_b-m_a)*m_generatore.Rand01();
			d+=pow(v[j],2);
		}

		d = sqrt(d); //Calcolo la distanza
		Ns+=m_integrand->Eval(d); //Valuto la funzione caratteristica
	}

	return pow((m_b-m_a),m_dim)*(Ns/n);
}

//////////////////////////////////////////////////OSCILLATORE ARMONICO////////////////////////////////////////////////////////////////////////////////

OscillatoreArmonico::OscillatoreArmonico(double omega){
	m_omega = omega;
}

VettoreLineare OscillatoreArmonico::Eval(double t, const VettoreLineare& x) const{

	VettoreLineare v(2);
	v.SetComponent(0,x.GetComponent(1));		       //Setta la velocità iniziale in posizione 0 (x.GetComponent(1) = v0)
	v.SetComponent(1,(-pow(m_omega,2)*x.GetComponent(0))); //Setta l'accelerazione in posizione 1 (x.GetComponent(0) = x0)
	return v;
}

//////////////////////////////////////////////////OSCILLATORE FORZATO////////////////////////////////////////////////////////////////////////////////

OscillatoreForzato::OscillatoreForzato(double omega0, double omega, double alpha){
	m_omega0 = omega0;
	m_omega  = omega;
	m_alpha  = alpha;
}

VettoreLineare OscillatoreForzato::Eval(double t, const VettoreLineare& x) const{

	VettoreLineare v(2);
	v.SetComponent(0,x.GetComponent(1)); //Setta la velocità iniziale in posizione 0 (x.GetComponent(1) = v0)
	v.SetComponent(1,-pow(m_omega0,2)*x.GetComponent(0)-m_alpha*x.GetComponent(1)+sin(m_omega*t)); ////Setta l'accelerazione in posizione 1 (x.GetComponent(0) = x0)

	return v;
}

//////////////////////////////////////////////////PENDOLO/////////////////////////////////////////////////////////////////////////////////////////////

Pendolo::Pendolo(double l){
	m_l = l;
}

VettoreLineare Pendolo::Eval(double t, const VettoreLineare& x) const{

	const double g=9.8; //Accelerazione di gravità

	VettoreLineare v(2);
	v.SetComponent(0,x.GetComponent(1));		   //Setta la velocità angolare iniziale in posizione 0 (x.GetComponent(1) = omega0)
	v.SetComponent(1,-(g/m_l)*sin(x.GetComponent(0))); //Setta l'accelerazione angolare in posizione 1 (x.GetComponent(0) = A)
	return v;
}

//////////////////////////////////////////////////SFERA///////////////////////////////////////////////////////////////////////////////////////////////

Sfera::Sfera(double eta, double rho, double rho0, double R){
	m_eta = eta;
	m_rho = rho;
	m_rho0 = rho0;
	m_R = R;
}

VettoreLineare Sfera::Eval(double t, const VettoreLineare& x) const{

	const double g=9.81; //Accelerazione di gravità

	VettoreLineare v(2);
	v.SetComponent(0,x.GetComponent(1)); //Setta la velocità iniziale in posizione 0 (x.GetComponent(1) = v)
	v.SetComponent(1,-9/2.*m_eta/(m_rho*pow(m_R,2))*x.GetComponent(1)+(1-m_rho0/m_rho)*g); //Setta l'accelerazione in posizione 1 (x.GetComponent(1) = v)
	return v;
}

//////////////////////////////////////////////////EULERO//////////////////////////////////////////////////////////////////////////////////////////////

VettoreLineare Eulero::Passo(double t, const VettoreLineare& x, double h, FunzioneVettorialeBase *f) const{

	VettoreLineare a(x.GetN());

	a = x + f->Eval(t,x)*h;

	return a;
}

//////////////////////////////////////////////////RUNGE KUTTA/////////////////////////////////////////////////////////////////////////////////////////

VettoreLineare RungeKutta::Passo(double t, const VettoreLineare& x, double h, FunzioneVettorialeBase *f) const{

	VettoreLineare k1(x.GetN()),k2(x.GetN()),k3(x.GetN()),k4(x.GetN());

	k1 = f->Eval(t,x);
	k2 = f->Eval(t+h/2,x+k1*(h/2));
	k3 = f->Eval(t+h/2,x+k2*(h/2));
	k4 = f->Eval(t+h,x+k3*h);

	return x + (k1+k2*2+k3*2+k4)*(h/6);
}

//////////////////////////////////////////////////ESPERIMENTO PRISMA//////////////////////////////////////////////////////////////////////////////////

EsperimentoPrisma::EsperimentoPrisma() :
	m_generatore(1),
	m_lambda1(579.1E-9),
	m_lambda2(404.7E-9),
	m_alpha(60.*M_PI/180.),
	m_sigmat(0.3E-3),
	m_A_input(2.7),
	m_B_input(60000E-18)
{
	m_generatore.SetA(1664525);    //Valore parametro a
	m_generatore.SetC(1013904223); //Valore parametro c
	m_generatore.SetM(pow(2,31));  //Valore parametro m

	//Calcolo degli indici di rifrazione attesi
	m_n1_input = sqrt(m_A_input+m_B_input/(m_lambda1*m_lambda1));
	m_n2_input = sqrt(m_A_input+m_B_input/(m_lambda2*m_lambda2));
	//Calcolo dei valori attesi degli angoli misurati
	m_t0_input = M_PI/2.; //theta 0 è arbitrario
	m_dm1_input = 2.*asin(m_n1_input*sin(0.5*m_alpha))-m_alpha;
	m_t1_input = m_t0_input+m_dm1_input;
	m_dm2_input = 2.*asin(m_n2_input*sin(0.5*m_alpha))-m_alpha;
	m_t2_input = m_t0_input+m_dm2_input;
}

void EsperimentoPrisma::Esegui(){ //Misura sperimentale:

	m_t0_misurato = m_generatore.RandGauss(m_t0_input,m_sigmat); //1. dell'angolo corrispondente al fascio non deflesso in assenza del prisma

	m_t1_misurato = m_generatore.RandGauss(m_t1_input,m_sigmat); //2. dell'angolo corrispondente alla deviazione minima della riga del giallo

	m_t2_misurato = m_generatore.RandGauss(m_t2_input,m_sigmat); //3. dell'angolo corrispondente alla deviazione minima della riga del viola
}

void EsperimentoPrisma::Analizza(){ //Analisi dati

	m_dm1_misurato = m_t1_misurato-m_t0_misurato;
	m_dm2_misurato = m_t2_misurato-m_t0_misurato;

	m_n1_misurato = sin((m_dm1_misurato+m_alpha)/2)/sin(m_alpha/2);
	m_n2_misurato = sin((m_dm2_misurato+m_alpha)/2)/sin(m_alpha/2);

	m_A_misurato = (pow(m_lambda2*m_n2_misurato,2)-pow(m_lambda1*m_n1_misurato,2))/(m_lambda2*m_lambda2-m_lambda1*m_lambda1);
	m_B_misurato = (pow(m_n2_misurato,2)-pow(m_n1_misurato,2))/(pow(m_lambda2,-2)-pow(m_lambda1,-2));
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

