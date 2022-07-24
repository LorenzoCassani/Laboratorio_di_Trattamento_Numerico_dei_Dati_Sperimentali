#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <climits>

//Definizioni delle costanti fisiche:
#define G 6.67E-11 // N*m^2/Kg^2
#define epsilon_0 8.85418781762E-12 // F/m
//N.B: attenzione a non inizializzare variabili con lo stesso nome di una delle costanti!

using namespace std;

/*************************************************CLASSI*********************************************************************************************/

#ifndef __Posizione_h__
#define __Posizione_h__

class Posizione{

public:
	//Costruttori
	Posizione();
	Posizione(double x, double y, double z);
	//Distruttore
	~Posizione(){}
	//Metodi
	double GetX() const;
	double GetY() const;
	double GetZ() const;	
   	void SetX(double); //Coordinate cartesiane
   	void SetY(double);
   	void SetZ(double);

	double R() const; //Coordinate sferiche
	double Phi() const;
	double Theta() const;
	double Rho() const; //Raggio delle coordinate cilindriche
	double Distanza(const Posizione&) const; //Distanza da un altro punto

private:
	double m_x, m_y, m_z;

};

#endif


//Classe per una generica particella
//Definisce i metodi che ogni particella deve implementare e delle funzioni di utilizzo generale
#ifndef __Particella_h__
#define __Particella_h__

class Particella{

public:
	//Costruttore
	Particella(double massa, double carica);
	//Distruttore
	~Particella(){}
	//Metodi
	double GetMassa() const {return m_massa;}
	double GetCarica() const {return m_carica;}
	//Esempio di metodo virtuale: vogliamo che le classi derivate possano sovrascriverlo
	virtual void Print() const;

protected:
	//Data member: vogliamo che le classi derivate ci possano accedere
	double m_massa;
	double m_carica;

};

#endif


//Implementazione di una particella elementare
//In questo caso tutte le proprietà della particella sono costanti, definite nel costruttore di default
#ifndef __Elettrone_h__
#define __Elettrone_h__

class Elettrone: public Particella{

public:
	//Costruttore
	Elettrone();
	//Distruttore
	~Elettrone(){}
	//Metodo
	virtual void Print() const;
};

#endif


//Implementazione di un generico corpo celeste in cui vogliamo poter variare arbitrariamente massa, raggio e che vogliamo identificare tramite un nome
#ifndef __CorpoCeleste_h__
#define __CorpoCeleste_h__

class CorpoCeleste: public Particella{

public:
	//Costruttore
	CorpoCeleste(string nome, double massa, double raggio);
	//Distruttore
	~CorpoCeleste(){}
	//Nuovi metodi
	void SetNome(string nome) {m_nome=nome;}
	void SetMassa(double massa) {m_massa=massa;}
	void SetRaggio(double raggio) {m_raggio=raggio;}
	string GetNome() const {return m_nome;}
	double GetRaggio() const {return m_raggio;}
	virtual void Print() const;
	Posizione GetPosizione() {return m_posizione;}
	double GetPosizioneX() {return m_posizione.GetX();}
      	double GetPosizioneY() {return m_posizione.GetY();}
      	double GetPosizioneZ() {return m_posizione.GetZ();}
	void SetPosizione(double x,double y,double z) {m_posizione.SetX(x);m_posizione.SetY(y);m_posizione.SetZ(z);}
	double PotenzialeGravitazionale(Posizione); //Per usare questo metodo dovrei avere la posizione come se l'universo fosse un s.d.r. cartesiano

protected:
	//Nuovi data membri
	string m_nome;
	double m_raggio;
	Posizione m_posizione;
};

#endif


#ifndef __CampoVettoriale_h__
#define __CampoVettoriale_h__

class CampoVettoriale: public Posizione{

public:
	//Costruttore
	CampoVettoriale(const Posizione&);
	//Distruttore
	~CampoVettoriale(){}
	//Metodi
	double GetFX() const {return m_x;}
	double GetFY() const {return m_y;}
	double GetFZ() const {return m_z;}
	void SetFX(double x) {m_x=x;}
	void SetFY(double y) {m_y=y;}
	void SetFZ(double z) {m_z=z;}
	double Modulo() const; //Calcola la lunghezza del vettore
	void Somma(const CampoVettoriale&); //Modifica il vettore assegnato sommandogli il vettore dato come argomento

protected:
    	double m_x,m_y,m_z,c;

};

#endif


#ifndef __PuntoMateriale_h__
#define __PuntoMateriale_h__

class PuntoMateriale: public Particella, Posizione{

public:
	//Costruttore
    	PuntoMateriale(double massa, double carica, double x, double y, double z);
	//Distruttore
	~PuntoMateriale(){}
	//Metodi
    	CampoVettoriale CampoElettrico(const Posizione&) const;
    	CampoVettoriale CampoGravitazionale(const Posizione&) const;

};

#endif


//Classe astratta per una generica funzione
#ifndef __FunzioneBase_h__
#define __FunzioneBase_h__

class FunzioneBase{

public:
	virtual double Eval(double x) const =0; //Pongo a zero il metodo virtual: è obbligatorio per le classi figlie implementarlo

};

#endif


//funzione f(x)=ax^2+bx+c
#ifndef __Parabola_h__
#define __Parabola_h__

class Parabola: public FunzioneBase{

public:
	//Costruttori
	Parabola();
	Parabola(double a, double b, double c);
	//Distruttore
	~Parabola(){}
	//Metodi
	void SetA(double a) {m_a=a;}
	void SetB(double b) {m_b=b;}
	void SetC(double c) {m_c=c;}
	double GetA() const {return m_a;}
	double GetB() const {return m_b;}
	double GetC() const {return m_c;}
	double Eval(double x) const;

private:
	double m_a, m_b, m_c;

};

#endif


//funzione x=tan(x) (NON RISOLUBILE ANALITICAMENTE)
#ifndef __ENRA_h__
#define __ENRA_h__

class ENRA: public FunzioneBase{

public:
	//Costruttore
	ENRA();
	//Distruttore
	~ENRA(){}
	//Metodi
	double Eval(double x) const;

};

#endif


//funzione f(x)=sin(x)
#ifndef __Seno_h__
#define __Seno_h__

class Seno: public FunzioneBase{

public:
	Seno();
	~Seno(){}
	double Eval(double x) const;

};

#endif


//funzione Gaussiana
#ifndef __Gaussiana_h__
#define __Gaussiana_h__

class Gaussiana: public FunzioneBase{

public:
	Gaussiana(double mu, double sigma);
	~Gaussiana(){}
	double Eval(double x) const;

private:
	double m_mu;
	double m_sigma;

};

#endif


//funzione caratteristica
#ifndef __Caratteristica_h__
#define __Caratteristica_h__

class Caratteristica: public FunzioneBase{

public:
	Caratteristica(double raggio) {m_raggio=raggio;}
	~Caratteristica(){}
	double Eval(double x) const;

private:
	double m_raggio;

};

#endif


//Classe astratta per la ricerca di zeri
#ifndef __Solutore_h__
#define __Solutore_h__

class Solutore{

public:
	void SetPrecisione(double prec) {m_prec = prec;}
	double GetPrecisione() {return m_prec;}
	virtual double CercaZeri(double xmin, double xmax)=0;
	virtual bool Trovato()=0;
	virtual double Incertezza()=0;

protected:
	double m_a, m_b; //Estremi della regione di ricerca
	double m_prec;   //Precisione della soluzione
	const FunzioneBase *m_f;

};

#endif


#ifndef __Bisezione_h__
#define __Bisezione_h__

class Bisezione: public Solutore{

public:
	//Costruttori
	Bisezione();
	Bisezione(double a, double b, double prec, const FunzioneBase* f);
	//Distruttore
	~Bisezione(){}
	//Metodi
	double CercaZeri(double xmin, double xmax);
	bool Trovato();
	double Incertezza();

};

#endif


#ifndef __Vettore_h__
#define __Vettore_h__

class Vettore{

public:
	Vettore();		 //Costruttore di default
	Vettore(unsigned int N); //Costruttore con dimensione del vettore
	~Vettore(){}		 //Distruttore

	unsigned int GetN() const {return m_N;}  //Accede alla dimensione del vettore
	void SetComponent(unsigned int, double); //Modifica la componente i-esima
	double GetComponent(unsigned int) const; //Accede alla componente i-esima
	Vettore(const Vettore&); //Copy constructor
	Vettore& operator=(const Vettore&); //Overloading dell'operatore di assegnazione

protected:
	unsigned int m_N; //Dimensione del vettore
	double* m_v;	  //Vettore di dati

};

#endif


#ifndef __VettoreDati_h__
#define __VettoreDati_h__

class VettoreDati: public Vettore{

public:
	VettoreDati(unsigned int N);
	VettoreDati(unsigned int N, const char* filename);
	~VettoreDati(){}

	void Print();
	void Print(const char* filename);

	void SelectionSort();
	void QuickSort();
	void QuickSort(unsigned int, unsigned int); //Utilizzato per ricorsione
	
	double Media();
	double Var();
	double StdDev();
	double Mediana();
	double CoeffCorr(VettoreDati&);	

};

#endif


#ifndef __VettoreLineare_h__
#define __VettoreLineare_h__

class VettoreLineare: public Vettore{

public:
	//Costruttore
	VettoreLineare(unsigned int N);
	//Distruttore
	~VettoreLineare(){}
	//Algebra vettoriale
	VettoreLineare operator+(const VettoreLineare&) const; //Somma di vettori
	double operator*(const VettoreLineare&) const;	       //Prodotto scalare
	VettoreLineare operator*(const double) const;	       //Prodotto per uno scalare

	VettoreLineare Versore();		 //Versore di un vettore
	double Modulo() const;			 //Modulo di un vettore
	double Cos(const VettoreLineare&) const; //Coseno dell'angolo tra due vettori
	
};

#endif



#ifndef __Random_h__
#define __Random_h__

class Random{

public:
	Random(int seed) {m_seed=seed;} //Costruttore con seme
	~Random(){}

	void SetA(unsigned int a) {m_a=a;} //Metodo per cambiare il valore di a
	void SetC(unsigned int c) {m_c=c;} //Metodo per cambiare il valore di c
	void SetM(unsigned int m) {m_m=m;} //Metodo per cambiare il valore di m

	double Rand01(); //Metodo per generare numeri casuali uniformemente nell'intervallo [0,1]
	double RandEsp(double rate); //Metodo per generare numeri casuali esponenziali tramite il metodo della trasformata
	double RandGauss(double mu, double sigma); //Metodo per generare numeri casuali gaussiani tramite il metodo Box-Muller

private:
	unsigned int m_a;
	unsigned int m_c;
	unsigned int m_m;
	unsigned int m_seed;

};

#endif


#ifndef __AcceptReject_h__
#define __AcceptReject_h__

class AcceptReject{

public:
	AcceptReject(FunzioneBase *f, double a, double b, double max);
	~AcceptReject();

	double RandGauss(); //Metodo per generare numeri casuali gaussiani tramite il metodo Accept-Reject
	
private:
	FunzioneBase *m_f;
	double m_a;
	double m_b;
	double m_max;
	Random *s;
	Random *t;

};

#endif


#ifndef __Integral_h__
#define __Integral_h__

class Integral{

public:
	Integral(double a, double b, FunzioneBase *f);
	Integral(double a, double b, double prec, FunzioneBase *f); //Costruttore per metodo dei trapezi
	~Integral(){}


	//Metodi per l'integrazione:

	double Midpoint(int nstep);
	double Simpson(int nstep);
	double Trapezi(); //Il metodo dei trapezi è a precisione fissata


	//Metodi Monte Carlo per l'integrazione:

	//MONODIMENSIONALE
	double Media(int n);
	double HitOrMiss(int n, double max);
	double GetErrorMedia() const {return m_errorMedia;}
	double GetErrorHitOrMiss() const {return m_errorHitOrMiss;}
	double NMedia(double precision) {return (int)(pow(m_kMedia/precision,2))+1;}
	double NHitOrMiss(double precision) {return (int)(pow(m_kHitOrMiss/precision,2))+1;}

	//MULTIDIMENSIONALE
	double SetDim(int dim) {m_dim=dim;}
	double HitOrMiss(int n);

private:
	double m_a, m_b;
	double m_sum;
	double m_h;
	int m_sign;
	double m_integral;
	double m_prec; //Precisione fissata per metodo dei trapezi
	FunzioneBase *m_integrand;

	//Data membri metodi Monte Carlo:
	Random m_generatore;

	//MONODIMENSIONALE
	double m_errorMedia;
	double m_errorHitOrMiss;
	double m_kMedia;
	double m_kHitOrMiss;

	//MULTIDIMENSIONALE
	int m_dim;
	double *v;

};

#endif


//Classe astratta, restituisce la derivata da valutare nel punto r al tempo t
#ifndef __FunzioneVettorialeBase_h__
#define __FunzioneVettorialeBase_h__

class FunzioneVettorialeBase{

public:
	virtual VettoreLineare Eval(double t, const VettoreLineare& x) const =0;

};

#endif


//Classe astratta per un integratore di equazioni differenziali
#ifndef __EquazioneDifferenzialeBase_h__
#define __EquazioneDifferenzialeBase_h__

class EquazioneDifferenzialeBase{

public:
	virtual VettoreLineare Passo(double t, const VettoreLineare& inizio, double h, FunzioneVettorialeBase *f) const =0;

};

#endif


#ifndef __OscillatoreArmonico_h__
#define __OscillatoreArmonico_h__

class OscillatoreArmonico: public FunzioneVettorialeBase{

public:
	OscillatoreArmonico(double omega0);
	~OscillatoreArmonico(){}
	VettoreLineare Eval(double t, const VettoreLineare& x) const;

	double GetOmega() const {return m_omega;}
	void SetOmega(double omega) {m_omega=omega;}

private:
	double m_omega; //Pulsazione propria

};

#endif


#ifndef __OscillatoreForzato_h__
#define __OscillatoreForzato_h__

class OscillatoreForzato: public FunzioneVettorialeBase{

public:
	OscillatoreForzato(double omega0, double omega, double alpha);
	~OscillatoreForzato(){}
	VettoreLineare Eval(double t, const VettoreLineare& x) const;

	double GetOmega0() const {return m_omega0;}
	double GetOmega() const {return m_omega;}
	double GetAlpha() const {return m_alpha;}
	void SetOmega0(double omega0) {m_omega0=omega0;}
	void SetOmega(double omega) {m_omega=omega;}
	void SetAlpha(double alpha) {m_alpha=alpha;}

private:
	double m_omega0; //Pulsazione iniziale
	double m_omega;  //Pulsazione
	double m_alpha;  //Alpha

};

#endif


#ifndef __Pendolo_h__
#define __Pendolo_h__

class Pendolo: public FunzioneVettorialeBase{

public:
	Pendolo(double l);
	~Pendolo(){}
	VettoreLineare Eval(double t, const VettoreLineare& x) const;

	double GetL() const {return m_l;}
	void SetL(double l) {m_l=l;}

private:
	double m_l; //Lunghezza del filo

};

#endif


#ifndef __Sfera_h__
#define __Sfera_h__

class Sfera: public FunzioneVettorialeBase{

public:
	Sfera(double eta, double rho, double rho0, double R);
	~Sfera(){}
	VettoreLineare Eval(double t, const VettoreLineare& x) const;

	double GetEta() const {return m_eta;}
	double GetRho() const {return m_rho;}
	double GetRho0() const {return m_rho0;}
	double GetR() const {return m_R;}
	void SetEta(double eta) {m_eta=eta;}
	void SetRho(double rho) {m_rho=rho;}
	void SetRho0(double rho0) {m_rho0=rho0;}
	void SetR(double R) {m_R=R;}

private:
	double m_eta;
	double m_rho;
	double m_rho0;
	double m_R;

};

#endif


#ifndef __Eulero_h__
#define __Eulero_h__

class Eulero: public EquazioneDifferenzialeBase{

public:
	Eulero(){}
	~Eulero(){}
	VettoreLineare Passo(double t, const VettoreLineare& x, double h, FunzioneVettorialeBase *f) const;

};

#endif


#ifndef __RungeKutta_h__
#define __RungeKutta_h__

class RungeKutta: public EquazioneDifferenzialeBase{

public:
	RungeKutta(){}
	~RungeKutta(){}
	VettoreLineare Passo(double t, const VettoreLineare& x, double h, FunzioneVettorialeBase *f) const;

};

#endif


#ifndef __EsperimentoPrisma_h__
#define __EsperimentoPrisma_h__

class EsperimentoPrisma{

public:
	EsperimentoPrisma();
	~EsperimentoPrisma(){};
	void Esegui();
	void Analizza();

	double Getlambda1() const {return m_lambda1;}
	double Getlambda2() const {return m_lambda2;}
	double Getalpha() const {return m_alpha;}
	double Getsigmat() const {return m_sigmat;}
	double GetAInput() const {return m_A_input;}
	double GetAMisurato() const {return m_A_misurato;}
	double GetBInput() const {return m_B_input;}
	double GetBMisurato() const {return m_B_misurato;}
	double Getn1Input() const {return m_n1_input;}
	double Getn1Misurato() const {return m_n1_misurato;}
	double Getn2Input() const {return m_n2_input;}
	double Getn2Misurato() const {return m_n2_misurato;}
	double Gett0Input() const {return m_t0_input;}
	double Gett0Misurato() const {return m_t0_misurato;}
	double Gett1Input() const {return m_t1_input;}
	double Gett1Misurato() const {return m_t1_misurato;}
	double Gett2Input() const {return m_t2_input;}
	double Gett2Misurato() const {return m_t2_misurato;}
	double Getdm1Input() const {return m_dm1_input;}
	double Getdm1Misurato() const {return m_dm1_misurato;}
	double Getdm2Input() const {return m_dm2_input;}
	double Getdm2Misurato() const {return m_dm2_misurato;}

private:
	//Generatore di numeri casuali da usare in Esegui()
	Random m_generatore;
	//Parametri dell'apparato sperimentale
	double m_lambda1, m_lambda2, m_alpha, m_sigmat;
	//Valori delle quantità misurabili:
	// m_input:    valore assunto come ipotesi nella simulazione
	// m_misurato: valore misurato, calcolato in Esegui() ed Analizza()
	double m_A_input   , m_A_misurato ;
	double m_B_input   , m_B_misurato ;
	double m_n1_input  , m_n1_misurato;
	double m_n2_input  , m_n2_misurato;
	double m_t0_input  , m_t0_misurato;
	double m_t1_input  , m_t1_misurato;
	double m_t2_input  , m_t2_misurato;
	double m_dm1_input , m_dm1_misurato;
	double m_dm2_input , m_dm2_misurato;

};

#endif

/*************************************************STRUCT*********************************************************************************************/

struct vettore{

	unsigned int N; //Dimensione
	double* v;	//Puntatore ai dati
};

/*************************************************FUNZIONI*******************************************************************************************/

//////////////////////////////////////////////////LEZIONE_1///////////////////////////////////////////////////////////////////////////////////////////

int numeroScambi();

void scambiaByValue(double,double); //NON FUNZIONA

void scambiaByRef(double&,double&);

void scambiaByPointer(double*,double*);


double media(double [],int);

double var(double [],int);


double modulo(double); //Modulo numero

double modulo(double,double); //Modulo vettore 2 componenti

double modulo(double [],int); //Modulo vettore n componenti

//////////////////////////////////////////////////LEZIONE_2///////////////////////////////////////////////////////////////////////////////////////////

struct vettore read(unsigned int N,const char* filename); //Crea una struttura vettore leggendo N numeri dal file di nome filename

void print(const struct vettore,const char* filename); //Funzione per stampare il contenuto di un vettore

void selection_sort(struct vettore); //Funzione per ordinare gli elementi di un vettore

//////////////////////////////////////////////////LEZIONE_5///////////////////////////////////////////////////////////////////////////////////////////

int sign(double x); //Funzione segno
