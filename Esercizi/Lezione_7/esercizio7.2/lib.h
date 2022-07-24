#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <climits>

#define G 6.67E-11 // N*m^2/Kg^2
#define epsilon_0 8.85418781762E-12 // F/m

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
	~Posizione ();
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
	~Particella();
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
	~Elettrone();
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
	~CorpoCeleste();
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
	~CampoVettoriale();
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
	~PuntoMateriale();
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
	~Parabola();
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


//funzione f(x)=sin(x)
#ifndef __Seno_h__
#define __Seno_h__

class Seno: public FunzioneBase{

public:
	Seno();
	~Seno();
	double Eval(double x) const;

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
	~ENRA();
	//Metodi
	double Eval(double x) const;

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
	~Bisezione();
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
	~Vettore();		 //Distruttore

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
	VettoreDati(unsigned int N,const char* filename);
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


#ifndef __Integral_h__
#define __Integral_h__

class Integral{

public:
	Integral(double a, double b, FunzioneBase *f);
	Integral(double a, double b, double prec, FunzioneBase *f); //Costruttore per metodo dei trapezi
	double Midpoint(int nstep);
	double Simpson(int nstep);
	double Trapezi(); //Il metodo dei trapezi è a precisione fissata

private:
	double m_a, m_b;
	double m_sum;
	double m_h;
	int m_sign;
	double m_integral;
	double m_prec; //Precisione fissata per metodo dei trapezi
	FunzioneBase *m_integrand;

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
