#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <ctime>

/*************************************************CLASSI*********************************************************************************************/

#ifndef __posizione_h__
#define __posizione_h__

class Posizione{

public:
	//Costruttori
	Posizione();
	Posizione(double x, double y, double z);
	//Distruttore
	~Posizione ();
	//Metodi
	double X() const; //Coordinate cartesiane
	double Y() const;
	double Z() const;
	double R() const; //Coordinate sferiche
	double Phi() const;
	double Theta() const;
	double Rho() const; //Raggio delle coordinate cilindriche
	double Distanza(const Posizione&) const; //Distanza da un altro punto

private:
	double m_R, m_Phi, m_Theta;

};

#endif

/*************************************************STRUCT*********************************************************************************************/

struct Vettore{

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

struct Vettore read(unsigned int N,const char* filename); //Crea una struttura Vettore leggendo N numeri dal file di nome filename

void print(const struct Vettore,const char* filename); //Funzione per stampare il contenuto di un vettore

void selection_sort(struct Vettore); //Funzione per ordinare gli elementi di un vettore


