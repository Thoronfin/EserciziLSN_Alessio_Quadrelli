#ifndef _Rw_h_
#define _Rw_h_


#include "random.h"


class Rw: public Random {

	private:
	
	//coordinate x,y,z della posizione nel random walk
	double x;
	double y;
	double z;
	double a; //lunghezza del passo
	Random rnd; //serve per la generazione dei numeri casuali 
	

	public:

	Rw(double passo); //costruttore
	~Rw(); //distruttore
	void Restart();	//inizializza a 0 la posizione
	void Passo(); //aggiorna le posizioni dopo un singolo passo in un reticolo cubico tridimensionale
	double ReticoloCubico(int N); //permette di effettuare N passi nel reticolo cubico
	void PassoIsotropo(); //aggiorna le posizioni dopo un singolo passo in direzione casuale nello spazio tridimensionale
	double Diffusione(int N); //permette di effettuare N passi nello spazio tridimensionale
		
	
};


#endif
