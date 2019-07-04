#ifndef _OptionPrice_h_
#define _OptionPrice_h_


#include "random.h"


class OptionPrice: public Random {

	private:
	
	
	double _S0; //S(0)
	double _K; //strike price
	double _r; //tasso di interesse
	double _sigma; //volatilità
	double _ST; //S(T)
	Random _rnd; //serve per la generazione dei numeri casuali 
	

	public:

	OptionPrice(double S0, double K, double r, double sigma); //costruttore
	~OptionPrice(); //distruttore
	void GBM(double inizio, double fine); //implementa il moto browinano geometrico S(fine) applicato a S(inizio)
	double Put(double inizio, double fine, int passi); //valutazione della put option da t=inzio a t=fine in un numero di  passi "passi"
	double Call(double inizio, double fine, int passi); //valutazione della call option da t=inzio a t=fine in un numero di  passi "passi"
	void Restart(); // inizializza S(t) a S0 per iniziare nuovamente il calcolo

		
	
};


#endif
