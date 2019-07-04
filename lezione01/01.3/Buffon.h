#ifndef _Buffon_h_
#define _Buffon_h_

#include "random.h"

class Buffon: public Random{

	private:

		double d; // distanza fra le linee orizzontali
		double L; //lunghezza ago
		int throws; // numero di lanci dell'ago 
		Random rnd; //generatore di numeri casuali

	public:

		Buffon(); //costruttore vuoto, inizializza i daa membri a zero
		Buffon(double, double, int); //costruttore con cui assegnare d,L, throws
		~Buffon(); //distruttore
		void SetThrows(int); // permette di variare il numero di lanci dell'ago
		double pi(); //restituisce il valore del pigreco calcolato con con un totale di throws lanci

};

#endif
