#ifndef _Integrator_h_
#define _Integrator_h_

#include "Funzioni.h"
#include "random.h"

//classe per il calcolo degli integrali montecarlo
class Integral: public Random{
	
	private:
	
	int _sign; //serve per tenere conto di eventuali scambi negli estremi di integrazione
	double _a, _b; // estremi dell'integrale
	double _sum; //somme parziali
	double _integral; //risultato integrale
	Random rnd; // per la generazione dei numeri casuali
	FunzioneBase * _integrand; //funzione da integrare



	public:

	Integral(double a, double b, FunzioneBase * function); //costruttore per integrare function tra a e b
	~Integral(); //distruttore
	double Media(int N); //calcola l'integrale montecarlo generando numeri uniformenete distribuiti tra _a e _b
	double MediaRetta(int N); //calcola l'integrale montecarlo generando numeri distribuiti tra 0 e 1 con p(x) = 2*(1-x)
	double GetResult() const; // restituisce il valore dell'integrale calcolato

};



#endif

