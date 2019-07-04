#ifndef _Funzioni_h_
#define _Funzioni_h_
#include <cmath>

//classe astratta contente le funzioni da integrare
class FunzioneBase {
	public:
		virtual double Eval(double x) const=0;

};


// Funzione da integrare in esercizio 02.1
class Coseno: public FunzioneBase {

	virtual double Eval (double x) const;

};

// G è l'integranda dell'esercizio 02.1 modificata
class G :public FunzioneBase{

	virtual double Eval (double x) const;

};





#endif //_Funzioni_h_

