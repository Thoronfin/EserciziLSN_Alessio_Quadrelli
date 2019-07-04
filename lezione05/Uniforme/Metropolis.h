#ifndef _Metropolis_h_
#define _Metropolis_h_

#include "random.h"
#include "Funzioni.h"
#include <vector>

using namespace std;

class Metropolis {

private:
	Random _rnd;
	FunzioneBase *_f; //distribuzione da campionare
	double _passo; // ampiezza del passo casuale
	int _n_step;
	int _step_fatti; //totale delle configurazioni generate
	int _accepted; //configurazioni proposte ed accettate
	int _dim; //dimensione del problema
	vector<double> _next;
	vector < vector<double> > _posizioni;

protected:


public:
	
  // costruttore
	Metropolis(FunzioneBase *f, vector<double> posizione_0, double passo, int dim);
  // distruttore 
	~Metropolis();
  // metodi
	void Next(); // genera il punto successivo
	bool Accept(); // controllo se la nuova posizione viene accettata
	void Step(); // aggiorna la posizione
	void Sample(int); //permette di fare N passi
	void Restart(); //rimette i contatori a zero e permette di ripartire dall'ultima posizione raggiunta
	void Initialization(int); // permette di effettuare un numero di passi desiderato prima di iniziare il vero e proprio campionamento
	void Print(const char*);// Salva su file i campionamenti effettuati
	vector<double> Raggi(); //calcolo dell distanza dall'origine
	double AcceptanceProbability(); //calcolo della probabilità che il passo proposto vega accettato

};

#endif
