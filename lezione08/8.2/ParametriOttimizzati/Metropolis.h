#ifndef _Metropolis_h_
#define _Metropolis_h_

#include "random.h"
#include "Funzioni.h"
#include <vector>

using namespace std;

class Metropolis {

private:
	Random _rnd;
	FunzioneBase *_f;
	double _passo; // ampiezza del passo casuale
	double _dim;//dimensione del rpoblema
	int _n_step; 
	int _step_fatti; //totale delle configurazioni campionate
	int _accepted; //configurazioni accettate
	vector<double> _next; //nuova configurazione proposta
	vector < vector<double> > _posizioni; //congigurazioni campionate

	double _estremo; //estrmi per il calcolo istogramma
	int _nbins; // numero di bin nell'istogramma
	vector<vector<double>> _psi; // serve per istogramma delle psi
	int _niterazioni; // numero di istogrammi su cui mediare

protected:


public:
	
  // constructors
	Metropolis(FunzioneBase *f, vector<double> posizione_0, double passo, int M, int _nbins);
  // destructor
	~Metropolis();
  // methods
	void Next(); // genera il punto successivo
	bool Accept(); // controllo se la nuova posizione viene accettata
	void Step(); // aggiorna la posizione
	void Sample(int); //permette di fare N passi
	void Restart(); //rimette i contatori a zero e permette di ripartire dall'ultima posizione raggiunta
	void Initialization(int); // permette di effettuare un numero di passi desiderato prima di iniziare il vero e proprio campionamento
	double AcceptanceProbability(); //calcolo della probabilità che il passo proposto vega accettato
	double Integration(FunzioneBase* hamiltonian); // calcolo integrale campionando con metropolis
	void Histogram(int j); // aggiornamento istogramma della funzione d'onda di trial
	void AnalisiHistogram(int N, const char*);// Data blocking su istogramma di psi
	

};

#endif
