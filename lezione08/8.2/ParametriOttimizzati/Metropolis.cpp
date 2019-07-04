#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Metropolis.h"
#include "Statistica.h"

using namespace std;

Metropolis::Metropolis(FunzioneBase *f, vector<double> posizione_0, double passo, int M, int bins) : _next(posizione_0.size()), _psi(bins, vector<double>(M)) {

	//Metropolis
	_dim=posizione_0.size();
	_accepted=0;
	_step_fatti=0;
	_n_step=0;
	_passo=passo;
	_posizioni.push_back( posizione_0);
	_f=f;

	//dati per istogramma
	_estremo=3; //estremi dell'intervallo
	_nbins=bins; //bins
	_niterazioni=M; //numero di istogrammi

	//inizializzazione matrice per realizzazione degli istogrammi
	for(int j=0; j<_niterazioni; j++){
		for(int i=0; i<_nbins;i++){
			_psi[i][j]=0;
		}
	}


	
	

	//inizializzazione generatore numeri casuali
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	   
	  
	if (Primes.is_open()){
	   Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
	   while ( !input.eof() ){
	      input >> property;
	       if( property == "RANDOMSEED" ){
		 input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
		 _rnd.SetRandom(seed,p1,p2);
	      }
	   }
	   input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

	

}


Metropolis::~Metropolis(){

	_rnd.SaveSeed();
}

void Metropolis::Next(){ //passo proposto

	for(int i=0; i<_dim; i++){
		_next[i] = _rnd.Rannyu(- _passo,_passo) + _posizioni[_step_fatti][i];
	}

}

bool Metropolis::Accept(){ //controllo accettazione del passo

	 //modulo quadro della funzione d'onda, poichè eval restituisce solo il valore di essa
	double p = pow(_f-> Eval(_posizioni[_step_fatti]),2);
	double p_new =pow( _f-> Eval(_next),2);
	double accept= min(1., p_new/p);
	

	if(accept>=1){
		
		_accepted++;
		return true;
	}
	else{
		double a=_rnd.Rannyu();
		if( a<= accept){
			
			_accepted++;
			return true;
		}
		return false;
	}
			

}


void Metropolis::Step(){ // algoritmo di Metropolis

	Next(); //nuova configurazione proposta
	bool accept= Accept(); //verifica dell'accetazione della configurazione

	if(accept == true){
		_posizioni.push_back(_next);
			
	}
	else{
		
		_posizioni.push_back(_posizioni[_step_fatti]);
		
	}

	
	_step_fatti ++;

}

void Metropolis::Sample(int n_step){
	
	_n_step=n_step;
	for(int i=0; i<_n_step; i++){
		Step();
	}

}

void Metropolis::Restart(){

	//vengono concellati tutte le posizioni di prova esclusa l'ultima che serve come inizio per il campionamento;
	_posizioni.erase(_posizioni.begin(), _posizioni.end()-1);


	//vengono posti a zero i contatori
	_accepted=0;
	_step_fatti=0;
}


void Metropolis::Initialization(int step_inizializzazione){ //metodo per equilibrare il Metropolis e per inizializzare la configurazione iniziale con l'ultima campionata
	
	for(int i=0; i<step_inizializzazione; i++){
		Step();	
	}

	Restart();
	
}


void Metropolis::AnalisiHistogram(int N, const char* filename){ //data blocking sugli istogrammi
		
	double dr=(2*_estremo)/(double)_nbins;
	AnalisiMatrici(N, _psi , dr, _estremo, filename);

}

void Metropolis::Histogram(int j){
	
	double dr=(2*_estremo)/(double)_nbins; //ampiezza del bin
	int index;
	double dato;

	//riempimento istogramma
	for(int i=0; i<_step_fatti; i++){
		dato = _posizioni[i][0]; //lettura dei dati

		if( abs(dato) < _estremo){ // si selezionano i dati solo nell'intervallo desiderato
			
			index= int((_estremo +dato)/dr); // calcolo indice del bin corrispondente al dato traslando di estrmo in modo che gli 									inidici partano da zero
			_psi[index][j] +=1; //incremento del numero di dati nel bin

		}
	}


	//fattore di normalizzazione
	double somma=0.;
	for(int i=0; i< _psi.size(); i++){
		somma += _psi[i][j];	
	}
	double norm= 1./(dr*somma);

	//normalizzazione istogramma
	for(int i=0; i<_psi.size(); i++){
		_psi[i][j] = _psi[i][j]*norm;
		
	}
	
	

}


double Metropolis::Integration(FunzioneBase* hamiltonian){ //metodo per il calcolo dell'integrale Monte Carlo dell'enegia

	
	double integrale=0;

	for(int i=0; i<_step_fatti; i++){

		integrale += hamiltonian->Eval(_posizioni[i]);
	}
 
	return integrale/_step_fatti;
}

double Metropolis::AcceptanceProbability(){ //calcolo della probabilità di accettazione delle configurazioni
	return (_accepted*1.)/_n_step;
}
