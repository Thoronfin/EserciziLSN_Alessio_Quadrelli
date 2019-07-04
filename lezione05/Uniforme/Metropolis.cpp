#include <cmath>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include "Metropolis.h"

using namespace std;

Metropolis::Metropolis(FunzioneBase *f, vector<double> posizione_0, double passo, int dim) : _next(dim) {

	_accepted=0;
	_f=f;
	_step_fatti=0;
	_n_step=0;
	_dim=dim;
	_passo=passo;
	_posizioni.push_back( posizione_0);
	

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

void Metropolis::Next(){ //nuova configurazione proposta

	for(int i=0; i<_dim; i++){
		_next[i] = _rnd.Rannyu(- _passo,_passo) + _posizioni[_step_fatti][i];
	}

}

bool Metropolis::Accept(){ //controllo accettazione della configurazione

	double p = _f-> Eval(_posizioni[_step_fatti]);
	double p_new = _f-> Eval(_next);
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

	Next(); //nuova configurazione
	bool accept= Accept(); //verifica se accettare o no la nuova configurazione

	
	if(accept == true){ //nuova configurazione accettata
		_posizioni.push_back(_next);
			
	}
	else{ //nuova configurazione rifiutata
		
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


void Metropolis::Initialization(int step_inizializzazione){ //campionamenti per equilibrazione
	
	for(int i=0; i<step_inizializzazione; i++){
		Step();	
	}

	Restart();
	
}


void Metropolis::Print(const char* filename){ //alvataggio su file del campionamento

	ofstream Risultati;
	Risultati.open(filename);

	if (Risultati.is_open()){
		for(int j=0; j<_n_step+1; j++){
			for(int i=0; i <_dim; i++){
				Risultati <<  _posizioni[j][i] << " " ;
			}
		Risultati << endl;
		}
	     	 
   	} else cerr << "PROBLEM: Unable to save the results" << endl;
  	Risultati.close();
	
}


vector<double> Metropolis::Raggi(){ //salvataggio delle configurazioni

	vector<double> raggi(_n_step,0);

	for(int j=0; j<_n_step; j++){
			for(int i=0; i <_dim; i++){
				raggi[j] += pow( _posizioni[j][i], 2);
			}
		raggi[j]= sqrt(raggi[j]);
		}

return raggi;
}

double Metropolis::AcceptanceProbability(){ //probabilità di accettazione
	return (_accepted*1.)/_n_step;
}
