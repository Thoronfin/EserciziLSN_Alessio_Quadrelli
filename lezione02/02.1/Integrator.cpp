#include "Integrator.h"
#include "random.h"
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;


Integral::Integral(double a, double b, FunzioneBase * function): Random() {

	_integrand = function;
	//ordino gli estremi in ordine crescente e se non dovessero già essere ordinati assegno a _sign -1 per tenerne conto
	_a = min(a,b);
	_b = max(a,b);
	if ( _a>= _b){
		_sign=-1;
	}
	else _sign=1;

	// perparo il generatore di numeri casuali
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
	  	 rnd.SetRandom(seed,p1,p2);
	 	}
	     }
	    input.close();
	  } else cerr << "PROBLEM: Unable to open seed.in" << endl;


}

//distruttore
Integral::~Integral() {
		rnd.SaveSeed();
}

//Calcola l'integrale montecarlo usando numeri distribuiti uniformemente tra [_a,_b)
double Integral::Media(int N){

	_sum=0;
	double intervallo= _b-_a;

	for (int i=1; i <= N; i++){
		double x=rnd.Rannyu(_a,_b);
		double f= _integrand -> Eval(x);
		_sum += f;
	}	

	_integral = _sign*(_sum/double(N))*intervallo;

	return _integral;

};

//Calcola l'integrale montecarlo usando numeri distribuiti con P(x)=2*(1-x) tra [0,1)
double Integral::MediaRetta(int N){

	_sum=0;
	double intervallo= _b-_a;

	for (int i=1; i <= N; i++){
		double x=rnd.Retta();
		double f= _integrand -> Eval(x);
		_sum += f;
	}	

	_integral = _sign*(_sum/double(N))*intervallo;

	return _integral;

};

//restituisce il risultato dell'integrale senza bisogno di calcolarlo nuovamente
double Integral::GetResult() const{
	return _integral;
}
