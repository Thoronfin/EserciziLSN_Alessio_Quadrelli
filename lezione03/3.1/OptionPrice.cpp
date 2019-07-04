#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "OptionPrice.h"
using namespace std;

OptionPrice::OptionPrice(double S0, double K, double r, double sigma): Random(){

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
	

	_S0=S0;
	_K=K;
	_r=r;
	_sigma=sigma;
	_ST=S0;

}

OptionPrice::~OptionPrice(){
	_rnd.SaveSeed();
}


void OptionPrice::GBM(double inizio, double fine){ //moto browniano geometrico
	
	_ST = _ST*exp( (_r- 0.5*pow(_sigma,2))*(fine-inizio) + _sigma*_rnd.Gauss(0.,1.)*sqrt(fine-inizio)   );

}

void OptionPrice::Restart(){

	_ST=_S0;
}

double OptionPrice::Put(double inizio, double fine, int passi){

	//moto browniano geometrico che può essere discretizzato in "passi" passi
	double incremento = (fine-inizio)/passi;
	for(double i=0; i< passi; i++){
		GBM(inizio + i*incremento, inizio + (i+1)*incremento); 	
	}

	//valutazione del profitto tenendo
	double a = exp(-_r*(fine-inizio))*max(0., _K -_ST);
	Restart();
	
	return a;

	

}

double OptionPrice::Call(double inizio, double fine, int passi){

	//moto browniano geometrico che può essere discretizzato in "passi" passi
	double incremento = (fine-inizio)/passi;
	for(double i=0; i< passi; i++){
		GBM(inizio + i*incremento, inizio + (i+1)*incremento); 	
	}

	//valutazione del profitto
	double a = exp(-_r*(fine-inizio))*max(0., _ST - _K);
	Restart();
	
	return a;

}


