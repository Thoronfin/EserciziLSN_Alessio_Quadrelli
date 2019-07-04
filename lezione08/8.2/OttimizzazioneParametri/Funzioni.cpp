#include <cmath>
#include<vector>
#include <iostream>
#include "Funzioni.h"

using namespace std;



Psi:: Psi(){ //costruttore di default della funzione d'onda

	_mu=0.;
	_sigma=0.;

}


Psi:: Psi(double mu, double sigma){ //costruttore della funzione d'onda con parametri assegnati

	_mu=mu;
	_sigma=sigma;

}

Psi::~Psi(){

}

// distribuzioni da campionare a meno della normalizzazione 
double Psi:: Eval(vector<double> x) const {

	double dx= pow((x[0]-_mu),2)/(2.*pow(_sigma,2)); // gaussiana centrata nella buca di destra
	double sx= pow((x[0]+_mu),2)/(2.*pow(_sigma,2)); // gaussiana centrata nella buca di sinistra
	
	return exp(-sx)+exp(-dx); //somma delle due gaussiane
	

}


double Hamiltoniana:: Eval(vector<double> x) const { //(H applicata a psi)/psi

	double _sigma= _psi->GetSigma();
	double dx= pow(x[0]-_psi->GetMu(),2);
	double sx= pow(x[0]+_psi->GetMu(),2);
	double exp_dx= exp(-dx/(2.*pow(_sigma,2))); // gaussiana centrata nella buca di destra
	double exp_sx= exp(-sx/(2.*pow(_sigma,2))); // gaussiana centrata nella buca di sinistra
	double fac= 1./pow(_sigma,2);
	double fac_2= pow(fac,2);
	//enegia cinetica ottenuta dal laplaciano di psi
	double Kin= (-0.5*  ( (-exp_sx*fac) + (-exp_dx*fac) + (exp_sx*fac_2* sx) + (exp_dx*fac_2* dx) ) )/_psi -> Eval(x);
	// energia potenziale: il potenziale agisce moltiplicativamente e nel rapporto le funzioni d'onda si cancellano
	double V= ( pow(x[0],4) - 2.5*pow(x[0],2) );
	return  Kin+ V;
	

}



