#include <cmath>
#include<vector>
#include <iostream>
#include "Funzioni.h"

using namespace std;

#define a0 1

// modulo quadro della funzione d'onda dell'orbitale 1s dell'atomo di idrogeno
double uno_s:: Eval(vector<double> x) const {

	if( x.size() !=3){
		cout << "errore nella funzion 1s, la dimensione del vettore poszione è sbagliata" << endl;
		return 0.;
	}
	else{
	double r= sqrt(pow(x[0],2) + pow(x[1],2) + pow(x[2],2) ); // dalle coordinate cartesiane ricavo quelle sferiche
	return pow( (pow(a0, -3./2.)/sqrt(M_PI))*exp(-r/a0) ,2);
	}

}

// modulo quadro della funzione d'onda 2p dell'orbitale dell'atomo di idrogeno
double due_p:: Eval(vector<double> x) const {

	if( x.size() !=3){
		cout << "errore nella funzion 1s, la dimensione del vettore poszione è sbagliata" << endl;
		return 0;
	}
	else{
	double r= sqrt(pow(x[0],2) + pow(x[1],2) + pow(x[2],2) );  // dalle coordinate cartesiane ricavo quelle sferiche
	return pow((pow(a0, -5./2.) / 8.)*sqrt(2/M_PI)*x[2]*exp(-r/(2.*a0)),2); // z=r*cos(theta)
	}
}
