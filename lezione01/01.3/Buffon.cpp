#include "Buffon.h"
#include "random.h"
#include <math.h>
#include <cmath>
#include <iostream>
#include <fstream>

using namespace std;

Buffon::Buffon(): Random() {

	d=0.;
	L=0.;
	throws=0;
}

Buffon::Buffon(double a,double b,int N): Random() {

	d=a;
	L=b;
	throws=N;

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

Buffon:: ~Buffon(){
	rnd.SaveSeed();
}

void Buffon:: SetThrows(int N){

	throws=N;

}

double Buffon:: pi() {
	 

	int hit=0;
	double yc;
	double sin_theta;
	double y1,y2;
	

	for(int i=0; i<throws; i++){

		//genero la coordinata y del centro dell'ago essendo il problema invariante per traslazioni in x 
		//e seno di un angolo distribuito uniformemente tra 0 e 2*pi  per poi ricavare le coordinate degli estermi dell'ago
		yc = rnd.Rannyu(-(d/2),(d/2));
		sin_theta= rnd.Sin();
		// ricavo le coordinate y degli estremi dell'ago
		y1= yc + (L/2)*sin_theta;
		y2= yc - (L/2)*sin_theta;

		//controllo se l'ago attraversa la linea retta ed eventualmente aggiorno il contatore dei successi
		if ( (y1 >=0 and y2<=0) or (y1 <=0 and y2>=0)){
			
			hit ++;
		}	
	}

	//restituzione del valore del pigreco
	return (2.*L*float(throws))/(d*float(hit));
}
