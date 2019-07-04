#include <cmath>
#include <fstream>
#include <string>
#include <iostream>
#include "Rw.h"
using namespace std;

Rw::Rw(double passo):Random(){

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
            	rnd.SetRandom(seed,p1,p2);
        }
        }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	//posizione iniziale è l'origine
	x=0;
	y=0;
	z=0;
	//inizializzazione lunghezza del passo
	a=passo;
}

Rw::~Rw(){
	rnd.SaveSeed();
}

void Rw::Restart(){
	x=0;
	y=0;
	z=0;
}

//generazione di un numero casuale da 0 a 5 che determina verso quale primo vicino avviene il passo
void Rw::Passo(){

	int r = rnd.Rannyu(0,6);

	if(r==0){x++;}
	if(r==1){y++;}
	if(r==2){z++;}
	if(r==3){x--;}
	if(r==4){y--;}
	if(r==5){z--;}
	
}

//vengono effettuati N passi nel reticolo cubico, poi la posizione viene inizializzatra nuovamente a 0
double Rw::ReticoloCubico(int N){

	for(int i=0;i<N;i++){
		Passo();
	}


	double distanza = a*sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	Restart();
	return distanza;
}

//viene effettuato un passo in una direzione casuale nello spazio tridimensionale
void Rw::PassoIsotropo(){

	double theta = rnd.Rannyu(0,M_PI);
	double phi = rnd.Rannyu(0,2*M_PI);

	x+=sin(theta)*cos(phi);
	y+=sin(theta)*sin(phi);
	z+=cos(phi);
       
}

//vengono effettuati N passi nello spazio tridimensionale, poi la posizione viene inizializzatra nuovamente a 0
double Rw::Diffusione(int N){
	for(int i=0;i<N;i++){
		PassoIsotropo();
	}

	double distanza = a*sqrt(pow(x,2)+pow(y,2)+pow(z,2));
	Restart();

	return distanza;
}
