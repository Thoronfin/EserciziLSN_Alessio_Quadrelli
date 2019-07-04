#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "ChiQuadrato.h"

using namespace std;

void SetRandomGenerator(Random* rnd); //serve per inizializzare il generatore di numeri casuali

int main (int argc, char *argv[]){

	//prepaazione per la generazione dei numeri casuali
	Random rnd;
	SetRandomGenerator(&rnd);


	int N_tot=1000000;  //numero totale di punti per il chi quadrato       
   	int M=100;     //numero di sotto-intervalli in cui è diviso l'intervallo [0,1)          


	//numeri casuali utilizzati per il calcolo dei chi quadrato
	vector<double> s(N_tot);
  	for(int j=0; j<N_tot; j++){
		s.at(j) = rnd.Rannyu();
  	}
  

	//Analizzo i dati ottenuti facendo media e deviazione standard della media 
	//sommando progressivamente le N medie ottenute da int(M/N) risultati dell'integrale
	//Poi i dati vengono salvati su file
	ChiQuadrato(M,N_tot,s,"../risultati01.1.3.out");

	rnd.SaveSeed(); 

return 0;
}

//funzione per l'inizializzazione del generatore di numeri casuali
void SetRandomGenerator(Random* rnd){

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
		 rnd->SetRandom(seed,p1,p2);
	      }
	   }
	   input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;


}
