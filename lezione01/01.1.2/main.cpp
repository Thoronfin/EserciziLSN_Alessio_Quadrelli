/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include "random.h"
#include "Statistica.h"
#include <cmath>

using namespace std;

void SetRandomGenerator(Random* rnd);

 
int main (int argc, char *argv[]){

	//preparazione del generatore di numerica casuali
	Random rnd;
	SetRandomGenerator(&rnd);
   
	int M=100000;  //numero totale di dati generati            
	int N=100;     //numero di dati su cui mediare          
	
	//genero tutti i 10^5 numeri casuali sottraendoci 0.5 ed elevando al quadrato
	vector<double> s(M);
	
	for(int j=0; j<M; j++){
		s.at(j) = pow(rnd.Rannyu()-0.5,2);
	}

	//calcolo della media e deviazione standard della media utilizzando progressivamente sempre più blocchi di dati
	Analisi(N,M,s, "../risultati01.1.2.out");

	rnd.SaveSeed();

return 0;
}

//funzione per inizializzazione del generatore di numeri casuali
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


/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
