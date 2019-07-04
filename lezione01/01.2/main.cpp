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
  
	//preparazione stream per salvare su file i risultati
	ofstream Risultati;
	Risultati.open("../risultati01.2.out");

	int F[4] = {1,2,10,100}; //numero dei dati usati per calcolare le medie
	int M=10000; //numero degli esperimenti effettuati per ognuno degli F[i] fissati
 
	// vettori usati per contenere le medie dei risultati degli esperimenti
	vector<double> dado_mean(M);
	vector<double> exp_mean(M);
	vector<double> lor_mean(M);
   
  
	for(int m=0;m<4;m++){
  	 // preparazione dei vettori che conterrano i singoli elementi (che saranno progressiavmente 1,2,10,100) generati casualmente, sui quali poi verrà effettuata la media
		vector<double> lor(F[m]);
		vector<double> exp(F[m]);
		vector<double> s(F[m]);
  
		for(int i=0;i<M;i++){
			// riempio i vettori del corrispondente numero di elementi casuali generati con le apposite distribuzioni
		   	for(int j=0; j<F[m]; j++){		
				s.at(j) = int(rnd.Rannyu(1,7));
				exp.at(j) = rnd.Exponential(1.);
				lor.at(j) = rnd.Lorentz(1.,0.);
		   	}
	   	//riempio i vettori che contengono tutti i risultati, cioè le medie
			dado_mean.at(i) = Somma(s,0,F[m])/F[m];
			exp_mean.at(i) = Somma(exp,0,F[m])/F[m];
			lor_mean.at(i) = Somma(lor,0,F[m])/F[m];

		   	//stampo i risultati su un file con il seguente ordine: lancio dei dadi con distribuzione uniforme, esponenziale e lorentziana.
		   if (Risultati.is_open()){
		   		Risultati <<dado_mean.at(i)<< " "<< exp_mean.at(i) << " " <<lor_mean.at(i)<< endl;
		   } 
			else cerr << "PROBLEM: Unable to open random.out" << endl;
		   
		   }	
   
	}
	Risultati.close();	

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

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
