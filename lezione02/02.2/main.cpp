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
#include <cmath>
#include <vector>
#include "Rw.h"
#include "Statistica.h"
#include "random.h"


using namespace std;
 
int main (int argc, char *argv[]){
  
	int N=100; //numero di blocchi
	int M=1000000;  //numero totali di random walk effettuati 10^4 realizzazione per i passi (1,2,3,...,100) in toale sono 10^6
	vector<double> d(M); //distanza percorsa
	vector<double> d_2(M); //distanza percorsa al quadrato

	Rw w(1); //il passi viene fissato ad 1
   	int posizione;
	
	//Reticolo cubico
	for(int i=0; i< N;i++){

		for(int j=0; j<(M/N); j++){
			posizione = j +i*(M/N);
			d[posizione]=w.ReticoloCubico(i);
			d_2[posizione]=pow(d[posizione],2);			
	   	}

	}

	//savataggio su file della radice quadrata del modulo della distanza percorsa
	Analisi(N,M,d_2,"../risultati02.2.1.out");


	//continuo
	for(int i=0; i< N;i++){

		for(int j=0; j<(M/N); j++){
			posizione = j +i*(M/N);
			d[posizione]=w.Diffusione(i);
			d_2[posizione]=pow(d[posizione],2);
			
			
	   	}

	}

	//savataggio su file della radice quadrata del modulo della distanza percorsa
	Analisi(N,M,d_2,"../risultati02.2.2.out");

	
  
return 0;
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
