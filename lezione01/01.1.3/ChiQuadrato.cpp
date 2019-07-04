#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>

#include "ChiQuadrato.h"
#include "Statistica.h"

using namespace std;

//calcolo dellla media progressiva dei chi e delle deviazione standard della media progressive


void ChiQuadrato(int M, int N_tot, vector<double> s, const char * filename){

	//numero di dati generati casualmente sui cui calcolo i chi quadrati singoli
	int n = int(N_tot/M);

	//vettore contenente i chi quadrati singoli
	vector<double> chi(M);

	//vettore che contiene o e l'estremo di ogni intervallo in cui è diviso [0,1)
	vector<double> intervalli(M+1);
	for(int i=0; i<M+1; i++){	
		intervalli.at(i) = (1.*i)/M;
	}


	//il vettore n_i serve a contenere il numero di punti che cadadono in ognuno degli M intervalli. 
 	vector<double> n_i(M);

	//Segue il calcolo dei chi quadrati singoli
	for(int k=0; k<M; k++){

		//Ad ogni iterazione n_i viene riempito di zeri per iniziare il conteggio del numero di punti in ogni intervallo
		for(int b=0; b<M; b++){
			n_i.at(b) = 0.;
		}
	
		//calcolo del numero di punti che cadono in ognuno degli M intervalli
		  for(int i=k*n; i<(k+1)*n; i++){
			for(int j=0; j<M+1; j++){
		    		if( (s.at(i) >= intervalli.at(j) ) and ( s.at(i) < intervalli.at(j+1) ) ){
					n_i.at(j) = n_i.at(j)+1  ;
				}
			}
		  }
		

		
		//calcolo del chi quadrato di ognuno dei blocchi da n dati
		double somma= 0;
		double attesi = ((1.*n)/(1.*M)); 
		for(int a=0; a<M; a++){
			somma +=  pow( n_i.at(a)- attesi ,2) / attesi;			
		} 
		chi.at(k) = somma;
		
	}

	ofstream Risultati;
	Risultati.open(filename);

	for(int i=0; i<chi.size(); i++){
		Risultati<< i << " " << chi.at(i) << endl;
	}
		
	

}
	
