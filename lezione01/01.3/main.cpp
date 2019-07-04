#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "Buffon.h"
#include "Statistica.h"

using namespace std;

int main (int argc, char *argv[]){
  
	Buffon Esperimento(6,4,50000); //preparo l'esperimento con due linee a distanza 6 e un ago di lunghezz 4, per ogni stima del pi uso 5*10^4 lanci
	int M=10000;  //numero totale di dati generati            
   	int N=100;     //numero di dati su cui mediare          
  
	// riempio il vettore s con M realizzazioni dell'espeimento di Buffon
	vector<double> s(M);
  	for(int j=0; j<M; j++){
		s.at(j) = Esperimento.pi();
  	}

	//data blocking sugli M esperinenti con N blocchi e salvataggio dei risultati su file
	Analisi(N,M,s, "../risultati01.3.out");
	

return 0;
}
