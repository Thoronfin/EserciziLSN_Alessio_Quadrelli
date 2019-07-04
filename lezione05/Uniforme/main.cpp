#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "Funzioni.h"
#include "Statistica.h"
#include "Metropolis.h"
using namespace std;


int main (int argc, char *argv[]){ 

	int M=1000000;  //numero dei campionamenti    
	int prove=1000; // numero dei campionamenti per equlibrazione    
   	int N=100;     //numero di dati su cui mediare 

	double passo=1.2; 
	 

	//preparo la distribuzione da campionare
	FunzioneBase* orbitale_uno_s = new uno_s(); 
	vector<double> posizione_0(3); 
	//posizione iniziale
	posizione_0.at(0) = 0.1;
	posizione_0.at(1) = 0.1;
	posizione_0.at(2) = 0.1;
	//campionamento
	int dim=posizione_0.size(); //dimensione del problema
	Metropolis walker(orbitale_uno_s,posizione_0, passo, dim);
	walker.Initialization(prove); //equilibrazione
	walker.Sample(M); //campionamento 
	//probabilità di accettazione
	cout << "probabilità di accettare un la nuova posizione nell'orbitale 1s: " << walker.AcceptanceProbability() << endl; 
	
	//salvataggio delle configurazioni campionate
	walker.Print("../risultati05.1.posizioni.out");
	//calcolo della distanza media dall'origine.
	vector<double> raggi=walker.Raggi();
	Analisi(N,M,raggi,"../risultati05.1.raggi.out");

	delete orbitale_uno_s;

	
	passo=2.97;

	FunzioneBase* orbitale_due_p = new due_p(); //preparo la distribuzione da campionare
	//posizione iniziale
	posizione_0.at(0) = 0.1;
	posizione_0.at(1) = -0.5;
	posizione_0.at(2) = 0.1;
	dim=posizione_0.size();
	//campionamento
	Metropolis walker1(orbitale_due_p,posizione_0, passo, dim);
	walker1.Initialization(prove);
	walker1.Sample(M);
	//probabilità di accettazione
	cout << "probabilità di accettare un la nuova posizione nell'orbitale 2p: " << walker1.AcceptanceProbability() << endl;
	//salvataggio delle configurazioni campionate
	walker1.Print("../risultati05.2.posizioni.out");
	//calcolo della distanza media dall'orignie
	raggi=walker1.Raggi();
	Analisi(N,M,raggi,"../risultati05.2.raggi.out");
	
  	delete orbitale_due_p;

return 0;
}


