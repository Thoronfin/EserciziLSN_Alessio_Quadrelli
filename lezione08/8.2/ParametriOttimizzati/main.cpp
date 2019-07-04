#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "Funzioni.h"
#include "Statistica.h"
#include "Metropolis.h"
using namespace std;


int main (int argc, char *argv[]){ 

	//Dati
	int M= 10000; // numero totale di integrali da calcolare
	int N= 100; // numero di blocchi
	int n_passi_metropolis =1000; //campionamenti del metropolis
	int nbins= 200; //bin per istogramma

	vector<double> Integrali(M); //vettore che contiene l'energia variazionale ottenuta integrando l'hamiltoniana

	//Distribuzione da campionare
	Psi* psi = new Psi(); 
	psi->SetMu(0.8);
	psi->SetSigma(0.6125);

	//Hamiltoniana applicata a psi fratto psi
	FunzioneBase * H = new Hamiltoniana(psi);

	//posizione iniziale e passo
	vector<double> posizione_0(1); 
	posizione_0[0] = 0.5;
	double passo=2.6; 

	//Inizializzazione
	Metropolis walker(psi,posizione_0, passo, M, nbins);

	//equilibrazione
	walker.Initialization(n_passi_metropolis);

	//Vengono realizzate M stime dell'energia del ground state ed M istogrammi del modulo quadro della funzione d'onda
	for(int i=0; i<M; i++){

		walker.Restart();
		walker.Sample(n_passi_metropolis);
		Integrali[i]=walker.Integration(H);
		walker.Histogram(i);
	}
	//probabilità di accettazione del passi
	cout << "probabilità di accettare la nuova posizione: " << walker.AcceptanceProbability() << endl;

	//salvataggio su file del data blocking degli istogrammi e dell'energia
	walker.AnalisiHistogram(N,"../../Configurazioni.out");
	Analisi(N,M,Integrali,"../../IntegraliOttimizzati.out");

return 0;

}


