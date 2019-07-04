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
	vector<double> Integrali(M); //vettore per contenere le stime dell'energia variazionale ottenute integrnado l'hamiltoniana

	//preparo la distribuzione da campionare
	Psi* psi = new Psi(); 
	psi->SetMu(0.79);
	psi->SetSigma(0.4);

	//hamiltoniana
	FunzioneBase * H = new Hamiltoniana(psi);

	//posizione iniziale
	vector<double> posizione_0(1); 
	posizione_0.at(0) = 0.1;

	//equilibrazione
	double passo=0.7; 
	int dim; //dimensione del problema
	dim=posizione_0.size();
	Metropolis walker(psi,posizione_0, passo, dim);

	//salvataggio dati su file
	ofstream Risultati;
	Risultati.open("../../OptimizationParameters.out");
	

	//Calcolo degli integrali
	for(double sigma=0.55; sigma <0.65; sigma += 0.0125){

		psi->SetSigma(sigma);
		for(double mu=0.75; mu<0.85; mu += 0.025){
		psi->SetMu(mu);
		//vettore degli integrali
		vector<double> Integrali(M);

		//equilibrazione
		walker.Initialization(n_passi_metropolis);

			for(int i=0; i<M; i++){

				walker.Restart();
				walker.Sample(n_passi_metropolis);
				Integrali[i]=walker.Integration(H);
				//cout << "probabilità di accettare la nuova posizione: " << walker.AcceptanceProbability() << endl;
				//cout << "Valore integrale: " << Integrali[i] << endl;
		 	}
		double valore=Media(Integrali);
		Risultati << mu << " " << sigma << " " << valore << endl;
		cout << mu << " " << sigma << " " << valore << endl << endl;


		}
	}
	Risultati.close();
	delete psi;
	delete H;	

return 0;
}


