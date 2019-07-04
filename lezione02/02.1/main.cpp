#include <iostream>
#include <fstream>
#include <cmath>
#include "random.h"
#include "Integrator.h"
#include "Funzioni.h"
#include "Statistica.h"
using namespace std;



int main (int argc, char *argv[]){


	int M=10000;  //numero totale delle volte usate per calcolare integrale         
   	int N=100;     //numero di dati su cui mediare
	int punti = 10000; // numero di dati utilizzati per valutare l'integrale

	//Esercizio 2.1	

	FunzioneBase* coseno = new Coseno(); //preparo la funzione da integrare
	Integral Integrale(0,1, coseno); //preparo la classe per calcolare l'integrale da 0 a 1 di (pi/2)*cos[(pi/2)*x]

	// riempio il vettore s con M calcoli dell'integrale realizzati con 10^4 punti generati uniformemente in [0,1)
	vector<double> s(M);
  	for(int j=0; j<M; j++){
		s.at(j) = Integrale.Media(punti);
  	}

	//Analizzo i dati ottenuti facendo media e deviazione standard della media 
	//sommando progressivamente le N medie ottenute da int(M/N) risultati dell'integrale
	//Poi i dati vengono salvati su file
	Analisi(N,M,s,"../risultati02.1.1.out");       


	// Esercizio 2.2

	FunzioneBase* g = new G(); //preparo la funzione da integrare g
	Integral Integrale1(0,1, g); //preparo la classe per calcolare l'integrale da 0 a 1 di (pi/2)*cos[(pi/2)*x]/ 2*(1-x)

	// riempio il vettore s con M calcoli dell'integrale realizzati con 10^4 punti generati con distribuzione p(x)= 2*(1-x) con x 	 nell'intervallo [0,1)
	
  	for(int j=0; j<M; j++){
		s.at(j) = Integrale1.MediaRetta(punti);
  	}

	//Analizzo i dati ottenuti facendo media e deviazione standard della media 
	//sommando progressivamente le N medie ottenute da int(M/N) risultati dell'integrale.
	//Poi i dati vengono salvati su file
	Analisi(N,M,s,"../risultati02.1.2.out");
  

return 0;
}
