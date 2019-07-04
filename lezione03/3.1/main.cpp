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
#include "OptionPrice.h"
#include "Statistica.h"
#include "random.h"


using namespace std;
 
int main (int argc, char *argv[]){
 
	double t=0; //t iniziale
	double T=1.; //t finale
	int passi=100; //numero di passi in cui si divide l'intervallo T-t
	double s_0=100; //prezzo al tempo t=0
	double K=100; //strike price
	double r=0.1; //tasso di interesse
	double sigma=0.25; //volatilità

	OptionPrice asset(s_0, K, r, sigma); //costruzione della classe per valutare call e put option

	int M= 10000; //nuero totale di dati generati
	int N =100; //numero di blocchi

	//I vettori contentono rispettivamente M realizzazioni di: put-option dividendo l'intervallo [0,1] in 100 parti, call-option dividendo 		l'intervallo [0,1] in 100 parti, put-option passando riettamente da t=0 a T=1, call-option passando riettamente da t=0 a T=1
	vector<double> put_passi(M);
	vector<double> call_passi(M);
	vector<double> put(M);
	vector<double> call(M);

	//Valutazione di call e put option nei due casi
	for(int i=0; i<M; i++){
		put_passi.at(i) = asset.Put(t,T,passi);
		call_passi.at(i) = asset.Call(t,T,passi);
		put.at(i) = asset.Put(t,T,1);
		call.at(i) = asset.Call(t,T,1);
		
	}

	//Calcoo  della media e della deviazione datndard della media progressive
	Analisi(N,M,put_passi, "../risultati03.put_passi.out");
	Analisi(N,M,call_passi, "../risultati03.call_passi.out");
	Analisi(N,M,put, "../risultati03.put_un_passo.out");
	Analisi(N,M,call, "../risultati03.call_un_passo.out");
	
 
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
