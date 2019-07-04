#include <vector>
#include <cmath>
#include<iostream>
#include<fstream>
#include<string>

#include "Statistica.h"

using namespace std;

double Somma(vector<double> v, int i, int j){

	double somma=0;
	for(int a=i; a<j; a++){
		somma += v.at(a);
	}
	return somma;
}

double Media(vector<double> v) {	
	double somma=0;
	for(unsigned int i=0; i<v.size(); i++) {
		somma+=v.at(i);
	}

	return somma/v.size();

}

double Varianza(vector<double> v) {
	double somma=0;	
	for(unsigned int i=0; i<v.size(); i++){
		somma+=pow(Media(v)-v.at(i), 2.);
	}
	
	return somma/(v.size()-1) ;
 }

void Analisi(int N, int M, vector<double> v, const char* namefile){
 // riempio il vettore mean con le medie ottenuti da blocchi di M/N numeri e std con le medie al quadrato

ofstream Risultati;
Risultati.open(namefile);

int L = int(M/N);
int part = L/N;

for (int k=0; k<N; k++){

	
	vector<double> mean(part);
	vector<double> std(part);

 	double a; // variabile di appoggio
  	for(int i=0;i<part; i++){

		a= sqrt(Somma(v, k*L + i*part, k*L + (i+1)*part)/part);
		mean.at(i)= a;
		std.at(i)= pow(a,2);
	}
  	
//calcolo la media progressiva del vettore mean e le deviazione standard della media progressiva
 	double mean_tot;
	double std_tot;
	double mean_std;

	
	mean_tot= pow(Media(mean),2);
	std_tot= Somma(std,0,part)/part;
	mean_std= sqrt ( (std_tot-mean_tot) / (part-1) ) ;


// Stampo su file i risultati nel seguente formato: numero di dati usati per il calcolo di media e deviazione standard della media, media progressiva, deviazione standard della media progressiva
	

   
	
	if (Risultati.is_open()){

   		Risultati << k  << " "<< sqrt(mean_tot) << " " << mean_std << endl;
   	} 
	else cerr << "PROBLEM: Unable to open the file.out" << endl;
   	
  }	
	

  Risultati.close();	
}

