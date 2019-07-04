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

void Analisi(int N, int M, vector<double> s, const char* namefile){

// riempio il vettore mean con le medie ottenuti da blocchi di M/N numeri e std con le medie al quadrato

	int L = int(M/N);
	vector<double> ai(N);
	vector<double> ai_2(N);
 	double a; // variabile di appoggio

  	for(int i=0; i<N;i++){
		a= Somma(s,i*L,(i+1)*L)/L;
		ai.at(i) = a;
		ai_2.at(i) = pow(a,2);
	}
  	
//calcolo la media progressiva del vettore mean e le deviazione standard della media progressiva
 	vector<double> mean_prog(N);
	vector<double> mean_2_prog(N);
	vector<double> mean_std(N);

	for(int i=0;i<N;i++){
		mean_prog.at(i) = Somma(ai,0,i+1)/(i+1);
		mean_2_prog.at(i) = Somma(ai_2,0,i+1)/(i+1);
		mean_std.at(i) = sqrt( ( mean_2_prog.at(i) - pow(mean_prog.at(i),2) ) / (i+1)  );
	}

// Stampo su file i risultati nel seguente formato: numero di dati usati per il calcolo di media e deviazione standard della media, media progressiva, deviazione standard della media progressiva
	ofstream Risultati;
	Risultati.open(namefile);

   for(int i=1;i<=N;i++){	
	
	if (Risultati.is_open()){
   		Risultati <<i*L<< " "<< mean_prog.at(i-1) << " " <<mean_std.at(i-1)<< endl;
   } 
	else cerr << "PROBLEM: Unable to open the file.out" << endl;
   	
   }	
	
   Risultati.close();	
}


void AnalisiMatrici(int N,vector<vector<double>> s,double fattore_scala, double traslazione_inidice, const char* namefile){

	//stream per il salvataggio dei dati su file
	ofstream Risultati;
	Risultati.open(namefile, ios::app);

	//parametri per il data blocking
	int nrighe=s.size(); // numero di vettori( righe della matrice) su cui si effettua il data blocking
	int M= s[0].size(); // lunghezza di ogni vettore su cui si effettua data blocking
	int L = int(M/N); // lunghezza di ogni blocco

	for(int j=0; j<nrighe; j++){ //data blocking su ogni riga della matrice

		vector<double> ai(N);
		vector<double> ai_2(N);
	 	double a; // variabile di appoggio

	  	for(int i=0; i<N;i++){
			a= Somma(s[j],i*L,(i+1)*L)/(L);
			ai.at(i) = a;
			ai_2.at(i) = pow(a,2);
		}
	  	
	//calcolo la media progressiva del vettore mean e le deviazione standard della media progressiva
	 	vector<double> mean_prog(N);
		vector<double> mean_2_prog(N);
		vector<double> mean_std(N);

		for(int i=0;i<N;i++){
			mean_prog.at(i) = (Somma(ai,0,i+1)/(i+1));
			mean_2_prog.at(i) = Somma(ai_2,0,i+1)/(i+1) ;
			mean_std.at(i) = sqrt( ( mean_2_prog.at(i) - pow(mean_prog.at(i),2) ) / (i+1)  );
		}

	// Stampo su file i risultati nel seguente formato:indice, media progressiva, deviazione standard della media progressiva con il massimo numero di blocchi
		
		if (Risultati.is_open()){
	   		Risultati << j*fattore_scala - traslazione_inidice << " "<< mean_prog.at(N-1) << " " <<mean_std.at(N-1)<< endl;
	  	 } 
		else cerr << "PROBLEM: Unable to open the file.out" << endl;
	 	
	 	
	}

   Risultati.close();	
}

