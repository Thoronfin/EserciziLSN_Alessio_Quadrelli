#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include "random.h"
#include "Statistica.h"
#include "Salesman.h"

using namespace std;

double signum (double x){
    if (x>=0){return 1.;}else{return -1.;}
}

Salesman::Salesman(int NumberOfCity, int NumberOfParents): _BestPercorso(NumberOfCity){

  	//inizializzazione generatore numeri casuali
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
	   
	  
	if (Primes.is_open()){
	   Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();

	ifstream input("seed.in");
	string property;
	if (input.is_open()){
	   while ( !input.eof() ){
	      input >> property;
	       if( property == "RANDOMSEED" ){
		 input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
		 _rnd.SetRandom(seed,p1,p2);
	      }
	   }
	   input.close();
	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

  _NumberOfCity = NumberOfCity;
  _NumberOfParents = NumberOfParents;

  //vettore contenente l'indece delle città
  vector<int> cities(_NumberOfCity);
  for (int i=0; i<_NumberOfCity; i++){
    cities.at(i) = i;
  }

  //generazione della popolazione iniziale permutando casualmente l'inidce delle città, escludendo la prima che è sempre fissata
  srand(10);
  for(int i=0; i<_NumberOfParents; i++){
    random_shuffle(cities.begin()+1, cities.end());
    _StringsOfCities.push_back(cities);}

   _NewGeneration= _StringsOfCities;

   _accepted=0;
  _attempted=0;
}


Salesman::~Salesman(){
  _rnd.SaveSeed();
}


void Salesman::Square(void){ //inizializzazione delle città ne quadrato di lato 1
  double L=1.;
  for (int i=0; i<_NumberOfCity;i++){

    double x=_rnd.Rannyu(0,L);
    double y=_rnd.Rannyu(0,L);
    vector<double> tram(3);
    tram[0]=i; tram[1]=x; tram[2]=y;
    _Positions.push_back(tram);
  }

}



void Salesman::Circle(void){ //inizializzazione delle città lungo la circonferenza di raggio 1
  double R=1.;
  for (int i=0; i<_NumberOfCity;i++){

    double x=_rnd.Rannyu(-R,R);
    double y= sqrt(R*R -x*x)*signum(_rnd.Rannyu(-1,1));
    vector<double> tram(3);
    tram[0]=i; tram[1]=x; tram[2]=y;
    _Positions.push_back(tram);
  }


}



vector<double> Salesman::Fitness(vector<vector<int>> positions){ //calcolo della distanza su tutta la popolazione di percorsi

vector<double> Distance(positions.size());

  for(int j=0; j<positions.size(); j++){
    double L=0;
    for (int i=0; i<positions[0].size()-1;i++){
      int index = positions[j][i];
      int index_plus = positions[j][i+1];
      L+=pow(_Positions[index][1]-_Positions[index_plus][1],2) + pow(_Positions[index][2]-_Positions[index_plus][2],2);
    }
    L+=pow(_Positions[positions[j][0]][1]-_Positions[positions[j][_NumberOfCity-1]][1],2)+pow(_Positions[positions[j][0]][2]-_Positions[positions[j][_NumberOfCity-1]][2],2);
   Distance[j] = L;

  }
return Distance;
}


void Salesman::Mutation(double p ){ //permutazione casuale con probabilità p di una coppia di città in ognuno dei percorsi
  
  double random=_rnd.Rannyu(0.,1.);
  
  if(random<p){
    
    for (int i=0; i<_NewGeneration.size(); i++){
      unsigned int index_start=int(_rnd.Rannyu(1,_NumberOfCity));
      unsigned int index_end = int(_rnd.Rannyu(1,_NumberOfCity)); 
      int tram=_NewGeneration[i][index_start];
      _NewGeneration[i][index_start]=_NewGeneration[i][index_end];
      _NewGeneration[i][index_end]=tram;
    }
  }
  
 
  
}

double Salesman::BestDistance(vector<double> v){ //calcolo del valore minimo di un vettore

	return *min_element(v.begin(),v.end());
}


int Salesman::BestIndex(vector<double> v){ //calcolo dell'indice corrispondente al valore minimo di un vettore
	int index_minimum;
return  index_minimum= min_element(v.begin(),v.end()) - v.begin();

}

void Salesman::Metropolis(double beta, double p, int n_mutazion){
	
	//configurazioni iniziali
	vector<double> old_fitness= Fitness(_StringsOfCities); 
  int index_minimum=BestIndex(old_fitness);
  _BestDistance.push_back(old_fitness[index_minimum]);
  
	//nuove proposte
	_NewGeneration=_StringsOfCities; 
	for(int i=0; i<n_mutazion; i++){
	  Mutation(p);
	}
	
	vector<double> new_fitness= Fitness(_NewGeneration);
	double ratio;

	//Metropolis
	for(int i=0; i<_NumberOfParents; i++){

		ratio=exp( - beta*(new_fitness[i]-old_fitness[i]) );
		_attempted++;

		if(ratio>=1){
			_accepted++;
			_StringsOfCities[i]=_NewGeneration[i];
		}
		else{
			
			double p=_rnd.Rannyu();
			if(p<ratio){
				_StringsOfCities[i]=_NewGeneration[i];
				_accepted++;
				
			}
		}
	
	}
	
	//cout << "Acceptance Rate: " << ((1.*_accepted)/(1.*_attempted))  *100 << " %" << endl;

}

void Salesman::Sample(int n,double beta, double p, int n_mutazion){ // n iterazioni del Metropolis

	for(int i=0; i<n; i++){
		Metropolis(beta,p, n_mutazion);

	}
	
	vector<double> old_fitness= Fitness(_StringsOfCities); //configurazioni iniziali
  int index_minimum=BestIndex(old_fitness);
  _BestPercorso=_StringsOfCities[index_minimum];

}


void Salesman::Restart(){ //azzeramento del vettore delle distanze migliori e dei contatori

	_BestDistance.erase(_BestDistance.begin(),_BestDistance.end()); 
	_accepted=0;
	_attempted=0;
 
}


void Salesman::PrintDistances(const char* filename, int campionamenti){

	ofstream Risultati;
	Risultati.open(filename);
	vector<double> appo(campionamenti);
	
	for(int j=0; j<_BestDistance.size()/campionamenti; j++ ){
	  
  	for(int i=0; i< campionamenti; i++){
  	  appo[i]=_BestDistance[j*campionamenti+i];
  	}
  	
  	int minimum= BestIndex(appo);
	  Risultati << appo[minimum]<< endl;
	}
  	
  	Risultati.close();
}



void Salesman::PrintPath(const char* filename){

ofstream Path;
Path.open(filename);

	vector<double> v = Fitness(_StringsOfCities);
	int index_minimum = min_element(v.begin(),v.end()) - v.begin();
	int index;
	for(int i=0; i < _NumberOfCity; i++){
		index=_StringsOfCities[index_minimum][i];
		for(int j=0; j<3; j++){
		  Path << _Positions[index][j] << " ";
		}
	Path << endl;
	}
	Path << _NumberOfCity << " " << _Positions[0][1] << " " << _Positions[0][2] << endl; //aggiungo la prima città poichè il commesso torna da dove è partito

Path.close();

}
