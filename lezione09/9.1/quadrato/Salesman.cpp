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

Salesman::Salesman(int NumberOfCity, int NumberOfParents): _FitnessDistances(NumberOfParents) {

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

 
  //vettore contentenete l'inidice delle città
  vector<int> cities(_NumberOfCity);
  for (int i=0; i<_NumberOfCity; i++){
    cities.at(i) = i;
  }

  //la popolazione iniziale è ottenuta dalla permutazione casuale del vettore contenente l'inidce delle città
  srand(10);
  for(int i=0; i<_NumberOfParents; i++){
    random_shuffle(cities.begin()+1, cities.end());
    _StringsOfCities.push_back(cities);}

}


Salesman::~Salesman(){
  _rnd.SaveSeed();
}


void Salesman::Print(void){ //stampa a video dei cromosomi

  for(int j=0; j<_NumberOfParents; j++){
		for(int i=0; i<_NumberOfCity;i++){
			cout << _StringsOfCities[j][i] << " ";
		}
    cout << endl;
	}

}

void Salesman::Square(void){
  double L=1.;
  for (int i=0; i<_NumberOfCity;i++){

    double x=_rnd.Rannyu(0,L);
    double y=_rnd.Rannyu(0,L);
    vector<double> tram(3);
    tram[0]=i; tram[1]=x; tram[2]=y;
    _Positions.push_back(tram);
  }

}



void Salesman::Circle(void){
  double R=1.;
  for (int i=0; i<_NumberOfCity;i++){

    double x=_rnd.Rannyu(-R,R);
    double y= sqrt(R*R -x*x)*signum(_rnd.Rannyu(-1,1));
    vector<double> tram(3);
    tram[0]=i; tram[1]=x; tram[2]=y;
    _Positions.push_back(tram);
  }


}


void Salesman::Fitness(void){ //calcolo dell fitness L^(2)
  for(int j=0; j<_NumberOfParents; j++){
    double L=0;
    for (int i=0; i<_NumberOfCity-1;i++){
      int index = _StringsOfCities[j][i];
      int index_plus = _StringsOfCities[j][i+1];
      L+=pow(_Positions[index][1]-_Positions[index_plus][1],2) + pow(_Positions[index][2]-_Positions[index_plus][2],2);
    }
    L+=pow(_Positions[_StringsOfCities[j][0]][1]-_Positions[_StringsOfCities[j][_NumberOfCity-1]][1],2)+pow(_Positions[_StringsOfCities[j][0]][2]-_Positions[_StringsOfCities[j][_NumberOfCity-1]][2],2);
    _FitnessDistances[j] = L;

  }
}

void Salesman::Crossover(int parent1, int parent2){
  
  unsigned int index_start = int(_rnd.Rannyu(1,_NumberOfCity)); //0.9999
  unsigned int index_end = int(_rnd.Rannyu(index_start,_NumberOfCity)); //0.9999

  vector<int> son(_NumberOfCity, 0);
  for(int i=index_start; i<index_end; i++){
    son.at(i) = _StringsOfCities[parent1][i];
    
  }

//Prima parte del figlio
  int actual_position=0;
  for(int i=0; i<index_start; i++){
    
    for(int j=actual_position; j<_NumberOfCity; j++){//cambiamento
      
      int counter = 0;
      for(int k=index_start; k<index_end; k++){

          if(_StringsOfCities[parent2][j] == son.at(k)){
             
              counter++;}

      }
       
      if(counter == 0){
          
          son.at(i) = _StringsOfCities[parent2][j];
          actual_position = j+1;//crucial
          break;
      }
    }
  }
//Seconda parte del figlio
  for(int i=index_end; i<_NumberOfCity; i++){
    //cout << i << endl;
    for(int j=actual_position; j<_NumberOfCity; j++){
        int counter = 0;
      for(int k=index_start; k<index_end; k++){
        if(_StringsOfCities[parent2][j] == son.at(k)) counter++;
      }
        if(counter == 0){ son.at(i) = _StringsOfCities[parent2][j];actual_position = j+1; break;}//crucial
    }
  }


  _NewGeneration.push_back(son);

}


void Salesman::Mutation(double p ){ //mutazione delle nuova generazione
 
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
    

    return;
    
    
}




//selection dà solo un genitore
int Salesman::Selection(void){
   
    vector<double> sorted_distances= _FitnessDistances;
    sort(sorted_distances.begin(), sorted_distances.end()); //distanze ordinate in ordine crescente

    double p=_rnd.Rannyu(0,_NumberOfCity); // vengono tenuti come possibili genitori quelli con le prime p distanze minime
    double minimum;
    vector<int> index;
    for(int i=0; i<p;i++){
        
       minimum=sorted_distances[i];
       vector<double>::iterator index_minimum= find(_FitnessDistances.begin(), _FitnessDistances.end(), minimum);
       index.push_back( index_minimum - _FitnessDistances.begin());
    }


  //degli elementi selezionati ne viene scelto uno casualmente
    double p2=int(_rnd.Rannyu(0,p));
    if ((index.size())==0){

        return min_element(_FitnessDistances.begin(),_FitnessDistances.end()) - _FitnessDistances.begin();
	
    }else{

    return index[p2];
     }
}


void Salesman::Run(double n_generation,double probability_mutation){
    
    for(int i=0; i <n_generation;i++){
	
        Fitness(); //calcolo della fitness
        HalfPopulationMean(); //valor medio della distanza della metà migliore della popolazione
	
        for(int j=0; j<_NumberOfParents; j++){ //generazione dei figli dalla selezione di due genitori
            Crossover( Selection(),  Selection());
            
        }
        
     
      Mutation(probability_mutation); //mutazione delle nuova generazione
    
      _StringsOfCities=_NewGeneration; //aggiornamento della generazione attuale
      _NewGeneration.erase(_NewGeneration.begin(),_NewGeneration.end());
      
  
        
      }
    
    
    
    return;
    
}

void Salesman::HalfPopulationMean(){
  
  //selezione della migliore metà della popolazione
  vector<double> sorted_distances = _FitnessDistances;
  sort(sorted_distances.begin(), sorted_distances.end()); //distanze ordinate in ordine crescente
  vector<double> appo(_NumberOfParents/2);
  //riempimento di un vettore con le distanze percorse dalla metà migliore della popolazione 
  for(int j=0; j< _NumberOfParents/2; j++){
    vector<double>::iterator index_minimum = find(_FitnessDistances.begin(), _FitnessDistances.end(), sorted_distances[j]);
    int indice = index_minimum - _FitnessDistances.begin();
    appo[j]=_FitnessDistances[ indice ];
  }
  
  //calcolo della media della generazione attuale
  _HalfPopolutionMean.push_back(Media(appo));
  //cout << Media(appo) << endl;
}


void Salesman::PrintMean(const char* filename, int incremento){
  ofstream Risultati;
  Risultati.open(filename);
  
  double mean=0.;
  double devstd=0.;
  for(int i=0; i<_HalfPopolutionMean.size()/incremento; i++){
    
      mean = Media(_HalfPopolutionMean, i*incremento, (i+1)*incremento) ;
      Risultati << (incremento*i)+1 << " " << mean << endl;
  }
  
  Risultati.close();
}

void Salesman::Restart(){

 //la popolazione viene inizializzata casualmente
  vector<int> cities(_NumberOfCity);
  for (int i=0; i<_NumberOfCity; i++){
    cities.at(i) = i;
  }

  for(int i=0; i<_NumberOfParents; i++){
    random_shuffle(cities.begin()+1, cities.end());
    _StringsOfCities[i]=cities;
   _FitnessDistances[i]=0.;
  }
 //le distanze vengono azzerate
 
//la nuova generazione viene cancellata
  //_FitnessDistances.erase(_FitnessDistances.begin(),_FitnessDistances.end());
  _NewGeneration.erase(_NewGeneration.begin(),_NewGeneration.end()); 
}


double Salesman:: BestDistance(){


	return *min_element(_FitnessDistances.begin(),_FitnessDistances.end());

}

void Salesman::BestPath(const char* filename){

ofstream Path;
Path.open(filename);

	int index_minimum = min_element(_FitnessDistances.begin(),_FitnessDistances.end()) - _FitnessDistances.begin();
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
