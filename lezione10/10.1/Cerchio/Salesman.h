#ifndef _Salesman_h
#define _Salesman_h

#include "random.h"
#include <vector>


using namespace std;

class Salesman{

private:

  Random _rnd;
  int _NumberOfCity; //numero delle città
  int _NumberOfParents; //popolazione
  vector<vector<int>> _StringsOfCities; //lista dei percorsi
  vector<vector<double>> _Positions; //lista delle posizioni (x,y) e indice i di ogni città
  vector<vector<int>> _NewGeneration; // lista dei percorsi proposta per il metropolis
  vector<double> _BestDistance; //Vettore della migliore distanza ottenute ad ogni passo del metropolis
  vector<int> _BestPercorso; //matrice le cui righe sono le migliori configurazioni ad ogni passo del Metropolis
  int _accepted; //nuove configurazioni accettate
  int _attempted; //configurazioni totali usate
  int _n_mutazioni; //numero di mutazioni su ogni genitore che viene usata per generare la nuova configurazione da sottoporre al Metropolis


public:

  Salesman(int NumberOfCity, int NumberOfParents);
  ~Salesman();

//Generazioni delle posizioni delle città
  void Square(void);
  void Circle(void);

 
  vector<double> Fitness(vector<vector<int>>); //calcolo del vettore delle distanze corrispondenti ad ognuna delle righe della matrice
  void Mutation(double p); //Mutazione con probabilita` p che serve per generare i nuovi percorsi per il metropolis
  void Metropolis(double beta, double p, int n_mutazion); //campionamento con l'algoritmo di metropolis a beta fiato
  void Sample(int n,double beta, double p, int n_mutazion); // n iterazioni del metropolis
  double BestDistance( vector<double>); // calcolo del minimo di un vettore di distanze
  void PrintDistances(const char*, int campionamenti); //stampa del vettore BestDistance
  void Restart(); //azzermaneto di contatori e del vettore BestDistances
  void PrintPath(const char*); //stampa del perocrso corrispondente alla minima distanza
  int BestIndex(vector<double>); //restituisce l'indice dell'elemento corrispondente al minimo del vettore
 
 

};

#endif
