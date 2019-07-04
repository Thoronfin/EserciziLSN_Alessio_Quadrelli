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
  vector<int> _BestPercorso; //vettore con l'indece della configurazione migliore di ogni passo del Metropolis
  int _accepted; //numero di configurazioni accettate
  int _attempted; //numero totale delle configurazioni
  int _n_mutazioni; //numero di mutazioni su ogni genitore che viene usata per generare la nuova configurazione da sottoporre al Metropolis


public:

  Salesman(int NumberOfCity, int NumberOfParents, int seme); //il seme serve per la genrazione di configurazioni iniziali diverse
  ~Salesman();

//Generazioni delle posizioni delle città
  void Square(void); //all'interno di un quadraro di lato 1
  void Circle(void); //lungo una circonferenza di raggio 1

 
  vector<double> Fitness(vector<vector<int>>); //calcolo del vettore delle distanze corrispondenti ad ognuna delle righe della matrice
  void Mutation(int n_elem_mutazione); //Mutazione che serve per generare i nuovi percorsi per il metropolis
  void Metropolis(double beta); //algoritmo di Metropolis
  void Sample(int n,double beta); // n iterazioni del metropolis
  void Equilibration(int n,double beta); //n passi di equilibrazione per il metropolis e azzermaneto di contatori e del vettore BestDistances
  double BestDistance( vector<double>); // calcolo del minimo del vettore
  void SetNumberMutation(int n); // permette di settare il numero di mutazioni(permutazioni di coppie di città) da effettuare su ogni percorso
  void PrintDistances(const char*); //stampa del vettore BestDistance
  void Restart(); //azzermaneto di contatori e del vettore BestDistances
  void PrintPath(const char*); //stampa del perocrso corrispondente alla minima distanza
  int BestIndex(vector<double>); //restituisce l'inidce dell'elento minimo del vettore
  double MinimumDistance(); //resisituisce la minima distanza percrosa
 
 

};

#endif
