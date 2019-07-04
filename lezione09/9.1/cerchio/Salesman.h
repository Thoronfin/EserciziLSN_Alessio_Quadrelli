#ifndef _Salesman_h
#define _Salesman_h

#include "random.h"
#include <vector>


using namespace std;

class Salesman{

private:

  Random _rnd; //generatore di numeri casuali
  int _NumberOfCity; //numero delle città cioè dei geni di ogni cromosoma
  int _NumberOfParents; //numero dei genitori ovvero numero di cromosomi
  vector<vector<int>> _StringsOfCities; //matrice contenente tutti i cromosomi
  vector<vector<double>> _Positions; //vettore contenente le posizioni (x,y) delle città
  vector<double> _FitnessDistances; //vettore contennete la fitness di ogni cromosoma
  vector<vector<int>> _NewGeneration; //matrice delle nuova generazione
  vector<double> _HalfPopolutionMean; //vettore usato per contenere il valor medio della distanza della metà migliore della popolazione


public:

  Salesman(int NumberOfCity, int NumberOfParents);
  ~Salesman();
//Generazioni delle posizioni in un quadrato di lato 1 o lungo una circonferenza di raggio 1
  void Square(void);
  void Circle(void);

  void Fitness(void); //calcolo della loss function di ogni cromosoma
  void Crossover(int parent1, int parent2); //metodo di cross-over tra due genitori
  void Mutation(double p); //mutazione con probabilità p
  int Selection(void); //metrodo di selezione di un singolo genitore

  void Run(double n_generation,double probability_mutation); //esegue l'evoluzione della poloazione per n generazioni con probailità di mutazione arbitraria
  void Restart();


  void BestPath(const char*); //salvataggio sul  file del miglio perrcorso
  void HalfPopulationMean(); //valor medio della distanza della metà migliore della popolazione
  void Print(void); //stampa a video dei cromosomi
  double BestDistance(); //resistuisce la miglior distanza raggiunta dalla generazione corrente
  void PrintMean(const char*, int n); //salvataggio sul  file del valor medio della distanza della metà migliore della popolazione
  

};

#endif
