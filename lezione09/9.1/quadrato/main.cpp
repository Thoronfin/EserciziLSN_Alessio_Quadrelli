#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "Statistica.h"
#include "Salesman.h"

using namespace std;

int main (){

	double Probabilita_mutazione=0.08; //0.08 per square
	int N_generazioni=20000; //numero di generazioni totali
	int incremento=50; //incremento delle generazioni ad ogni passo
	
	//commesso viaggiatore con 30 città e una popolazione di 40 elementi
	Salesman Walker(30,40);
	//città in un quadrato di lato 1
	Walker.Square();

	ofstream Distanze;
	Distanze.open("../../Distanze_quadrato.dat"); // file per il salvataggio della miglior distanza ad ogni 'incremento' generazioni
	double best_distance;

	

	int index=0;
	for(int i=1; i<=N_generazioni+1; i+= incremento){
  		
		  Walker.Run(incremento, Probabilita_mutazione); //evoluzione per i generazioni
		  best_distance=Walker.BestDistance(); //cromosoma del gene con distanza minima dopo l'evoluzione
		  Distanze << i << " " << best_distance << endl;
		  index++;

		if(i == N_generazioni+1){
			Walker.BestPath("../../BestPath_quadrato.dat");  //salavataggio del miglior percorso
		  	Walker.PrintMean("../../Mean_quadrato.dat", incremento); //salvataggio della media della distanza percorsa ad ogni generazione dalla metà migliore della popolazione
		}
		
	}


  return 0;

}
