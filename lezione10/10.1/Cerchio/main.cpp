#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include "random.h"
#include "Salesman.h"

using namespace std;

int main (){

	
	int campionamenti=2000; //numero di campionamenti per ogni beta
  	double Probaility_mutation= 1; //probabilità di mutazione
 	int n_mutazioni=1; //numero delle mutazioni

	double Beta_iniziale=0.; // T infinita
	double Beta_finale=40.;  // T bassa
	double DeltaBeta=0.006; //
	
	//commesso viaggiatore con 30 città e una popolazione di 1 elemento
	Salesman Walker(30,1);
	//città in un cerchio di raggio 1
	Walker.Circle();

	for(double Beta=Beta_iniziale; Beta<Beta_finale; Beta+=DeltaBeta){ //campionamento con il Metropolis per ogni beta
		cout << "beta: " << Beta << endl;
		
		Walker.Sample(campionamenti, Beta,Probaility_mutation, n_mutazioni); //campionamento di 2000 configurazioni a beta fissato
		cout << endl;

	}

	Walker.PrintDistances("../Distanze_cerchio.dat", campionamenti); //salvataggio della migliore distanza ottenuta da 2000 campionamenti
	Walker.PrintPath("../BestPath_cerchio.dat"); //salavataggio del miglior percorso


  return 0;

}
