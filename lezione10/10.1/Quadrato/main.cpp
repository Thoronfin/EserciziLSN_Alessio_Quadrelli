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

	
	int campionamenti=2000; //numero di campionamenti per ogni beta
  	double Probaility_mutation= 1.; //probabilità di mutazione
	int n_mutazioni=1; //numero delle mutazioni

	double Beta_iniziale=3.;
	double Beta_finale=40.;
	double DeltaBeta=0.009;
	
	//commesso viaggiatore con 30 città e una popolazione di 40 elementi
	Salesman Walker(30,1);
	//città in un quadrato di lato 1
	Walker.Square();
\

	for(double Beta=Beta_iniziale; Beta<Beta_finale; Beta+=DeltaBeta){
		cout << "beta: " << Beta << endl;
		
		Walker.Sample(campionamenti, Beta, Probaility_mutation, n_mutazioni);
		cout << endl;

	}

	Walker.PrintDistances("../Distanze_quadrato.dat",campionamenti);  //salvataggio della migliore distanza ottenuta da 2000 campionamenti
	Walker.PrintPath("../BestPath_quadrato.dat"); //salavataggio del miglior percorso

	

	

	


  return 0;

}
