#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <algorithm>
#include "random.h"
#include "Statistica.h"
#include "Salesman.h"
#include "mpi.h"

using namespace std;

int main (int argc,char* argv[]){
  
  MPI::Init(argc,argv);
  int size = MPI::COMM_WORLD.Get_size();
  int rank = MPI::COMM_WORLD.Get_rank();
  
  int index_minimum; //nodo che contiene la configurazione di minima distanza
  double minimaldistance; //valore delle minima distanza
  double distance; //distanza calcolata da ogni processo
  vector<double> distances(size); //vettore per conterenere le distanze di ognuno dei processi
  double appoggio; //variabile di appoggio
  
  //Salvataggio dei risultati
  ofstream Risultati;
  Risultati.open("../Distanze_cerchio.dat");
  
  Salesman Walker(30,1, rank); //30 citta e una poplazione di un elemento
  Walker.Circle(); //iniziliazzazione città
  Walker.SetNumberMutation(1); //1 permutazioni di coppia ad ogni passo
  
  

    int campionamenti=2000;
    double Beta_iniziale=0.;
    double Beta_finale=20.;
    double DeltaBeta=0.006;
    
    for(double Beta=Beta_iniziale; Beta<Beta_finale; Beta+=DeltaBeta){
     
      Walker.Sample(campionamenti, Beta); //ogni nodo campiona
      distance = Walker.MinimumDistance(); //ogni nodo scrive la propria miglio distanza raggiunta finora
      
      
                        for(int i=1; i<size; i++){ //ogni nodo (escluso il nodo zero) invia al nodo zero la propria miglior distanza
                          
                          if(rank==i){
                            
                            MPI::COMM_WORLD.Send(&distance,1,MPI::DOUBLE,0,1);
                            
                          }
                          else{
                            if(rank==0){ // il nodo zero riceve tutte le distanze
                              MPI::COMM_WORLD.Recv(&appoggio,1,MPI::DOUBLE,i,1);
                              distances[i] = appoggio;
                            }
                          }
                          
                        }
                        
                        // il nodo zero calcola la minima distanza e il nodo che ha raggiunto la minima.  Il nodo inidcato salva la distanza 				su file
                        if(rank == 0){ 
                          cout << Beta << endl;
                          distances[0]=distance;
                          
                          index_minimum= min_element(distances.begin(),distances.end()) - distances.begin();
                          minimaldistance = distances[index_minimum];
                          
                          Risultati << minimaldistance << endl;
                        }
                        
      
    }
      
      
      
      //Walker.Restart();
     
    //Il nodo zero invia a tutti i nodi l'indice del miglior nodo che ha identificato 
    MPI_Bcast(&index_minimum,1,MPI::INTEGER,0,MPI::COMM_WORLD);
  

    //il nodo corrispondente all'inidice salva su file il miglior percorso
    if( rank==index_minimum){
              Walker.PrintPath("../BestPath_cerchio.dat");
    }
       
      
        
    
    Risultati.close();
 
      
  MPI::Finalize();

  return 0;

}
