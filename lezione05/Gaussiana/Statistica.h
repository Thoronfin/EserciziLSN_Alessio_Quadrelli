#ifndef _Statistica_h_
#define _Statistica_h_

#include <vector>

using namespace std;

double Somma(vector<double> v, int i, int j); //somma gli elementi di v dall'i-esimo incluso, al j-esimo escluso
double Media(vector<double> v);
double Varianza(vector<double> v);
void Analisi(int N, int M, vector <double> v, const char* namefile); //data blocking ottenuto dividendo gli M elementi del vettore dei dati in N blocchi. I risultati vengono salvati su file 


#endif
