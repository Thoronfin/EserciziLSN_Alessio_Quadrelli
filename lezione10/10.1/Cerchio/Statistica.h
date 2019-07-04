#ifndef _Statistica_h_
#define _Statistica_h_

#include <vector>

using namespace std;

double Somma(vector<double> v, int i, int j); //somma gli elementi di v dall'i-esimo incluso, al j-esimo escluso
double Media(vector<double> v);
double Varianza(vector<double> v);
void Analisi(int N, int M, vector <double> v, const char* namefile);
void AnalisiMatrici(int N, vector<vector <double>> v, double fattore_scala, double traslazione_inidice, const char* namefile); //data blocking di una matrice stampando media e deviazione standard della media realizzata con tutti i blocchi


#endif
