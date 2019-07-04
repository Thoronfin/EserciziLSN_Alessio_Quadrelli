#ifndef _Statistica_h_
#define _Statistica_h_

#include <vector>

using namespace std;

double Somma(vector<double> v, int i, int j); //somma gli elementi di v dall'i-esimo incluso, al j-esimo escluso
double Media(vector<double> v);
double Varianza(vector<double> v);
void Analisi(int N, int M, vector <double> v, const char* namefile); // esegue il data blocking su un vettore di M dati composta da N sequenze di risultati del random walk con il medesimo passo, formando blocchi da M/N^2 dati e restituendo solo la media e la deviazione standard della media con tutti i blocchi disponibili.


#endif
