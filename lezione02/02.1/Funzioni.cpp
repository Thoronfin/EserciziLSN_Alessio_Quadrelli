#include <cmath>
#include "Funzioni.h"

using namespace std;

// funzione da integrae in esercizio 02.1
double Coseno:: Eval(double x) const {
	return (M_PI/2.)*cos(x*(M_PI/2.));

}
// funzione da integrae in esercizio 02.1 divisa per 2*(1-x)
double G:: Eval(double x) const {
	return ((M_PI/2.)*cos(x*(M_PI/2))) / (2.*(1.-x)) ;

}
