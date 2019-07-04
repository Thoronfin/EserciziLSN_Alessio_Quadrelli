#ifndef _Funzioni_h_
#define _Funzioni_h_
#include <vector>
using namespace std;

//classe astratta contente le distribuzioni da campionare
class FunzioneBase {
	public:
		virtual double Eval(vector<double> x) const=0;

};


// Distribuzioni da campionare
class uno_s: public FunzioneBase {

	virtual double Eval (vector<double> x) const;

};

class due_p: public FunzioneBase {

	virtual double Eval (vector<double> x) const;

};



#endif

