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
class Psi: public FunzioneBase {

	private:
		double _mu;
		double _sigma;
	public:

		Psi(double mu, double sigma);
		Psi();
		~Psi();	
		virtual double Eval (vector<double> x) const;
		void SetMu( double mu){_mu=mu;}
		double GetMu(){ return _mu;}
		void SetSigma(double sigma){_sigma=sigma;}
		double GetSigma(){ return _sigma;}

};

//hamiltoniana del sistema
class Hamiltoniana: public FunzioneBase {

	private:
		Psi *_psi;

	public:
		Hamiltoniana(Psi * psi){_psi=psi;}
		~Hamiltoniana(){ }	
		virtual double Eval (vector<double> x) const;


};




#endif

