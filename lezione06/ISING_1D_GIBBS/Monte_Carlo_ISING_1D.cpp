#include <iostream>
#include <fstream>
#include <ostream>
#include <cmath>
#include <iomanip>
#include "Monte_Carlo_ISING_1D.h"

using namespace std;

int main(){


  Input(); //Inizialization

  //non è necessaria l'equlibrazione poichè a T=2 dopo poche configurazioni il sistema è all'equlibrio
  while(temp >= 0.5){
   beta = 1./temp;
  cout << "Temperatura: " << temp << endl;
     for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
      accepted = 0;
      Reset(iblk);   //Reset block averages
      for(int istep=1; istep <= nstep; ++istep){

        Move(metro);
        Measure();
        Accumulate(); //Update block averages

      }
      Averages(iblk);   //Print results for current block
    }
    temp -= 0.05;
  }
    ConfFinal(); //Write final configuration

    return 0;
}


void Input(void){

  ifstream ReadInput;

  cout << "Classic 1D Ising model             " << endl;
  cout << "Monte Carlo simulation             " << endl << endl;
  cout << "Nearest neighbour interaction      " << endl << endl;
  cout << "Boltzmann weight exp(- beta * H ), beta = 1/T " << endl << endl;
  cout << "The program uses k_B=1 and mu_B=1 units " << endl;

//Read seed for random numbers
  int p1, p2;
  ifstream Primes("Primes");
  Primes >> p1 >> p2 ;
  Primes.close();

  ifstream input("seed.in");
  input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
  rnd.SetRandom(seed,p1,p2);
  input.close();

//Read input informations
  ReadInput.open("input.dat");

  ReadInput >> temp;
  beta = 1.0/temp;
  cout << "Temperature = " << temp << endl;

  ReadInput >> nspin;
  cout << "Number of spins = " << nspin << endl;

  ReadInput >> J;
  cout << "Exchange interaction = " << J << endl;

  ReadInput >> h;
  cout << "External field = " << h << endl << endl;

  ReadInput >> metro; // if=1 Metropolis else Gibbs

  ReadInput >> nblk;

  ReadInput >> nstep;

  if(metro==1) cout << "The program perform Metropolis moves" << endl;
  else cout << "The program perform Gibbs moves" << endl;
  cout << "Number of blocks = " << nblk << endl;
  cout << "Number of steps in one block = " << nstep << endl << endl;
  ReadInput.close();


//Prepare arrays for measurements
  iu = 0; //Energy
  ic = 1; //Heat capacity
  im = 2; //Magnetization
  ix = 3; //Magnetic susceptibility

  n_props = 4; //Number of observables

//initial configuration                                   //Configurazione iniziale a T=infinito, cioè spin random
  for (int i=0; i<nspin; ++i){

    if(rnd.Rannyu() >= 0.5) s[i] = 1;
    else s[i] = -1;
  }

//Evaluate energy etc. of the initial configuration
  Measure();

//Print initial values for the potential energy and virial
  cout << "Initial energy = " << walker[iu]/(double)nspin << endl;
}


void Move(int metro){
  int o;
  double p, energy_old, energy_new, sm;
  double energy_up, energy_down;


  for(int i=0; i<nspin; ++i){
  //Select randomly a particle (for C++ syntax, 0 <= o <= nspin-1)
    o = (int)(rnd.Rannyu()*nspin);

    if(metro==1){ //Metropolis

      Metropolis(o);
    }
    else{ //Gibbs sampling

      Gibbs(o);
    }
  }
}

double Boltzmann(int sm, int ip){
  double ene = -J * sm * ( s[Pbc(ip-1)] + s[Pbc(ip+1)] ) - h * sm;
  return ene;
}

void Measure(){
  int bin;
  double u = 0.0, m = 0.0, chi = 0.0, s_tot = 0.0, M = 0.0,u_2 = 0.0;

//cycle over spins
  for (int i=0; i<nspin; ++i){
    u += -J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]);   //Energia interna (H)
    //u_2 += pow(-J * s[i] * s[Pbc(i+1)] - 0.5 * h * (s[i] + s[Pbc(i+1)]),2);    // Energia al quadrato
    s_tot += s[i];

  }
  chi = pow(s_tot,2);   //Suscettività
  M = s_tot;    //Magnetizzazione

  walker[iu] = u;
  walker[ic] = pow(u,2);
  walker[im] = M;
  walker[ix] = chi;
  S = M;
}


void Reset(int iblk){ //Reset block averages

  if(iblk == 1){
    for(int i=0; i<n_props; ++i){
      glob_av[i] = 0;
      glob_av2[i] = 0;
    }
  }

  for(int i=0; i<n_props; ++i){
    blk_av[i] = 0;
  }
  blk_norm = 0;
  attempted = 0;
  accepted = 0;
}


void Accumulate(void){ //Update block averages

  for(int i=0; i<n_props; ++i){
    blk_av[i] = blk_av[i] + walker[i];
  }
  blk_norm = blk_norm + 1.0;
}


void Averages(int iblk){ //Print results for current block

  ofstream Ene, Heat, Mag, Chi, ResultsU, ResultsC, ResultsM, ResultsX;
  const int wd=12;

//  cout << "Block number " << iblk << endl;
//  cout << "Acceptance rate " << accepted/attempted << endl << endl;

  
  Ene.open("output.ene.0",ios::app);
  stima_u = blk_av[iu]/blk_norm/(double)nspin; //Energy
  glob_av[iu]  += stima_u;
  glob_av2[iu] += stima_u*stima_u;  
  err_u=Error(glob_av[iu],glob_av2[iu],iblk);
  Ene << setw(wd) << iblk <<  setw(wd) << stima_u << setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;
  Ene.close();

  Heat.open("output.heat.0",ios::app);
 stima_u2 =blk_av[ic]/blk_norm/(double)nspin; // <H^2>
  stima_c = pow(beta,2) * (stima_u2*(double)nspin-pow(stima_u,2)*pow((double)nspin,2) )/ (double)nspin; //Colore specifico in unità naturali
  glob_av[ic]  += stima_c;
  glob_av2[ic] += stima_c*stima_c;
  err_c=Error(glob_av[ic],glob_av2[ic],iblk);
  Heat << setw(wd) << iblk <<  setw(wd) << stima_c << setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;
  Heat.close();

  Mag.open("output.mag.0",ios::app);
  stima_m = blk_av[im]/blk_norm/(double)nspin; //Magnetizzazione
  glob_av[im]  += stima_m;
  glob_av2[im] += stima_m*stima_m;
  err_m=Error(glob_av[im],glob_av2[im],iblk);
  Mag << setw(wd) << iblk <<  setw(wd) << stima_m << setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;
  Mag.close();

  Chi.open("output.chi.0",ios::app);
  stima_x = beta*blk_av[ix]/blk_norm/(double)nspin; //Suscettività
  glob_av[ix]  += stima_x;
  glob_av2[ix] += stima_x*stima_x;
  err_x=Error(glob_av[ix],glob_av2[ix],iblk);
  Chi << setw(wd) << iblk <<  setw(wd) << stima_x << setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
  Chi.close();

  if(iblk == nblk){
    ResultsU.open("../ene.result.Gibbs",ios::app);
    ResultsU << setw(wd) << temp <<  setw(wd) << glob_av[iu]/(double)iblk << setw(wd) << err_u << endl;

    ResultsC.open("../heat.result.Gibbs",ios::app);
    ResultsC << setw(wd) << temp <<  setw(wd) << glob_av[ic]/(double)iblk << setw(wd) << err_c << endl;

    ResultsM.open("../mag.result.Gibbs",ios::app);
    ResultsM << setw(wd) << temp <<  setw(wd) << glob_av[im]/(double)iblk << setw(wd) << err_m << endl;

    ResultsX.open("../chi.result.Gibbs",ios::app);
    ResultsX << setw(wd) << temp <<  setw(wd) << glob_av[ix]/(double)iblk << setw(wd) << err_x << endl;
  }
//  cout << "----------------------------" << endl << endl;
}


void ConfFinal(void){
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");
  for (int i=0; i<nspin; ++i){
    WriteConf << s[i] << endl;
  }
  WriteConf.close();

  rnd.SaveSeed();
}

int Pbc(int i){  //Algorithm for periodic boundary conditions

  if(i >= nspin) i = i - nspin;
  else if(i < 0) i = i + nspin;
  return i;
}

double Error(double sum, double sum2, int iblk){
  return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}

void Metropolis(int o){

  attempted++;
  double deltaE = 2.* J * s[o] * ( s[Pbc(o-1)] + s[Pbc(o+1)] ) +h*s[o]*2;
//  cout << deltaE << endl;
  if(deltaE < 0) s[o] = s[o]*(-1), accepted++;
  else{
    double a = rnd.Rannyu();
    double boltzweight = exp(-1/temp * deltaE);
//    cout << endl << "a: " << a << endl;
		if( a <= boltzweight){
//			cout << "False poi True" << endl;
			s[o] = s[o]*(-1);
      accepted++;
    }
  }

}

void Gibbs(int o){

 double deltaE = 2.* J * ( s[Pbc(o-1)] + s[Pbc(o+1)] ) + 2 * h;
 double p_up = 1./(1.+exp(-deltaE*beta));
 double a = rnd.Rannyu(0,1.);
   
   if(a <= p_up){
	s[o]=1;
   }
   else{
	s[o]=-1;
   }


}


