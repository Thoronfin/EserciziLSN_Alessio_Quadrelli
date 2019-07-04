#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <iomanip>
#include "MolDyn_NVE.h"
#include "random.h"
#include <vector>
using namespace std;

int main(){ 
  Input();   //Inizialization
      
  //equlilibrazione del sistema  
  int counter = 1;
  int Nrestart=0;
  ofstream Term;
  Term.open("output_restart.dat",ios::app);
  
  cout<<"La temperatura di target è: "<<T_target<<endl;
  
 do{

    Move();
    counter++;
  
    if(counter%10 == 0){
	Restart();
	Nrestart++;
	Term << stima_temp << endl;
    }
  }
  while(counter < 10000);
  
  Term.close();
  cout<<"Le velocità del sistema sono state riscalate: "<<Nrestart<<" volte"<<endl;
  cout<<"Temperatura sistema ha raggiunto il valore: "<<temp<<endl;


 //simulazione del sistema all'equlibrio
  int nconf = 1;
  for(int iblk=1; iblk <= nblk; ++iblk){ //Simulation
  cout << "BLOCCO NUMERO: " << iblk << endl << endl; 
 
    Reset(iblk);   //Reset block averages
    for(int istep=1; istep <= nstep; ++istep)
    {
      Move();
      
      if(nstep%10 == 0){
	Measure();
        Accumulate(); //Update block averages
        //ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
      }
    }
    Averages(iblk);   //Print results for current block
  }
  //ConfFinal(); //Write final configuration

  return 0;
}


 

void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
 
  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  //inizializzazione seme del generatore di numeri casuali
  int seed[4];
	int p1, p2;
	ifstream Primes("Primes");
   	if (Primes.is_open()){
		Primes >> p1 >> p2 ;
	} else cerr << "PROBLEM: Unable to open Primes" << endl;
	Primes.close();
	ifstream input("seed.in");
	string property;
	if (input.is_open()){
		while ( !input.eof() ){
	        input >> property;
        if( property == "RANDOMSEED" ){
            	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            	rnd.SetRandom(seed,p1,p2);
        }
        }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nblk;
  ReadInput >> nstep;

   //Tail corrections for potential energy and pressure
  vtail = (8.0*pi*rho)/(9.0*pow(rcut,9)) - (8.0*pi*rho)/(3.0*pow(rcut,3));
  ptail = (32.0*pi*rho)/(9.0*pow(rcut,9)) - (16.0*pi*rho)/(3.0*pow(rcut,3));
  cout << "Tail correction for the potential energy = " << vtail << endl;
  cout << "Pail correction for the virial           = " << ptail << endl; 

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of blocks = " << nblk << endl << endl;
  cout << "Number of steps in each block = " << nstep << endl << endl;
  ReadInput.close();
  T_target = temp;

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  iw = 4; //pressure
  iv_0=5; //potenziale senza correzioni di coda
  ie_0=6; //energia totale senza correzioni di coda
  iw_0=7; //viriale senza correzioni di coda
  n_props = 8; //Number of observables

//measurement of g(r)
  igofr = n_props;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rnd.Rannyu(0.,1.) - 0.5;
     vy[i] = rnd.Rannyu(0.,1.) - 0.5;
     vz[i] = rnd.Rannyu(0.,1.) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3. * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = x[i] - vx[i] * delta;
     yold[i] = y[i] - vy[i] * delta;
     zold[i] = z[i] - vz[i] * delta;

    
   }

   
   return;
}



void Restart(void) { //riscalamento progressivo della velocità partendo da T fornita nel file di input fino alla temperatura di equilibrio
	
    double sumv[3] = {0.0, 0.0, 0.0};

   for(int i=0; i<npart; ++i){ //Verlet integration scheme

     vx[i] = Pbc(x[i] - xold2[i])/(2.0 * delta);
     vy[i] = Pbc(y[i] - yold2[i])/(2.0 * delta);
     vz[i] = Pbc(z[i] - zold2[i])/(2.0 * delta);
   }

   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
    for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   Square_Velocity();
   temp = (1.0/ 3.0) * t;
   stima_temp=temp; 

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;
   }
     return;
}

void Velocity(void){ //calculates velocities of the particles

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    vx[i] = Pbc(x[i] - xold2[i])/(2.0 * delta);
    vy[i] = Pbc(y[i] - yold2[i])/(2.0 * delta);
    vz[i] = Pbc(z[i] - zold2[i])/(2.0 * delta);
  }

}


void Square_Velocity(void) { //calcolo delle somma delle velocità quadratiche delle particelle  

	
    t=0.; // somma delle velocità quadratiche medie
    Velocity(); //aggiornamento delle velocità

   for (int i=0; i<npart; ++i){
     t += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
 
   t=t/double(npart);
   return;
}

void Temperature(void) {
  ofstream Temp;
  Temp.open("output_temp.dat",ios::app);
  stima_temp = (1. / 3.) * t; //Temperature= (1/3)*( <v^2>/N)
  walker[it]=stima_temp;
  Temp << stima_temp << endl;
  Temp.close();
  return;
}

void Kinetic_Energy(void) {
  ofstream Ekin;
  Ekin.open("output_ekin.dat",ios::app);
  stima_kin = t/2.; //Kinetic energy
  Ekin << stima_kin  << endl;
  Ekin.close();
  walker[ik]=stima_kin;
  return;
}

void Total_Energy(void) {
  ofstream Etot;
  Etot.open("output_etot.dat",ios::app);
  stima_etot = stima_kin+stima_pot; //Total energy
  walker[ie]=stima_etot;
  Etot << stima_etot << endl;
  Etot.close();
  return;
}




void Termini_a_due_corpi(void){ //calcolo delle grandezze a due corpi: Energia potenziale, pressione e g(r)

  ofstream Epot;
  Epot.open("output_epot.dat",ios::app);
  ofstream Pres;
  Pres.open("output_pres.dat",ios::app);

  
  v = 0.;
  virial=0;
  double dx, dy, dz, dr;
  double vij;

  //reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k) walker[k]=0.0;

  
  //cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( x[i] - x[j] );
     dy = Pbc( y[i] - y[j] );
     dz = Pbc( z[i] - z[j] );

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     //update of the histogram of g(r)
     walker[igofr+int(dr/bin_size)]+=2;

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
       v += vij;

       vij = 48.*(1./pow(dr,12) - 0.5/pow(dr,6));
       virial += vij;
     }
    }          
  }


  stima_pot = v/(double)npart + vtail; //Potential energy
  stima_pot = v/(double)npart;
  Epot << stima_pot  << endl;
  Epot.close();

  stima_pres = rho*stima_temp+( ((1./3.)*virial) + ptail*(double(npart) ))/vol;
  stima_pres = rho*stima_temp+ virial/(3.*vol);
  Pres << stima_pres << endl;
  Pres.close();

  walker[iv] = v/(double)npart;
  walker[iw] = virial/3.;
  return;
}




void Measure(){ //Properties measurement
  
  Square_Velocity();
  Termini_a_due_corpi();
  Kinetic_Energy();
  Temperature();
  Total_Energy();
 
  return;
}


void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme

    //r(t-dt) lo salvo per il calcolo delle velocità nella funzione Velocity
    xold2[i] = xold[i];
    yold2[i] = yold[i];
    zold2[i] = zold[i];

    //r(t+dt)
    xnew = Pbc( 2.0 * x[i] - xold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - yold[i] + fy[i] * pow(delta,2) );
    znew = Pbc( 2.0 * z[i] - zold[i] + fz[i] * pow(delta,2) );

  
    //r(t)
    xold[i] = x[i];
    yold[i] = y[i];
    zold[i] = z[i];

    x[i] = xnew;   //r(t+dt)
    y[i] = ynew;
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}



void Reset(int iblk) //Reset block averages
{
   
   if(iblk == 1)
   {
       for(int i=0; i<n_props; ++i)
       {
           glob_av[i] = 0;
           glob_av2[i] = 0;
       }
   }

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = 0;
   }
   blk_norm = 0;

}

void Accumulate(void) //Update block averages
{

   for(int i=0; i<n_props; ++i)
   {
     blk_av[i] = blk_av[i] + walker[i];
   }
   blk_norm = blk_norm + 1.0;
}


void Averages(int iblk) //Print results for current block
{
    
   double r, gdir;
   ofstream Gofr, Gave, Epot, Ekin, Etot, Pres, Temp,Epot_0,Pres_0, Etot_0 ;
   double stima_epot_0,stima_pres_0;
   const int wd=20;
    
    cout << "Block number " << iblk << endl;
    
    
    Epot.open("../ave_epot.0",ios::app);
    Ekin.open("../ave_ekin.0",ios::app);
    Etot.open("../ave_etot.0",ios::app);
    Pres.open("../ave_pres.0",ios::app);
    Temp.open("../ave_temp.0",ios::app);
    Gofr.open("../ave_gofr.0",ios::app);
    Gave.open("../ave_gave.0",ios::app);

    Epot_0.open("../ave_epot_0.0", ios::app);
    Etot_0.open("../ave_etot_0.0", ios::app);
    Pres_0.open("../ave_pres_0.0", ios::app);
    
    //sistema con correzioni di coda
   
    stima_pot = blk_av[iv]/blk_norm + vtail; //Potential energy
    glob_av[iv] += stima_pot;
    glob_av2[iv] += stima_pot*stima_pot;
    err_epot=Error(glob_av[iv],glob_av2[iv],iblk);

    stima_kin = blk_av[ik]/blk_norm; //Kinetic energy
    glob_av[ik] += stima_kin;
    glob_av2[ik] += stima_kin*stima_kin;
    err_ekin=Error(glob_av[ik],glob_av2[ik],iblk);

    stima_etot = blk_av[ie]/blk_norm; //Total energy
    glob_av[ie] += stima_etot;
    glob_av2[ie] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie],glob_av2[ie],iblk);


    stima_temp = blk_av[it]/blk_norm; //Potential energy
    glob_av[it] += stima_temp;
    glob_av2[it] += stima_temp*stima_temp;
    err_temp=Error(glob_av[it],glob_av2[it],iblk);

    
    stima_pres = rho * temp + (blk_av[iw]/blk_norm + ptail * (double)npart) / vol; //Pressure
    glob_av[iw] += stima_pres;
    glob_av2[iw] += stima_pres*stima_pres;
    err_pres=Error(glob_av[iw],glob_av2[iw],iblk);

    

  

//Potential energy per particle
    Epot << setw(wd) << iblk <<  setw(wd) << stima_pot << setw(wd) << glob_av[iv]/(double)iblk << setw(wd) << err_epot << endl;
//Kinetic energy per particle
    Ekin << setw(wd) << iblk <<  setw(wd) << stima_kin << setw(wd) << glob_av[ik]/(double)iblk << setw(wd) << err_ekin << endl;
//Total energy per particle
    Etot << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie]/(double)iblk << setw(wd) << err_etot << endl;
//Temperature
    Temp << setw(wd) << iblk <<  setw(wd) << stima_temp << setw(wd) << glob_av[it]/(double)iblk << setw(wd) << err_temp << endl;
//Pressure
    Pres << setw(wd) << iblk <<  setw(wd) << stima_pres << setw(wd) << glob_av[iw]/(double)iblk << setw(wd) << err_pres << endl;


   //sistema senza correzioni di cora
   
    stima_epot_0 = blk_av[iv]/blk_norm ; //Potential energy
    glob_av[iv_0] += stima_epot_0;
    glob_av2[iv_0] += stima_epot_0*stima_epot_0;
    err_epot=Error(glob_av[iv_0],glob_av2[iv_0],iblk);


    stima_etot = blk_av[ie]/blk_norm + vtail; //Total energy
    glob_av[ie_0] += stima_etot;
    glob_av2[ie_0] += stima_etot*stima_etot;
    err_etot=Error(glob_av[ie_0],glob_av2[ie_0],iblk);


    stima_pres_0 = rho * temp + (blk_av[iw]/blk_norm)/ vol; //Pressure
    glob_av[iw_0] += stima_pres_0;
    glob_av2[iw_0] += stima_pres_0*stima_pres_0;
    err_pres=Error(glob_av[iw_0],glob_av2[iw_0],iblk);

//Potential energy per particle
    Epot_0 << setw(wd) << iblk <<  setw(wd) << stima_epot_0 << setw(wd) << glob_av[iv_0]/(double)iblk << setw(wd) << err_epot << endl;
//Total energy per particle
    Etot_0 << setprecision(10) << setw(wd) << iblk <<  setw(wd) << stima_etot << setw(wd) << glob_av[ie_0]/(double)iblk << setw(wd) << err_etot << endl;
//Pressure
    Pres_0 << setw(wd) << iblk <<  setw(wd) << stima_pres_0 << setw(wd) << glob_av[iw_0]/(double)iblk << setw(wd) << err_pres << endl;


	//g(r)  

    for (int k=igofr; k<igofr+nbins; ++k) {
	r = (k-igofr)*bin_size; //faccio variare r
	deltaV = (4./3.)*M_PI*(pow(r+bin_size,3)-pow(r,3));  //calcolo il delta(r)
	stima_gdir = (1./(rho*deltaV*(double)npart))*(blk_av[k]/blk_norm);   
	glob_av[k] += stima_gdir;
    	glob_av2[k] += stima_gdir*stima_gdir;
    	err_gdir=Error(glob_av[k],glob_av2[k],iblk);
	Gofr    << setw(wd) << r <<  setw(wd) << stima_gdir << endl;  //stampo i valori di g(r) di ciascun blocco
	//cout << iblk << endl;
	if(iblk == nblk){  //prendo il valore del blocco finale
   	 
    	  Gave << setw(wd) << r <<  setw(wd) << stima_gdir << setw(wd) << err_gdir<< endl; //stampo l'ultimo blocco di g(r) al variare del raggio con rispettivo errore
		//cout << "gave scritto" << endl;
			
	}	

    }
    Gofr<<endl;
    cout << "----------------------------" << endl << endl;

    Epot.close();
    Ekin.close();
    Etot.close();
    Temp.close();
    Pres.close();
    Gofr.close();
    Gave.close();

    Epot_0.close();
    Etot_0.close();
    Pres_0.close();
}

double Error(double sum, double sum2, int iblk)
{
    return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)iblk);
}



void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;

  cout << "Print final configuration to file config.final " << endl << endl;
  WriteConf.open("config.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();
  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
}

