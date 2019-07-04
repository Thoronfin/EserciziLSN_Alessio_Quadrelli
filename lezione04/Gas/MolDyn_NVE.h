#ifndef __MolDyn_NVE__
#define __MolDyn_NVE__


#include "random.h"
//parameters, observables
const int m_props=200;
double walker[m_props];
int n_props, nblk;
int iv,ik,it,ie,iw,igofr,iv_0,ie_0,iw_0;
double stima_pot, stima_kin, stima_etot, stima_temp, stima_pres;
double virial;
double T_target;

//generator
Random rnd;

// averages
double acc,att;
double blk_av[m_props],blk_norm;
double glob_av[m_props],glob_av2[m_props];
double err_temp, err_ekin, err_epot, err_etot,err_pres,err_gdir;

//configuration
const int m_part=108;
double x[m_part],y[m_part],z[m_part];
double xold[m_part],yold[m_part],zold[m_part];
double xold2[m_part],yold2[m_part],zold2[m_part];
double vx[m_part],vy[m_part],vz[m_part];
double t; //variabile per il calcolo della velocità quadratica media del sistema
double v; //variabile dell'energia potenziale--->utile per la stima dell'energia totale

//g(r)
int nbins;
double bin_size;

// thermodynamical state
int npart;
double energy,temp,vol,rho,box,rcut;

// simulation
int nstep, iprint, seed;
double delta;

//pigreco
const double pi=3.1415927;
double vtail,ptail;
double deltaV, stima_gdir;

//functions
void Input(void);
void Move(void); //Algoritmo di Verlet
void Restart(void); // Algoritmo per la termalizzazione del sistema
void ConfFinal(void);
void ConfXYZ(int);
void Measure(void); //misura delle grandezze del sistema: energia cinetica, energia potenizale, energia totale, temperatura, pressione
double Force(int, int);
double Pbc(double);
void Termini_a_due_corpi(void); //calcolo delle grendezze che sono il risultato di interazioni a due corpi: energia poteniale e pressione
void Square_Velocity(void); //calcolo delle somma delle velocità quadratiche delle particelle
void Kinetic_Energy(void); //calcolo dell'energia cinetica
void Temperature(void); //calcolo della temperatura
void Total_Energy(void); //calcolo dell'energia totlae
void Reset(int);
void Accumulate(void);
void Averages(int);
double Error(double, double, int);

#endif // __MolDyn_NVE__


