#include <cmath>
#include "Potentials.hpp"

double Pi   = acos(-1.0);
double Hbar = 1.0;

double (*nint)(double n);

double abc(double n){
  return 0.0;
};

double pbc(double n){
  return floor(n+0.5);
};

// --- POTENZIALI DI INTERAZIONE --------------------------------- //

// PARTICELLE NON INTERAGENTI (default)
double NonInt(double x[],double y[],int dim, double l[], double ecut,double cost[]){
  return 0.0;
};

// LENNARD - JONES  ----- cost[0] = profondita' della buca, cost[1] = raggio hard core
double LenJ(double x[],double y[],int dim, double l[], double ecut,double cost[]){
  double r2=0.0;
  double u;
  double r;
  for(int i=0;i<dim;i++){
    r = x[i]-y[i];
    r = r-l[i]*nint(r/l[i]);
    r2 = r2+r*r;
  };
  u = 4.0*cost[0]*(pow(cost[1]*cost[1]/r2,6.0)-pow(cost[1]*cost[1]/r2,3.0))-ecut;
  return u;
};

// DIPOLI VERTICALI ----- cost[0] = forza dei dipoli
double VDip(double x[],double y[],int dim,double l[],double ecut,double cost[]){
  double r2=0.0;
  double u;
  double r;
  for(int i=0;i<dim;i++){
    r = x[i]-y[i];
    r = r-l[i]*nint(r/l[i]);
    r2 = r2+r*r;
  };
  u = cost[0]*cost[0]*pow(r2,-1.5);
  return u;
};

// DIPOLI GENERICI ----- cost[0] = forza dei dipoli, cost[1] = cos(\Theta)
double Dip(double x[],double y[],int dim,double l[],double ecut,double cost[]){
  double u;
  double X,Y,r;
  X = (x[0]-y[0]);
  X = X-l[0]*nint(X/l[0]);
  Y = (x[1]-y[1]);
  Y = Y-l[1]*nint(Y/l[1]);
  r = sqrt(X*X+Y*Y);
  if(r==0) return 1.0e10;
  double cosen = fabs(X)*cost[1]/r;
  u = cost[0]*cost[0]*(1.0 - 3.0*cosen*cosen)/pow(r,3.0);
  return u;
};

// --- POTENZIALI ESTERNI ------------------------------------- //

// NESSUN CAMPO ESTERNO (default)
double Zero(double x[],int dim){
  return 0.0;
};

// TRAPPOLA ARMONICA - m = 1.0, \omega = 1.0 (niente pbc!)
double Harm(double x[],int dim){
  double r2 = 0.0;
  double u;
  for(int i=0;i<dim;i++)
    r2 = r2 +x[i]*x[i];
  u = 0.5*r2;
  return u;
};
