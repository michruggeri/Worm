#ifndef POTENTIALS__CPP
#define POTENTIALS_CPP

extern double Pi;
extern double Hbar;
extern double (*nint)(double n);

double abc(double n);
double pbc(double n);

double NonInt(double x[],double y[],int dim, double l[], double ecut,double cost[]);
double LenJ(double x[],double y[],int dim, double l[], double ecut,double cost[]);
double VDip(double x[],double y[],int dim,double l[],double ecut,double cost[]);
double Dip(double x[],double y[],int dim,double l[],double ecut,double cost[]);

double Zero(double x[],int dim);
double Harm(double x[],int dim);
double Cigar(double x[],int dim);

#endif
