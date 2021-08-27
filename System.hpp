#ifndef SYSTEM__HPP
#define SYSTEM__HPP

#include <fstream>
#include <iostream>
#include <cstring>
#include <iomanip>
#include <cmath>


#include "Link.hpp"
#include "Random.hpp"
#include "Potentials.hpp"

class System{
  private:
    int SLICES;
    int MAXPART;
    int SECTOR;
    int MBAR;
    int MBRI;
    int BLOCKS;
    int EQUI;
    int SPARSE;
    int STEPS;
    int RESTART;
    int NMOV;
    int WLINK;
    int DIM;
    int NG;
    int NK;
    int NX;
    int NY;
    int NUM;
    int TMAX;

    int *NLNK;
    int *FIRST_EMPTY;
    int *NEXT_TAU;
    int *PREV_TAU;
    int *IRA;
    int *MASHA;
    
    int **INDEX;
    int **NEXT_EMPTY;
    int **PREV;
    int **NEXT;

    double TEMP;
    double TAU;
    double CI0;
    double MASS;
    double LAMBDA;
    double VOL;
    double SIGMA;
    double RMAX;
    double DR;
    double DX1;
    double DX2;
    double DK;
    double BEFF;
    double ECUT;

    double *SIDE;
    double *DX;
    double *TRY;
    double *ACC;
    double *WEIGHT;
    double *COST;

    double **KINETIC;
    double **POTENTIAL;
    double **ENERGY;
    double **SUPER;

    double ***GOFR;
    
    double ****RHOXY;
    double ****GOFXY;
    double ****SKX_RE;
    double ****SKX_IM;

    Link** LINKS;
    Random RNG;

  public:
    System(char* input_name);
    ~System();

    void Read_Input(char* filename);
    void Restart(int res);

    void Set_Slices(int slices);
    void SetMASSPart(int part);
    void SetSECTOR(int sect);
    void Set_T(double t);
    void SetTAU(double tau);
    void SetMBAR(int mbar);
    void SetMBRI(int mbri);
    void SetC0(double c);
    void SetBLOCKS(int blocks);
    void SetEQUI(int equi);
    void SetSPARSE(int sparse);
    void SetSTEPS(int steps);
    void SetMASS(double m);
    void SetVOL(double v);
    void Set_L(int i,double l);

    int Sector() const;
    int Slices() const;
    int MPart()  const;
    int MBar()   const;
    int MBri()   const;
    int Blocks() const;
    int Equi()   const;
    int Sparse() const;
    int Steps()  const;
    int WLink()  const;
    int NP()     const;
    int NLnk(int i)  const;
    int Roulette(double vect[],int dim,double sig) const;

    double T()   const;
    double Tau() const;
    double C0()  const;
    double M()   const;
    double Vol() const;
    double L(int i) const;
    double Sigma()  const;
    double V(int time,int part,double x0[],double x1[]) const;
    double V(int time,int part,Link l) const;
    double Vpot(int time,int part,double x0[]) const;
    double Vpot(int time,int part,Link l) const;
    double Pot(Link l1, Link l2) const;
    double (*Uint)(double x[],double y[],int dim,double l[],double ecut,double cost[]);
    double (*Uext)(double x[],int dim);
    double U0(double r) const;
    double Gauss(double d2,double sigma) const;

    typedef double (System::*GetKin)(int b);
    double GetKin_0(int b);
    double GetKin_SuperW(int b);
    double GetKin_SuperA(int b);
    double GetPot() const;

    Link& Links(int i,int j) const;

    void PrintSector() const;
    void PrintWorm()   const;
    void PrintSyst()   const;
    void PrintAcc()   const;
    void PrintPart(int p) const;
    void WriteOut(char filename[],char empties[]) const;

    void AddLink(int time,int index,double *x0,double *x1);
    void RemLink(int time,int index);
    void Bridge(double x0[],double xf[],double xnew[][3],int emme);
    
    void Advance();
    void Recede();
    void Open();
    void Close();
    void Swap();
    void Wiggle();
    void Shift();

    void Move();
    void Measure(int block);
    void Average(int block);
    void Run();
    void TestPot() const;
    void GetGofR(int block);
    void NormGofR(int block);
    void GetSofK(int block);
    void GetRho(int block);
    void GetGxy(int block);
    void PrintG(int block) const;
    void PrintS(int block) const;
    void PrintRho(int block) const;
    void PrintGxy(int block) const;
};

    extern System::GetKin aaa;
#endif
