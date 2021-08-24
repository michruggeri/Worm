#include <iostream>
#include <cstdlib>
#include <fstream>
#include <cmath>

#include "ran2.hpp"
#include "Random.hpp"

double PG = acos(-1.0);

Random::Random(long seed):_seed(seed){};

Random::Random(const Random &model):_seed(model._seed){};

void Random::SetSeed(long newseed) {
  _seed = newseed;
};

int Random::Int_Random(int n=1, int m=0) const {
  if(n<m){
    int s;
    s = n;
    n = m;
    m = s;
  };
  return int(m+ran2(Seed())*(n-m));
};

double Random::Double_Random(double n=1.0, double m=0.0) const {
  if(n<m){
    double s;
    s = n;
    n = m;
    m = s;
  };
  return m+ran2(Seed())*(n-m);
};

double Random::Gauss_Random(double av=0.0,double sigma=1.0) const {

  double A = Double_Random(0.0,1.0);
  double B = Double_Random(0.0,1.0);
  return sqrt(-2.0*log(A))*cos(2.0*PG*B)*sigma+av;
};

const Random & Random::operator=(const Random &random) {
  _seed = random._seed;
  return *this;
};
