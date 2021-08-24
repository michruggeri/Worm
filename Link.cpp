#include <cstdlib>
#include <iostream>

#include "Link.hpp"

Link::Link(int dim):_Dim(dim){
  _X0 = new double[Dim()];
  _X1 = new double[Dim()];
  for(int i=0;i<dim;i++){
    _X0[i]=0.0;
    _X1[i]=0.0;
  };
};

Link::Link(int dim,double *x0,double *x1):_Dim(dim){
  _X0 = new double[Dim()];
  _X1 = new double[Dim()];  
  for(int i=0;i<dim;i++){
    _X0[i]=x0[i];
    _X1[i]=x1[i];
  };
};

Link::Link(const Link &model):_Dim(model.Dim()){
  _X0 = new double[Dim()];
  _X1 = new double[Dim()];
  for(int i=0;i<Dim();i++){
    _X0[i]=model.X0(i);
    _X1[i]=model.X1(i);
  };
};

Link::~Link(){
  delete [] _X0;
  delete [] _X1;
};

void Link::SetDim(int d){;
  _Dim = d;
  return;
};

void Link::SetX0(int i,double x){
  _X0[i] = x;
  return;
};

void Link::SetX1(int i,double x){
  _X1[i] = x;
  return;
};

double Link::X0(int i) const{
  return _X0[i];
};

double Link::X1(int i) const{
  return _X1[i];
};

int Link::Dim() const{
    return _Dim;
};

const Link& Link::operator=(const Link &link){
  SetDim(link.Dim());
  delete [] _X0;
  delete [] _X1;
  _X0 = new double[Dim()];
  _X1 = new double[Dim()];
  for(int i=0;i<Dim();i++){
    _X0[i]=link.X0(i);
    _X1[i]=link.X1(i);
  };
  return *this;
};
