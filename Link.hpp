#ifndef LINK__HPP
#define LINK__HPP

class Link{
  private:
    int     _Dim;
    double  *_X0,*_X1;
    
  public:
    Link(int dim=2);
    Link(int dim,double *x0,double *x1);
    Link(const Link &model);
    ~Link();

    void SetDim(int d);
    void SetX0(int i,double x);
    void SetX1(int i,double x);

    double X0(int i)const;
    double X1(int i)const;
    
    int Dim()const;

    const Link& operator=(const Link &link);
};

#endif
