#include <fstream>
#include <iostream>
#include <cstring>
#include <iomanip>
#include <cmath>

#include "System.hpp"

System::GetKin KinAndSup;

System::System(char* input_name){
  COST        = new double[5];
  NK = 0;
  Uext = Zero;
  Uint = NonInt;
  nint = abc;
  KinAndSup = &System::GetKin_0;
  std::cout << "\nWorm canonical simulation, algorithm based on Boninsegni & Prokof'ef, Phys. Rev. E 74, 036701 (2006)\n\n";
  Read_Input(input_name);
  std::cout << "Initialializing arrays... ";
  LINKS       = new Link*[SLICES];
  NLNK        = new int [SLICES];
  NEXT_TAU    = new int [SLICES];
  PREV_TAU    = new int [SLICES];
  FIRST_EMPTY = new int [SLICES];
  NEXT_EMPTY  = new int*[SLICES];
  INDEX       = new int*[SLICES];
  NEXT        = new int*[SLICES];
  PREV        = new int*[SLICES];
  KINETIC     = new double*[3];
  POTENTIAL   = new double*[3];
  ENERGY      = new double*[3];
  SUPER       = new double*[3];
  GOFR        = new double**[3];
  SKX_RE      = new double***[3];
  SKX_IM      = new double***[3];
  RHOXY       = new double***[3];
  GOFXY       = new double***[3];
  BEFF = 0.0;
  NX = 1+2*NG;
  NY = 1+2*NG;
  DX1 = L(0)/double(NX);
  DX2 = L(1)/double(NY);
  for(int i=0;i<3;i++){
    KINETIC[i]   = new double[BLOCKS];
    POTENTIAL[i] = new double[BLOCKS];
    ENERGY[i]    = new double[BLOCKS];
    SUPER[i]     = new double[BLOCKS];
    GOFR[i]      = new double*[BLOCKS];
    SKX_RE[i]    = new double**[BLOCKS];
    SKX_IM[i]    = new double**[BLOCKS];
    RHOXY[i]     = new double**[BLOCKS];
    GOFXY[i]     = new double**[BLOCKS];
    for(int j=0;j<BLOCKS;j++){
      KINETIC[i][j]   = 0.0;
      POTENTIAL[i][j] = 0.0;
      ENERGY[i][j]    = 0.0;
      SUPER[i][j]     = 0.0;
      GOFR[i][j] = new double[NG];
      for(int k=0;k<NG;k++)
        GOFR[i][j][k] = 0.0;
      SKX_RE[i][j] = new double*[NK];
      SKX_IM[i][j] = new double*[NK];
      for(int k=0;k<NK;k++){
        SKX_RE[i][j][k] = new double[TMAX];
        SKX_IM[i][j][k] = new double[TMAX];
        for(int t=0;t<TMAX;t++){
          SKX_RE[i][j][k][t] = 0.0;
          SKX_IM[i][j][k][t] = 0.0;
        };
      };
      RHOXY[i][j] = new double*[NX];
      GOFXY[i][j] = new double*[NX];
      for(int k=0;k<NX;k++){
        RHOXY[i][j][k] = new double[NY];
        GOFXY[i][j][k] = new double[NY];
        for(int l=0;l<NY;l++){
          RHOXY[i][j][k][l] = 0.0;
          GOFXY[i][j][k][l] = 0.0;
        };
      };
    };
  };
  WEIGHT = new double[BLOCKS];
  for(int i=0;i<BLOCKS;i++)
    WEIGHT[i] = 0.0;
  for(int i=0;i<SLICES;i++){
    NEXT_TAU[i]=(i+1)%Slices();
    PREV_TAU[i]=(i-1+Slices())%Slices();
    NLNK[i]  = 0;
    LINKS[i]      = new Link[MAXPART];
    NEXT_EMPTY[i] = new int[MAXPART];
    INDEX[i]      = new int[MAXPART];
    NEXT [i]      = new int[MAXPART];
    PREV [i]      = new int[MAXPART];
    for(int j=0;j<MAXPART;j++)
      LINKS[i][j] = Link(DIM);
  };
  NMOV = 7;
  TRY = new double[NMOV];
  ACC = new double[NMOV];
  for(int i=0;i<NMOV;i++){
    TRY[i] = 0.0;
    ACC[i] = 0.0;
  };
  IRA   = new int[2];
  MASHA = new int[2];
  for(int i=0;i<2;i++){
    IRA[i]   = -1;
    MASHA[i] = -1;
  };
  WLINK = 0;
  std::cout << "Done!\n";
  Restart(RESTART);
  PrintSyst();
  std::cout << std::scientific << std::showpos << std::setprecision(6);
  WriteOut("coord.dat","index.dat");
};

System::~System(){
  for(int i=0;i<SLICES;i++){
    delete [] LINKS[i];
    delete [] NEXT[i];
    delete [] PREV[i];
    delete [] INDEX[i];
    delete [] NEXT_EMPTY[i];
  };
  delete [] LINKS;
  delete [] NEXT;
  delete [] PREV;
  delete [] INDEX;
  delete [] NEXT_EMPTY;
  delete [] SIDE;
  delete [] DX;
  delete [] NLNK;
  delete [] TRY;
  delete [] ACC;
  delete [] IRA;
  delete [] MASHA;
  delete [] NEXT_TAU;
  delete [] PREV_TAU;
  delete [] KINETIC;
  delete [] POTENTIAL;
  delete [] ENERGY;
  delete [] SUPER;
  delete [] GOFR;
  delete [] SKX_RE;
  delete [] SKX_IM;
  delete [] RHOXY;
  delete [] GOFXY;
  delete [] WEIGHT;
  delete [] COST;
};

void System::Read_Input(char* filename){
  std::string name;
  int idata;
  double rdata;
  std::ifstream file_in;
  file_in.open(filename);
  if(!file_in){ 
    std::cout << "The input file is missing! Please specify an existing file name!\n\n";
    exit(-1);
  };
  std::cout << "Reading input from " << filename << "... ";
  while(name!="END"){
    file_in >> name;
    if(name=="dim"){
      file_in >> idata;
      DIM = idata;
      SIDE = new double[DIM];
    };
    if(name=="n_slices"){
      file_in >> idata;
      Set_Slices(idata);
      TMAX = idata;
    };
    if(name=="max_part"){
      file_in >> idata;
      SetMASSPart(idata);
    };
    if(name=="restart"){
      file_in >> idata;
      RESTART = idata;
      if(idata == 0)
        SetSECTOR(0);
    };
    if(name=="n_part"){
      file_in >> idata;
      NUM = idata;
      SetMASSPart(NUM+1);
    };
    if(name=="T"){
      file_in >> rdata;
      Set_T(rdata);
      SetTAU(1.0/(T()*double(Slices())));
    };
    if(name=="tau"){
      file_in >> rdata;
      SetTAU(rdata);
      Set_T(1.0/(Tau()*double(Slices())));
    };
    if(name=="mbar"){
      file_in >> idata;
      SetMBAR(idata);
    };
    if(name=="mbridge"){
      file_in >> idata;
      SetMBRI(idata);
    };
    if(name=="c0"){
      file_in >> rdata;
      SetC0(rdata);
    };
    if(name=="m"){
      file_in >> rdata;
      SetMASS(rdata);
      LAMBDA=Hbar*Hbar*0.5/M();
    };
    if(name=="lambda"){
      file_in >> rdata;
      LAMBDA = rdata;
      MASS = Hbar*Hbar*0.5/LAMBDA;
    };
    if(name=="tmax"){
      file_in >> idata;
      TMAX = idata;
    };
    if(name=="density"){
      file_in >> rdata;
      SetVOL(double(NUM)/rdata);
      for(int i=0;i<DIM;i++)
        Set_L(i,pow(VOL,1./double(DIM)));
      DK = 2.0*Pi/L(0);
      RMAX=L(0)*0.5*sqrt(DIM);
    };
    if(name=="volume"){
      file_in >> rdata;
      SetVOL(rdata);
      for(int i=0;i<DIM;i++)
        Set_L(i,pow(rdata,1./double(DIM)));
      DK = 2.0*Pi/L(0);
      RMAX=L(0)*0.5*sqrt(DIM);
    };
    if(name=="side"){
      SIDE = new double[DIM];
      double v = 1.0;
      for(int i=0;i<DIM;i++){
        file_in >> rdata;
        Set_L(i,rdata);
        v = v * L(i);
      };
      SetVOL(v);
      RMAX=L(0)*0.5*sqrt(DIM);
      DK = 2.0*Pi/L(0);
    };
    if(name=="shift"){
      DX = new double[DIM];
      for(int i=0;i<DIM;i++){
        file_in >> rdata;
        DX[i] = rdata;
      };
    };
    if(name=="dr"){
      file_in >> rdata;
      DR = rdata;
      NG = int(RMAX*0.5/DR);
    };
    if(name=="nk"){
      file_in >> idata;
      NK = idata;
    };
    if(name=="seed"){
      file_in >> idata;
      RNG = Random(idata);
    };
    if(name=="lennardj"){
      Uint = LenJ;
      file_in >> rdata;
      COST[0] = rdata;
      file_in >> rdata;
      COST[1] = rdata;
    };
    if(name=="vdip"){
      Uint = VDip;
      file_in >> rdata;
      COST[0] = rdata;
    };
    if(name=="dipo"){
      Uint = Dip;
      file_in >> rdata;
      COST[0] = rdata;
      file_in >> rdata;
      COST[1] = cos(rdata/180.0*Pi);
    };
    if(name=="harm"){
      Uext = Harm;
    };
    if(name=="abc"){
      nint = abc;
      KinAndSup = &System::GetKin_SuperA;
      NG=int(SIDE[0]/DR);
      ECUT = 0.;
    };
    if(name=="pbc"){
      nint = pbc;
      KinAndSup = &System::GetKin_SuperW;
      ECUT = U0(SIDE[0]*0.5);
    };
    if(name=="nosup"){
      KinAndSup = &System::GetKin_0;
    };
    if(name=="blocks"){
      file_in >> idata;
      SetBLOCKS(idata);
      file_in >> idata;
      SetEQUI(idata);
      file_in >> idata;
      SetSPARSE(idata);
      file_in >> idata;
      SetSTEPS(idata);
    };
  };
  SIGMA = sqrt(Tau()*2.0*LAMBDA); 
  file_in.close();
  CI0 = CI0/(Vol()*double(Slices()*MBar()));
  std::cout << "Done!" << std::endl;
  return;
};

void System::Restart(int res){
  if(res==0)
    std::cout << "\nStarting from scratch... ";
  else if(res==1) 
    std::cout << "\nReading configurations from coord.dat, index.dat... ";
  else{
    std::cout << "\nSomething very wrong is happening in System::Restart(int res)!\n";
    exit(-1);
  }
  for(int i=0;i<Slices();i++){
      FIRST_EMPTY[i] = 0;
      for(int j=0;j<MPart();j++){
        NEXT_EMPTY[i][j] = (j+1)%MPart();
        INDEX[i][j] = -1;
        NEXT [i][j] = -1;
        PREV [i][j] = -1;
        for(int k=0;k<DIM;k++){
          LINKS[i][j].SetX0(k,0.0);
          LINKS[i][j].SetX1(k,0.0);
        };
      };
    };
  if(res==0){
    int m = Slices()+1;
    int n = NUM;
    double v[n][m][3];
    double vz[DIM];
    for(int i=0;i<DIM;i++){
      vz[i]=RNG.Double_Random(-0.5*L(i),0.5*L(i));
    };
    for(int i=0;i<NUM;i++){
      for(int j=0;j<DIM;j++){
        vz[j]=RNG.Double_Random(-0.5*L(j),0.5*L(j));
      };
      Bridge(vz,vz,v[i],m);
    };
    for(int i=0;i<Slices();i++){
      for(int j=0;j<NUM;j++){
        NEXT [i][j] = j;
        PREV [i][j] = j;
        AddLink(i,j,v[j][i],v[j][i+1]);
      };
    };
  };
  if(res==1){
    std::ifstream file_conf;
    std::ifstream file_ind;
    int t,ind,t0,ii;
    double x0;
    
    file_ind.open("index.dat");
    if(!file_ind){
      std::cout << "index.dat is missing! Restarting is impossible!\n\n";
      exit(-1);
    };
    file_ind >> t0;
    file_ind >> ii;
    SetSECTOR(ii);
    file_ind >> IRA[0]   >> IRA[1];
    file_ind >> MASHA[0] >> MASHA[1];
    file_ind >> WLINK;
    t = t0;
    for(int i=0;i<Slices();i++)
      for(int j=0;j<MAXPART;j++)
        file_ind >> NEXT_EMPTY[i][j];
    for(int i=0;i<Slices();i++){
      file_ind >> t >> NLNK[t] >> FIRST_EMPTY[t];
      for(int j=0;j<NLNK[t];j++){
        file_ind >> ind;
        file_ind >> NEXT[t][ind] >> PREV[t][ind];
        INDEX[t][j] = ind;
      };
      t = NEXT_TAU[t];
    };      
    file_ind.close();

    file_conf.open("coord.dat");
    if(!file_conf){
      std::cout << "coord.dat is missing! Restarting is impossible!\n\n";
      exit(-1);
    };
    for(int i=0;i<NLnk(t0);i++){
      t=t0;
      ind = INDEX[t][i];
      for(int j=0;j<Slices();j++){
        if(ind!=-1){
          file_conf >> t;
          for(int d=0;d<DIM;d++){
            file_conf >> x0;
            Links(t,ind).SetX0(d,x0);
          };
          for(int d=0;d<DIM;d++){
            file_conf >> x0;
            Links(t,ind).SetX1(d,x0);
          };
          ind = NEXT[t][ind];
          t   = NEXT_TAU[t];
        };
      };
    };
    file_conf.close();
  };
  PrintWorm();
  std::cout << "Done!\n";
  return;
};

void System::Set_Slices(int slices){
  SLICES = slices;
  return;
};

void System::SetMASSPart(int part){
  MAXPART = part;
  return;
};

void System::SetSECTOR(int sect){
  SECTOR = sect;
  return;
};

void System::Set_T(double t){
  TEMP = t;
  return;
};

void System::SetTAU(double tau){
  TAU = tau;
  return;
};

void System::SetMBAR(int mbar){
  MBAR = mbar;
  return;
};

void System::SetMBRI(int mbri){
  MBRI = mbri;
  return;
};

void System::SetC0(double c){
  CI0 = c;
  return;
};

void System::SetBLOCKS(int blocks){
  BLOCKS = blocks;
  return;
};

void System::SetEQUI(int equi){
  EQUI = equi;
  return;
};

void System::SetSPARSE(int sparse){
  SPARSE = sparse;
  return;
};

void System::SetSTEPS(int steps){
  STEPS = steps;
  return;
};

void System::SetMASS(double m){
  MASS = m;
  return;
};

void System::SetVOL(double v){
  VOL = v;
  return;
};

void System::Set_L(int i,double l){
  SIDE[i] = l;
  return;
};

void System::PrintSector() const {
  std::cout << Sector() << "\n";
  return;
};

void System::PrintWorm() const {
  if(Sector()==0) return;
  std::ofstream wfile;
  wfile.open("worm.dat");
  int tim0 = MASHA[0];
  int ind0 = MASHA[1];
  int tim = tim0;
  int ind = ind0;
  for(int i=0;i<WLink();i++){
    wfile << Links(tim0,ind0).X0(0) << "\t" << Links(tim0,ind0).X0(1) << "\t";
    if(tim0==MASHA[0]&&ind0==MASHA[1]) wfile << "M";
    wfile << "\n";
    if(tim0==IRA  [0]&&ind0==IRA  [1])
      wfile << Links(tim0,ind0).X1(0) << "\t" << Links(tim0,ind0).X1(1) << "\tI\n";
    tim = NEXT_TAU[tim0];
    ind = NEXT[tim0][ind0];
    tim0=tim;
    ind0=ind;
  };
  wfile.close();
  return;
};

void System::PrintPart(int p) const {
  int tim0 = 0;
  int ind0 = INDEX[0][p];
  int tim = tim0;
  int ind = ind0;
  for(int i=0;i<Slices();i++){
    std::cout << tim0 << "\t" << ind0 << "\t" << PREV[tim0][ind0] << "\t" << NEXT[tim0][ind0] << "\n";
    std::cout << Links(tim0,ind0).X0(0) << "\t" << Links(tim0,ind0).X0(1) << "\t";
    if(tim0==MASHA[0]&&ind0==MASHA[1]) std::cout << "M";
    std::cout << "\n";
    std::cout << Links(tim0,ind0).X1(0) << "\t" << Links(tim0,ind0).X1(1) << "\t";
    if(tim0==IRA  [0]&&ind0==IRA  [1]) std::cout << "I";
    std::cout << "\n";
    tim = NEXT_TAU[tim0];
    ind = NEXT[tim0][ind0];
    tim0=tim;
    ind0=ind;
  };
  return;
};

void System::WriteOut(char name_conf[],char name_ind[]) const{
  std::ofstream file_conf,file_ind;

  int t0 = MASHA[0];
  int t,ind;
  if(Sector()==0) t0 = 0;
  file_ind.open(name_ind);
  file_ind << t0 << "\t" << Sector() << "\n";
  file_ind << IRA[0]   << "\t" << IRA[1]   << "\n";
  file_ind << MASHA[0] << "\t" << MASHA[1] << "\n";
  file_ind << WLINK << "\n";
  t = t0;
  for(int i=0;i<Slices();i++)
    for(int j=0;j<MAXPART;j++)
      file_ind << NEXT_EMPTY[i][j] << "\t";
  file_ind << "\n";
  for(int i=0;i<Slices();i++){
    file_ind << t << "\t" << NLnk(t) << "\t" << FIRST_EMPTY[t] << "\t";
    for(int j=0;j<NLnk(t);j++){
      ind = INDEX[t][j];
      file_ind << ind << "\t" << NEXT[t][ind] << "\t" << PREV[t][ind] << "\t";
    };
    file_ind << "\n";
    t = NEXT_TAU[t];
  };
  file_ind.close();
  file_conf.open(name_conf);
  file_conf << std::scientific << std::setprecision(6) << std::showpos;
  for(int i=0;i<NLnk(t0);i++){
    t = t0;
    ind = INDEX[t][i];
    for(int j=0;j<Slices();j++){
      if(ind!=-1){
        file_conf << t << "\t";
        for(int d=0;d<DIM;d++){
          double pos = Links(t,ind).X0(d) - L(d)*nint(Links(t,ind).X0(d)/L(d));
          file_conf << pos << "\t";
        };
        for(int d=0;d<DIM;d++){
          double pos = Links(t,ind).X1(d) - L(d)*nint(Links(t,ind).X1(d)/L(d));
          file_conf << pos << "\t";
        }
        file_conf << "\n";
        ind = NEXT[t][ind];
        t   = NEXT_TAU[t];
      };
    };
    file_conf << "\n";
  };
  file_conf.close();
  return;
};

int System::Sector() const {
  return SECTOR;
};

int System::Slices() const {
  return SLICES;
};

int System::MPart() const {
  return MAXPART;
};

int System::MBar() const {
  return MBAR;
};

int System::MBri() const {
  return MBRI;
};

int System::Blocks() const {
  return BLOCKS;
};

int System::Equi() const {
  return EQUI;
};

int System::Sparse() const {
  return SPARSE;
};

int System::Steps() const {
  return STEPS;
};

int System::WLink() const {
  return WLINK;
};

int System::NP() const {
  if (Sector() == 0) return NLNK[0];
  else return NLNK[IRA[0]];
};

int System::NLnk(int i) const {
  return NLNK[i];
};

int System::Roulette(double vect[],int dim,double sig) const{
  double somma = 0.0;
  double rand  = RNG.Double_Random(0.0,1.0);
  for(int i=0;i<dim;i++){
    somma += vect[i];
    if(somma >= rand*sig)return i;
  };
  return -1;
};

double System::T() const {
  return TEMP;
};

double System::Tau() const {
  return TAU;
};

double System::C0() const {
  return CI0;
};

double System::M() const {
  return MASS;
};

double System::Vol() const {
  return VOL;
};

double System::L(int i) const{
  return SIDE[i];
};

double System::Sigma() const{
  return SIGMA;
};

double System::V(int time,int part,double x0[],double x1[]) const{
  double v = 0.0;
  double y0[DIM],y1[DIM];
  int np = NLNK[time];
  v = v + 0.5*Uext(x0,DIM)+0.5*Uext(x1,DIM);
  for(int n=0;n<np;n++){
    int ind = INDEX[time][n];
    if(ind!=part){
      for(int d=0;d<DIM;d++){
        y0[d] = Links(time,ind).X0(d);
        y1[d] = Links(time,ind).X1(d);
      };
      v = v + 0.5*Uint(x0,y0,DIM,SIDE,ECUT,COST)+0.5*Uint(x1,y1,DIM,SIDE,ECUT,COST);
    };
  };
  return v;
};

double System::Vpot(int time,int part,double x0[]) const{
  double v = 0.0;
  double y0[DIM];
  int np = NLNK[time];
  v = v + Uext(x0,DIM);
  for(int n=0;n<np;n++){
    int ind = INDEX[time][n];
    if(ind!=part){
      for(int d=0;d<DIM;d++){
        y0[d] = Links(time,ind).X0(d);
      };
      v = v + 0.5*Uint(x0,y0,DIM,SIDE,ECUT,COST);
    };
  };
  return v;
};

double System::V(int time,int part,Link l) const{
  double v =0.0;
  double a[DIM],b[DIM];
  for(int i=0;i<DIM;i++){
    a[i] = l.X0(i);
    b[i] = l.X1(i);
  };
  v = V(time,part,a,b);
  return v;
};

double System::Vpot(int time,int part,Link l) const{
  double v =0.0;
  double a[DIM];
  for(int i=0;i<DIM;i++){
    a[i] = l.X0(i);
  };
  v = Vpot(time,part,a);
  return v;
};

double System::U0(double r) const{
  double u = 0.0;
  double x0[DIM];
  double x1[DIM];
  for(int i=0;i<DIM;i++){
    x0[i] = 0.0;
    x1[i] = 0.0;
  };
  x1[0] = r;
  u = Uint(x0,x1,DIM,SIDE,0.0,COST);
  return u;
};

double System::Gauss(double d2,double dsigd)const{
  double esp = -double(DIM)*0.5;
  return pow(dsigd*Pi,esp)*exp(-d2/dsigd);
};

Link& System::Links(int i,int j) const {
  return LINKS[i][j];
};

void System::AddLink(int time,int index,double *x0,double *x1){
  NLNK[time]++;
  if(NLNK[time]>MPart()) std::cerr << "Troppe particelle!!\n";
  FIRST_EMPTY[time] = NEXT_EMPTY[time][index];
  NEXT_EMPTY[time][index] = -1;
  for(int i=0;i<DIM;i++){
    Links(time,index).SetX0(i,x0[i]);
    Links(time,index).SetX1(i,x1[i]);
  };
  INDEX[time][NLNK[time]-1] = index; 
  return;
};

void System::RemLink(int time,int ind){
  for(int i=0;i<NLnk(time);i++)
    if(INDEX[time][i]==ind){
      INDEX[time][i] = INDEX[time][NLNK[time]-1];
      INDEX[time][NLNK[time]-1] = -1;
      NEXT_EMPTY[time][ind] = FIRST_EMPTY[time];
      FIRST_EMPTY[time] = ind;
      NLNK[time]--;
      return;
    };
  std::cerr << "what is this i dont even \n";
  std::cerr << time << "\t" << ind << "\n";
  return;
};

void System::Bridge(double x0[],double xf[],double xnew[][3],int emme){
  emme--;
  for(int i=0;i<DIM;i++){
    xnew[0][i]    = x0[i];
    xnew[emme][i] = xf[i];
  };
  int l3=emme;
  for(int i=1;i<emme;i++){
    int l1 = i-1;
    int l2 = i;
    double a = double(l3-l2);
    double b = double(l3-l1);
    double s = Sigma()*sqrt(a/b);
    for(int j=0;j<DIM;j++){
      double d = xnew[l3][j]-xnew[l1][j];
      d = d - L(j)*nint(d/L(j));
      d = d/b;
      xnew[i][j] = xnew[l1][j]+RNG.Gauss_Random(d,s);
    };
  };
  return;
};

void System::Advance(){
  if(Sector()==0) return;
  int m,tim,tim0,ind,ind0;
  double Pacc,Rand;
  m = RNG.Int_Random(1,MBar());
  double xnew[m+1][2];
  xnew[0][0] = LINKS[IRA[0]  ][IRA[1]  ].X1(0);
  xnew[0][1] = LINKS[IRA[0]  ][IRA[1]  ].X1(1);
  for(int i=0;i<m;i++){
    for(int j=0;j<DIM;j++){
      xnew[i+1][j] = RNG.Gauss_Random(xnew[i][j],Sigma());
    };
  };
  Pacc = 0.0;
  tim0 =-1;
  ind0 =-1;
  tim = MASHA[0];
  ind = MASHA[1];
  for(int i=0;i<m;i++){
    tim0=tim;
    ind0=ind;
    Pacc = Pacc + (V(tim0,ind0,Links(tim0,ind0))-V(tim0,ind0,xnew[i],xnew[i+1]))*Tau();
    ind = NEXT[tim0][ind0];
    tim = NEXT_TAU[tim0];
  };
  Pacc = exp(Pacc);
  Rand = RNG.Double_Random(0.0,1.0);
  TRY[0]++;
  if(Rand < Pacc){
    ACC[0]++;
    tim0 =-1;
    ind0 =-1;
    tim = MASHA[0];
    ind = MASHA[1];
    NEXT[IRA[0]][IRA[1]] = ind;
    PREV[tim][ind] = IRA[1];
    for(int i=0;i<m;i++){
      tim0=tim;
      ind0=ind;
      RemLink(tim0,ind0);
      AddLink(tim0,ind0,xnew[i],xnew[i+1]);
      tim = NEXT_TAU[tim0];
      ind = NEXT[tim0][ind0];
      NEXT[tim0][ind0] = ind;
      PREV[tim ][ind ] = ind0;
    };
    IRA[0]=tim0;
    IRA[1]=ind0;
    MASHA[0] = tim;
    MASHA[1] = ind;
    NEXT[IRA[0]  ][IRA[1]  ] = -1;
    PREV[MASHA[0]][MASHA[1]] = -1;
  };
  return;
};

void System::Recede(){
  if(Sector()==0) return;
  int m;
  m = RNG.Int_Random(1,MBar());
  if(m>=WLink()) return;
  double Pacc,Rand;
  double xnew[m+1][2];
  int tim0=-1,tim;
  int ind0=-1,ind;
  tim = IRA[0];
  ind = IRA[1];
  xnew[0][0] = LINKS[MASHA[0]][MASHA[1]].X0(0);
  xnew[0][1] = LINKS[MASHA[0]][MASHA[1]].X0(1);
  for(int i=1;i<m+1;i++){
    for(int j=0;j<DIM;j++){
      xnew[i][j] = RNG.Gauss_Random(xnew[i-1][j],Sigma());
    };
  };
  Pacc = 0.0;
  for(int i=0;i<m;i++){
    tim0 = tim;
    ind0 = ind;
    Pacc = Pacc + (V(tim0,ind0,Links(tim0,ind0))-V(tim0,ind0,xnew[i+1],xnew[i]))*Tau();
    ind = PREV[tim0][ind0];
    tim = PREV_TAU[tim0];
  };
  Pacc = exp(Pacc);
  Rand = RNG.Double_Random(0.0,1.0);
  TRY[1]++;
  if(Rand < Pacc){
    ACC[1]++;
    tim = IRA[0];
    ind = IRA[1];
    NEXT[IRA[0]  ][IRA[1]  ] = MASHA[1]; 
    PREV[MASHA[0]][MASHA[1]] = IRA[1];
    for(int i=0;i<m;i++){
      tim0=tim;
      ind0=ind;
      RemLink(tim0,ind0);
      AddLink(tim0,ind0,xnew[i+1],xnew[i]);
      ind = PREV[tim0][ind0];
      tim = PREV_TAU[tim0];
    };
    IRA[0]   = tim;
    IRA[1]   = ind;
    MASHA[0] = tim0;
    MASHA[1] = ind0;
    NEXT[IRA[0]  ][IRA[1]  ]=-1;
    PREV[MASHA[0]][MASHA[1]]=-1;
  };
  return;
};

void System::Close(){
  if(Sector()==0) return;
  int m = RNG.Int_Random(1,MBar()+1);
  int tst,ist;
  tst = MASHA[0];
  ist = MASHA[1];
  for(int i=0;i<m;i++){
    ist = NEXT[tst][ist];
    tst = NEXT_TAU[tst];
  };
  double xnew[m+1][3];
  double XIra[DIM],XSt[DIM];
  for(int i=0;i<DIM;i++){
    XIra[i] = Links(IRA[0],IRA[1]).X1(i);
    XSt [i] = Links(tst   ,ist   ).X0(i);
  };
  double d2 = 0.0;
  double r;
  for(int i=0;i<DIM;i++){
    r = XIra[i]-XSt[i];
    r = r-L(i)*nint(r/L(i));
    d2 = d2+r*r;
  };
  double dsigd=2.0*Sigma()*Sigma();
  if(d2/dsigd>4.0) return;
  double rho = Gauss(d2,dsigd);
  int M = m+1;
  Bridge(XIra,XSt,xnew,M);
  double Pacc=0.0;
  int nb=NLnk(IRA[0])*Slices();
  double den = C0()*double(MBar())*double(nb);
  int tim0,tim,ind0,ind;
  tim = MASHA[0];
  ind = MASHA[1];
  for(int i=0;i<m;i++){
    tim0 = tim;
    ind0 = ind;
    Pacc = Pacc +(V(tim0,ind0,Links(tim0,ind0))- V(tim0,ind0,xnew[i],xnew[i+1]))*Tau();
    ind = NEXT[tim0][ind0];
    tim = NEXT_TAU[tim0];
  };
  Pacc = rho*exp(Pacc)/den;
  double Rand;
  Rand = RNG.Double_Random(0.0,1.0);
  TRY[3]++;
  if(Rand < Pacc){
    ACC[3]++;
    SetSECTOR(0);
    WLINK  = 0;
    tim = MASHA[0];
    ind = MASHA[1];
    NEXT[IRA[0]][IRA[1]] = ind;
    PREV[tim][ind]= IRA[1];
    for(int i=0;i<m;i++){
      tim0 = tim;
      ind0 = ind;
      RemLink(tim0,ind0);
      AddLink(tim0,ind0,xnew[i],xnew[i+1]);
      tim = NEXT_TAU[tim0];
      ind = NEXT[tim0][ind0];
      PREV[tim] [ind]  = ind0;
      NEXT[tim0][ind0] = ind;
    };
    IRA[0]=-1;
    IRA[1]=-1;
    MASHA[0]=-1;
    MASHA[1]=-1;
  };
  return;
};

void System::Open(){
  if(Sector()==1) return;
  int m,Start,Pol,Index;
  m     = RNG.Int_Random(1,MBar()+1);
  Start = RNG.Int_Random(0,Slices());
  Pol   = RNG.Int_Random(0,NLnk(Start));
  Index = INDEX[Start][Pol];
  int newirat,newirai,newmashat,newmashai;
  int ind0,ind,tim0,tim;
  newirat = PREV_TAU[Start];
  newirai = PREV[Start][Index];
  newmashat=Start;
  newmashai=Index;
  double XIra[DIM],XMasha[DIM];
  for(int i=0;i<DIM;i++){
    XIra  [i] = Links(newirat,newirai).X1(i);
  };
  double xnew[m+1][DIM];
  ind = newmashai;
  tim = newmashat;
  for(int i=0;i<m;i++){
    ind = NEXT[tim][ind];
    tim = NEXT_TAU[tim];
  };
  for(int j=0;j<DIM;j++)
    xnew[m][j] = RNG.Gauss_Random(Links(tim,ind).X0(j),Sigma());
  for(int i=1;i<m+1;i++)
    for(int j=0;j<DIM;j++)
      xnew[m-i][j] = RNG.Gauss_Random(xnew[m-i+1][j],Sigma());
  for(int j=0;j<DIM;j++)
    XMasha[j] = xnew[0][j];
  double d2 = 0.0;
  double dsigd = 2.0*Sigma()*Sigma();
  double r;
  for(int i=0;i<DIM;i++){
    r = XIra[i]-XMasha[i];
    r = r-L(i)*nint(r/L(i));
    d2 = d2+r*r;
  };
  if(d2/dsigd>4.0) return;
  double rho = Gauss(d2,dsigd);
  ind = Index;
  tim = Start;
  int nb=NLnk(tim)*Slices();
  double num=C0()*double(MBar())*double(nb);
  double Pacc = 0.0;
  for(int i=0;i<m;i++){
    tim0=tim;
    ind0=ind;
    Pacc = Pacc + (V(tim0,ind0,Links(tim0,ind0))-V(tim0,ind0,xnew[i],xnew[i+1]))*Tau();
    ind=NEXT[tim0][ind0];
    tim=NEXT_TAU[tim0];
  };
  Pacc = exp(Pacc)*num/rho;
  double Rand;
  Rand = RNG.Double_Random(0.0,1.0);
  TRY[2]++;
  if(Rand<Pacc){
    ACC[2]++;
    ind = Index;
    tim = Start;
    for(int i=0;i<m;i++){
      tim0=tim;
      ind0=ind;
      RemLink(tim0,ind0);
      AddLink(tim0,ind0,xnew[i],xnew[i+1]);
      ind=NEXT[tim0][ind0];
      tim=NEXT_TAU[tim0];
    };
    MASHA[0]=newmashat;
    MASHA[1]=newmashai;
    IRA[0]  =newirat;
    IRA[1]  =newirai;
    NEXT[IRA[0]  ][IRA[1]  ] = -1;
    PREV[MASHA[0]][MASHA[1]] = -1;
    tim=MASHA[0];
    ind=MASHA[1];
    WLINK=1;
    SetSECTOR(1);
    while(!(tim==IRA[0] && ind==IRA[1])){
      tim0=tim;
      ind0=ind;
      WLINK++;
      tim = NEXT_TAU[tim0];
      ind = NEXT[tim0][ind0];
    };
    if(WLINK%Slices()!=0) std::cerr << "Oh noes!!\n";
  };
  return;
};

void System::Swap(){
  if(Sector()==0) return;
  int m;
  m = RNG.Int_Random(1,MBar()+1);
  int tfinal = IRA[0];
  for(int i=0;i<m;i++)
    tfinal = NEXT_TAU[tfinal];
  if(NLNK[tfinal]==0) {
    return;
  };
  int alpha,zeta,tz,target;
  double pal[NLNK[tfinal]];
  int  indal[NLNK[tfinal]];
  int narr=0;
  double dsigd = 4.0*LAMBDA*Tau()*m;
  double sigira = 0.0,sigzeta =0.0;
  for(int i=0;i<NLNK[tfinal];i++){
    int ial = INDEX[tfinal][i];
    double d2 = 0.0;
    double r;
    for(int i=0;i<DIM;i++){
      r = Links(IRA[0],IRA[1]).X1(i)-Links(tfinal,ial).X1(i);
      r = r-L(i)*nint(r/L(i));
      d2 = d2+r*r;
    };
    if(d2<4.0/dsigd){
      pal[narr] = Gauss(d2,dsigd);
      indal[narr] = ial;
      narr++;
    };
  };
  if(narr==0) return;
  for(int i=0;i<narr;i++){
    sigira = sigira + pal[i];
  };
  int asd = Roulette(pal,narr,sigira);
  target = indal[asd];
  alpha = target;
  zeta   = alpha;
  tz = tfinal;
  for(int i=0;i<m;i++){
    zeta = PREV[tz][zeta];
    tz = PREV_TAU[tz];
    if(zeta==-1){
      return;
    };
  };
  if(tz!=IRA[0]) std::cerr << "WTF\n";
  if(NEXT[tz][zeta]==-1) return;
  if(tz==MASHA[0] && zeta==MASHA[1]){
    return;
  };
  narr=0;
  for(int i=0;i<NLNK[tfinal];i++){
    int ial = INDEX[tfinal][i];
    double d2 = 0.0;
    double r;
    for(int i=0;i<DIM;i++){
      r = Links(tz,zeta).X1(i)-Links(tfinal,ial).X1(i);
      r = r-L(i)*nint(r/L(i));
      d2 = d2+r*r;
    };
    if(d2<4.0){
      pal[narr] = Gauss(d2,dsigd);
      narr++;
    };
  };
  if(narr==0) return;
  for(int i=0;i<narr;i++){
    sigzeta = sigzeta + pal[i];
  };
  double xnew[m+1][3];
  double XIra[DIM],XAlpha[DIM];
  for(int i=0;i<DIM;i++){
    XIra  [i] = Links(IRA[0],IRA[1]).X1(i);
    XAlpha[i] = Links(tfinal ,alpha  ).X1(i);
  };
  int M = m+1;
  Bridge(XIra,XAlpha,xnew,M);
  int ind,ind0,tim,tim0;
  tim = NEXT_TAU[IRA[0]];
  ind = NEXT[IRA[0]][zeta];
  double Pacc = 0.0;
  for(int i=0;i<m;i++){
    tim0=tim;
    ind0=ind;
    Pacc = Pacc + (V(tim0,ind0,Links(tim0,ind0))-V(tim0,ind0,xnew[i],xnew[i+1]))*Tau();
    tim = NEXT_TAU[tim0];
    ind = NEXT[tim0][ind0];
  };
  Pacc = exp(Pacc)*sigira/sigzeta;
  double Rand;
  Rand = RNG.Double_Random(0.0,1.0);
  TRY[4]++;
  if(Rand < Pacc){
    tim = NEXT_TAU[IRA[0]];
    ind = NEXT[IRA[0]][zeta];
    if(ind==-1) return;
    ACC[4]++;
    NEXT[IRA[0]][IRA[1]]=ind;
    PREV[tim][ind]=IRA[1];
    for(int i=0;i<m;i++){
      tim0=tim;
      ind0=ind;
      RemLink(tim0,ind0);
      AddLink(tim0,ind0,xnew[i],xnew[i+1]);
      tim = NEXT_TAU[tim0];
      ind = NEXT[tim0][ind0];
      NEXT[tim0][ind0] = ind;
      PREV[tim] [ind]  = ind0;
    };
    NEXT[tz][zeta] = -1;
    IRA[0]=tz;
    IRA[1]=zeta;
    tim=MASHA[0];
    ind=MASHA[1];
    WLINK=1;
    while(!(tim==IRA[0] && ind==IRA[1])){
      tim0=tim;
      ind0=ind;
      WLINK++;
      tim = NEXT_TAU[tim0];
      ind = NEXT[tim0][ind0];
    };
    if(WLINK%Slices()!=0) std::cerr << "Oh snap!!\n";
  };
  return;
};

void System::Wiggle(){
  if(Sector()==1) return;
  int part  = RNG.Int_Random(0,NLnk(0));
  int tstart = RNG.Int_Random(0,Slices());
  int istart = INDEX[tstart][part];
  int m = RNG.Int_Random(1,MBri());
  int tend = tstart;
  int iend = istart;
  for(int i=0;i<m;i++){
    iend = NEXT[tend][iend];
    tend = NEXT_TAU[tend];
  };
  double xs[DIM],xe[DIM],xnew[m+2][3];
  for(int i=0;i<DIM;i++){
    xs[i] = Links(tstart,istart).X0(i);
    xe[i] = Links(tend  ,iend  ).X1(i);
  };
  int m1 = m+2;
  Bridge(xs,xe,xnew,m1);
  TRY[5]++;
  double Pacc = 0.0;
  int tim,tim0,ind,ind0;
  tim = tstart;
  ind = istart;
  for(int i=0;i<m+1;i++){
    tim0=tim;
    ind0=ind;
    Pacc = Pacc + (V(tim0,ind0,Links(tim0,ind0))-V(tim0,ind0,xnew[i],xnew[i+1]))*Tau();
    tim = NEXT_TAU[tim0];
    ind = NEXT[tim0][ind0];
  };
  Pacc = exp(Pacc);
  double Rand;
  Rand = RNG.Double_Random(0.0,1.0);
  if(Rand < Pacc){
    ACC[5]++;
    tim = tstart;
    ind = istart;
    for(int i=0;i<m+1;i++){
      tim0=tim;
      ind0=ind;
      RemLink(tim0,ind0);
      AddLink(tim0,ind0,xnew[i],xnew[i+1]);
      tim = NEXT_TAU[tim0];
      ind = NEXT[tim0][ind0];
    };
  };
  return;
};

void System::Shift(){
  if(Sector()==1) return;
  if(NLnk(0) ==0) return;
  TRY[6]++;
  int part  = RNG.Int_Random(0,NLnk(0));
  int len = 1;
  int tim0,tim,ind0,ind;
  tim0 = 0;
  ind0 = INDEX[tim0][part];
  tim = NEXT_TAU[tim0];
  ind = NEXT[tim0][ind0];
  while(!(tim==tim0 && ind==ind0)){
    len++;
    ind = NEXT[tim][ind];
    tim = NEXT_TAU[tim];
  };
  double dx[DIM];
  for(int i=0;i<DIM;i++)
    dx[i]=RNG.Double_Random(-DX[i],DX[i]);
  double xnew[len][DIM];
  tim=0;
  ind=INDEX[tim][part];
  for(int i=0;i<len;i++){
    tim0=tim;
    ind0=ind;
    for(int d=0;d<DIM;d++){
      xnew[i][d]=dx[d]+Links(tim0,ind0).X0(d);
    };
    ind = NEXT[tim0][ind0];
    tim = NEXT_TAU[tim0];
  };
  double Pacc = 0.0;
  tim=0;
  ind=INDEX[tim][part];
  for(int i=0;i<len;i++){
    tim0=tim;
    ind0=ind;
    int ii = (i+1)%len;
    Pacc = Pacc + (V(tim0,ind0,Links(tim0,ind0))-V(tim0,ind0,xnew[i],xnew[ii]))*Tau();
    tim = NEXT_TAU[tim0];
    ind = NEXT[tim0][ind0];
  };
  Pacc = exp(Pacc);
  double Rand;
  Rand = RNG.Double_Random(0.0,1.0);
  if(Rand<Pacc){
    ACC[6]++;
    tim=0;
    ind=INDEX[tim][part];
    for(int i=0;i<len;i++){
      tim0=tim;
      ind0=ind;
      int ii = (i+1)%len;
      RemLink(tim0,ind0);
      AddLink(tim0,ind0,xnew[i],xnew[ii]);
      tim = NEXT_TAU[tim0];
      ind = NEXT[tim0][ind0];
    };
  };
  return;
};

void System::Move(){
  int dice;
  dice = RNG.Int_Random(0,NMOV);
  if(dice==0) Advance();
  if(dice==1) Recede();
  if(dice==2) Open();
  if(dice==3) Close();
  if(dice==4) Swap();
  if(dice==5) Wiggle();
  if(dice==6) Shift();
  return;
};

void System::Measure(int block){
  if(Sector()==1) return;
  WEIGHT[block]++;
  if(NLnk(0)!=0){
    KINETIC  [0][block] += (*this.*KinAndSup)(block)/double(NLnk(0));
    POTENTIAL[0][block] += GetPot()/double(NLnk(0));
  };
  GetGofR(block);
  GetSofK(block);
  GetRho(block);
  GetGxy(block);
  return;
};

double System::GetKin_0(int b){
  if(NLnk(0)==0) return 0.0;
  double k = NLnk(0)*double(DIM)*0.5/Tau();
  double d2 = 0.0;
  double r,asup,mi,wn[DIM];
  int ind,tim;
  asup = 0.0;
  mi   = 0.0;
  for(int i=0;i<DIM;i++)
    wn[i] = 0.0;
  for(int i=0;i<NLnk(0);i++){
    ind = INDEX[0][i];
    tim = 0;
    for(int l=0;l<Slices();l++){
      for(int d=0;d<DIM;d++){
        r = Links(tim,ind).X1(d)-Links(tim,ind).X0(d);
        r = r - L(d)*nint(r/L(d));
        d2 = d2 - r*r;
      };
      ind = NEXT[tim][ind];
      tim = NEXT_TAU[tim];
    };
  };
  d2 = d2/(double(Slices())*4.0*LAMBDA*Tau()*Tau());
  k = k + d2;
  return k;
};

double System::GetKin_SuperW(int b){
  if(NLnk(0)==0) return 0.0;
  double k = NLnk(0)*double(DIM)*0.5/Tau();
  double d2 = 0.0;
  double r,asup,mi,wn[DIM];
  int ind,tim;
  asup = 0.0;
  mi   = 0.0;
  for(int i=0;i<DIM;i++)
    wn[i] = 0.0;
  for(int i=0;i<NLnk(0);i++){
    ind = INDEX[0][i];
    tim = 0;
    for(int l=0;l<Slices();l++){
      for(int d=0;d<DIM;d++){
        r = Links(tim,ind).X1(d)-Links(tim,ind).X0(d);
        r = r - L(d)*nint(r/L(d));
        wn[d] = wn[d]+r;
        d2 = d2 - r*r;
      };
      ind = NEXT[tim][ind];
      tim = NEXT_TAU[tim];
    };
  };
  d2 = d2/(double(Slices())*4.0*LAMBDA*Tau()*Tau());
  k = k + d2;
  for(int i=0;i<DIM;i++)
    SUPER[0][b] += wn[i]*wn[i]/(2.0*LAMBDA*double(Slices())*Tau()*double(NLnk(0)));
  return k;
};

double System::GetKin_SuperA(int b){
  if(NLnk(0)==0) return 0.0;
  double k = NLnk(0)*double(DIM)*0.5/Tau();
  double d2 = 0.0;
  double r,asup,mi,wn[DIM];
  int ind,tim;
  asup = 0.0;
  mi   = 0.0;
  for(int i=0;i<DIM;i++)
    wn[i] = 0.0;
  for(int i=0;i<NLnk(0);i++){
    ind = INDEX[0][i];
    tim = 0;
    for(int l=0;l<Slices();l++){
      asup = asup - Links(tim,ind).X1(0)*Links(tim,ind).X0(1) + Links(tim,ind).X1(1)*Links(tim,ind).X0(0);
      mi   = mi   + Links(tim,ind).X1(0)*Links(tim,ind).X0(0) + Links(tim,ind).X1(1)*Links(tim,ind).X0(1);
      for(int d=0;d<DIM;d++){
        r = Links(tim,ind).X1(d)-Links(tim,ind).X0(d);
        r = r - L(d)*nint(r/L(d));
        d2 = d2 - r*r;
      };
      ind = NEXT[tim][ind];
      tim = NEXT_TAU[tim];
    };
  };
  d2 = d2/(double(Slices())*4.0*LAMBDA*Tau()*Tau());
  k = k + d2;
  asup = asup * 0.5;
  mi   = mi/double(Slices());
  SUPER[0][b] += 4.0*asup*asup/mi/(Tau()*double(Slices()));
  return k;
};

double System::GetPot() const{
  if(NLnk(0)==0) return 0.0;
  double u = 0.0;
  double kount = 0.0;
  int ind;
  for(int t=0;t<Slices();t=t+10){
    kount = kount + 1.0;
    for(int i=0;i<NLnk(0);i++){
      ind = INDEX[t][i];
      u = u + Vpot(t,ind,Links(t,ind));
    };
  };
  u = u /kount;
  return u;
};

void System::GetGofR(int block){
  int ind1,ind2,bin;
  double d,r;
  for(int t=0;t<Slices();t++){
    for(int p=0;p<NLnk(t);p++){
      ind1 = INDEX[t][p];
      for(int q=p+1;q<NLnk(t);q++){
        ind2 = INDEX[t][q];
        d = 0.0;
        for(int b=0;b<DIM;b++){
          r = Links(t,ind1).X0(b)-Links(t,ind2).X0(b);
          r = r-L(b)*nint(r/L(b));
          d+= r*r;
        };
        d = sqrt(d);
        bin = std::min(int(d/DR),NG);
        GOFR[0][block][bin]++;
      };
    };
  };
  return;
};

void System::GetSofK(int block){
  int ind,t1;
  double x1,cos1[Slices()],sin1[Slices()],k;
  for(int a=0;a<NK;a++){
    k = a*DK;
    for(int t=0;t<Slices();t++){
      cos1[t] = 0.0;
      sin1[t] = 0.0;
    };
    for(int t=0;t<Slices();t++){
      for(int p=0;p<NLnk(t);p++){ 
        ind = INDEX[t][p];
        x1 = Links(t,ind).X0(0);
        cos1[t] += cos(k*x1);
        sin1[t] += sin(k*x1);
      };
      
    };
    for(int t=0;t<Slices();t++)
      for(int s=0;s<TMAX;s++){
        t1 = (t+s)%Slices();
        SKX_RE[0][block][a][s] = SKX_RE[0][block][a][s] + 
               (cos1[t]*cos1[t1]+sin1[t]*sin1[t1])/double(Slices()*NLnk(0)); ;
        SKX_IM[0][block][a][s] = SKX_IM[0][block][a][s] - 
               (cos1[t]*sin1[t1]-sin1[t]*cos1[t1])/double(Slices()*NLnk(0)); ;
      };
  };
  return;
};

void System::GetRho(int block){
  int ind,ix,iy;
  double x,y;
  for(int t=0;t<Slices();t++){
    for(int p=0;p<NLnk(t);p++){
      ind = INDEX[t][p];
      x = Links(t,ind).X0(0);
      y = Links(t,ind).X0(1);
      ix = int((0.5*L(0)+x)/DX1);
      iy = int((0.5*L(1)+y)/DX2);
      if(ix<NX && iy<NY)
        if(ix>-1 && iy>-1)
          RHOXY[0][block][ix][iy]=RHOXY[0][block][ix][iy]+1.0/double(Slices());
    };
  };
  return;
};

void System::GetGxy(int block){
  int ind1,ind2,ix,iy;
  double x,y;
  for(int t=0;t<Slices();t++){
    for(int p=0;p<NLnk(t);p++){
      ind1 = INDEX[t][p];
      for(int q=p+1;q<NLnk(t);q++){
        ind2 = INDEX[t][q];
        x = Links(t,ind1).X0(0)-Links(t,ind2).X0(0);
        x = x - L(0)*nint(x/L(0));
        y = Links(t,ind1).X0(1)-Links(t,ind2).X0(1);
        y = y - L(1)*nint(x/L(1));
        ix = int(floor((0.5*L(0)+x)/DX1)+0.5);
        iy = int(floor((0.5*L(1)+y)/DX2)+0.5);
        if(ix<NX && iy<NY)
          if(ix>-1 && iy>-1)
            GOFXY[0][block][ix][iy]=GOFXY[0][block][ix][iy]+1.0/double(Slices());
      };
    };
  };
  return;
};


void System::PrintG(int block) const{
  std::ofstream file_g;
  file_g.open("gofr.dat");
  file_g << std::scientific << std::showpos << std::setprecision(6);
  for(int i=0;i<NG;i++)
    file_g << i*DR << "\t" << GOFR[0][block][i] << "\t" << GOFR[1][block][i] << "\t" << GOFR[2][block][i] <<"\n";
  file_g.close();
  return;
};

void System::PrintS(int block) const{
  std::ofstream file_sr,file_si;
  file_sr.open("skx_re.dat");
  file_si.open("skx_im.dat");
  file_sr << std::scientific << std::showpos << std::setprecision(6);
  file_si << std::scientific << std::showpos << std::setprecision(6);
  for(int i=0;i<NK;i++){
    for(int j=0;j<TMAX;j++){
      file_sr << i*DK << "\t" << j << "\t" << SKX_RE[0][block][i][j] << "\t"
              << SKX_RE[1][block][i][j] << "\t" << SKX_RE[2][block][i][j] <<"\n";
      file_si << i*DK << "\t" << j << "\t" << SKX_IM[0][block][i][j] << "\t"
              << SKX_IM[1][block][i][j] << "\t" << SKX_IM[2][block][i][j] <<"\n";
    };
    if(TMAX>1){
      file_sr << "\n";
      file_si << "\n";
    };
  };
  file_sr.close();
  file_si.close();
  return;
};

void System::PrintRho(int block) const{
  std::ofstream file_rho;
  file_rho << std::scientific << std::showpos << std::setprecision(6);
  file_rho.open("rho_xy.dat");
  for(int i=1;i<NX;i++){
    for(int j=1;j<NY;j++)
      file_rho << -L(0)*0.5+i*DX1 << "\t" << -L(1)*0.5+j*DX2 << "\t"
               << RHOXY[1][block][i][j] << "\t" << RHOXY[2][block][i][j] << "\n";
    file_rho << "\n";
  };
  file_rho.close();
  return;
};

void System::PrintGxy(int block) const{
  std::ofstream file_rho;
  file_rho << std::scientific << std::showpos << std::setprecision(6);
  file_rho.open("g_xy.dat");
  for(int i=1;i<NX;i++){
    for(int j=1;j<NY;j++)
      file_rho << -L(0)*0.5+i*DX1 << "\t" << -L(1)*0.5+j*DX2 << "\t"
               << GOFXY[1][block][i][j] << "\t" << GOFXY[2][block][i][j] << "\n";
    file_rho << "\n";
  };
  file_rho.close();
  return;
};

void System::Average(int block){
  double kin2=0.0,pot2=0.0,tot2=0.0,num2=0.0,sup2=0.0,swn2=0.0;
  double gor2[NG],skr2[NK][TMAX],ski2[NK][TMAX],rho2[NX][NY],gxy2[NX][NY];
  NormGofR(block);
  for(int g=0;g<NG;g++)
    gor2[g] = 0.0;
  for(int k=0;k<NK;k++)
    for(int t=0;t<TMAX;t++){
      skr2[k][t] = 0.0;
      ski2[k][t] = 0.0;
    };
  for(int i=0;i<NX;i++)
    for(int j=0;j<NY;j++){
      rho2[i][j] = 0.0;
      gxy2[i][j] = 0.0;
  };
  ENERGY   [0][block] = ENERGY[0][block]+KINETIC[0][block]+POTENTIAL[0][block];
  if(WEIGHT[block]!=0.0){
    BEFF++;
    KINETIC  [0][block] = KINETIC  [0][block]/WEIGHT[block];
    POTENTIAL[0][block] = POTENTIAL[0][block]/WEIGHT[block];
    ENERGY   [0][block] = ENERGY   [0][block]/WEIGHT[block];
    SUPER    [0][block] = SUPER    [0][block]/WEIGHT[block];
    for(int i=0;i<NG;i++)
      GOFR[0][block][i] = GOFR[0][block][i]/WEIGHT[block];
    for(int i=0;i<NK;i++){
      for(int t=0;t<TMAX;t++){
      SKX_RE[0][block][i][t] = SKX_RE[0][block][i][t]/WEIGHT[block];
      SKX_IM[0][block][i][t] = SKX_IM[0][block][i][t]/WEIGHT[block];
      };
    };
    for(int i=0;i<NX;i++){
      for(int j=0;j<NY;j++){
        RHOXY[0][block][i][j] = RHOXY[0][block][i][j]/WEIGHT[block];
        GOFXY[0][block][i][j] = GOFXY[0][block][i][j]/WEIGHT[block];
      };
    };
  };
  double norm  = 0.0;
  for(int h=0;h<=block;h++){
    norm  += WEIGHT[h];

    KINETIC  [1][block] += KINETIC  [0][h]*WEIGHT[h];
    POTENTIAL[1][block] += POTENTIAL[0][h]*WEIGHT[h];
    ENERGY   [1][block] += ENERGY   [0][h]*WEIGHT[h];
    SUPER    [1][block] += SUPER [0][h]*WEIGHT[h];
    
    kin2 += KINETIC  [0][h]*KINETIC  [0][h]*WEIGHT[h];
    pot2 += POTENTIAL[0][h]*POTENTIAL[0][h]*WEIGHT[h];
    tot2 += ENERGY   [0][h]*ENERGY   [0][h]*WEIGHT[h];
    swn2 += SUPER    [0][h]*SUPER    [0][h]*WEIGHT[h];
    
    for(int g=0;g<NG;g++){
      GOFR[1][block][g] += GOFR[0][h][g]*WEIGHT[h];
      gor2[g] += GOFR[0][h][g]*GOFR[0][h][g]*WEIGHT[h];
    };
    for(int k=0;k<NK;k++){
      for(int t=0;t<TMAX;t++){
        SKX_RE[1][block][k][t] += SKX_RE[0][h][k][t]*WEIGHT[h];
        skr2[k][t] += SKX_RE[0][h][k][t]*SKX_RE[0][h][k][t]*WEIGHT[h];
        SKX_IM[1][block][k][t] += SKX_IM[0][h][k][t]*WEIGHT[h];
        ski2[k][t] += SKX_IM[0][h][k][t]*SKX_IM[0][h][k][t]*WEIGHT[h];
      };
    };
    for(int i=0;i<NX;i++){
      for(int j=0;j<NY;j++){
        RHOXY[1][block][i][j] += RHOXY[0][h][i][j]*WEIGHT[h];
        rho2[i][j] += RHOXY[0][h][i][j]*RHOXY[0][h][i][j]*WEIGHT[h];
        GOFXY[1][block][i][j] += GOFXY[0][h][i][j]*WEIGHT[h];
        gxy2[i][j] += GOFXY[0][h][i][j]*GOFXY[0][h][i][j]*WEIGHT[h];
      };
    };
  };

  KINETIC  [1][block] = KINETIC  [1][block]/norm;
  POTENTIAL[1][block] = POTENTIAL[1][block]/norm;
  ENERGY   [1][block] = ENERGY   [1][block]/norm;
  SUPER    [1][block] = SUPER    [1][block]/norm;
  kin2 = kin2/norm;
  pot2 = pot2/norm;
  tot2 = tot2/norm;
  num2 = num2/norm;
  sup2 = sup2/norm;
  swn2 = swn2/norm;
  for(int g=0;g<NG;g++){
    gor2[g] = gor2[g]/norm;
    GOFR[1][block][g] = GOFR[1][block][g]/norm;
  };
  for(int g=0;g<NK;g++){
    for(int t=0;t<TMAX;t++){
    skr2[g][t] = skr2[g][t]/norm;
    SKX_RE[1][block][g][t] = SKX_RE[1][block][g][t]/norm;
    ski2[g][t] = ski2[g][t]/norm;
    SKX_IM[1][block][g][t] = SKX_IM[1][block][g][t]/norm;
    };
  };
  for(int i=0;i<NX;i++){
    for(int j=0;j<NY;j++){
      rho2[i][j] = rho2[i][j]/norm;
      RHOXY[1][block][i][j] = RHOXY[1][block][i][j]/norm;
      gxy2[i][j] = gxy2[i][j]/norm;
      GOFXY[1][block][i][j] = GOFXY[1][block][i][j]/norm;
    };
  };
  if(block>0){
      KINETIC  [2][block] = sqrt(kin2-KINETIC  [1][block]*KINETIC  [1][block])/sqrt(BEFF);
      POTENTIAL[2][block] = sqrt(pot2-POTENTIAL[1][block]*POTENTIAL[1][block])/sqrt(BEFF);
      ENERGY   [2][block] = sqrt(tot2-ENERGY   [1][block]*ENERGY   [1][block])/sqrt(BEFF);
      SUPER    [2][block] = sqrt(swn2-SUPER    [1][block]*SUPER    [1][block])/sqrt(BEFF);
      for(int g=0;g<NG;g++){
        GOFR[2][block][g] = sqrt(gor2[g]-GOFR[1][block][g]*GOFR[1][block][g])/sqrt(BEFF);
      };
      for(int g=0;g<NK;g++){
        for(int t=0;t<TMAX;t++){
        SKX_RE[2][block][g][t] = sqrt(skr2[g][t]-SKX_RE[1][block][g][t]*SKX_RE[1][block][g][t])/sqrt(BEFF);
        SKX_IM[2][block][g][t] = sqrt(ski2[g][t]-SKX_IM[1][block][g][t]*SKX_IM[1][block][g][t])/sqrt(BEFF);
        };
      };
      for(int i=0;i<NX;i++){
        for(int j=0;j<NY;j++){
          RHOXY[2][block][i][j] = sqrt(rho2[i][j]-RHOXY[1][block][i][j]*RHOXY[1][block][i][j])/sqrt(BEFF);
          GOFXY[2][block][i][j] = sqrt(gxy2[i][j]-GOFXY[1][block][i][j]*GOFXY[1][block][i][j])/sqrt(BEFF);
        };
      };
  };
  return;
};

void System::Run(){
  double sect   = 0.0;
  double sectb  = 0.0;
  double scount = 0.0;
  if(Equi()>0)
    std::cout << "\nBeginning equilibration\n";
  for(int e=0;e<Equi();e++){
    for(int s=0;s<Steps();s++){
      for(int a=0;a<Sparse();a++)
        Move();
    };
    std::cout << "Equilibration block " << e+1 << "/" << Equi() << "\tdone\n";
  };
  if(Equi()>0)
    std::cout << "Equilibration done!\n";
  for(int b=0;b<Blocks();b++){
    sectb = 0.0;
    for(int s=0;s<Steps();s++){
      for(int a=0;a<Sparse();a++)
        Move();
      Measure(b);      
      sect += double(Sector());
      sectb+= double(Sector());
      scount++;
    };
    Average(b);
    std::cout << "\nBlock " << b+1 << "/" << Blocks() << " completed!\n";
    int t0 = 0;
    if(Sector() == 1) t0 = IRA[0];
    std::cout << "N\t" << NLnk(t0) << "\tSECT\t" << sect/scount << std::endl;
    std::cout << "Sector                  \t" << sectb/double(STEPS) << "\t" << sect/scount << "\n";
    std::cout << "Kinetic_energy          \t" << KINETIC  [0][b] << "\t" << KINETIC  [1][b]
                                      << "\t" << KINETIC  [2][b] << "\t" << WEIGHT[b] << "\n";
    std::cout << "Potential_energy        \t" << POTENTIAL[0][b] << "\t" << POTENTIAL[1][b]
                                      << "\t" << POTENTIAL[2][b] << "\t" << WEIGHT[b] << "\n";
    std::cout << "Total_energy            \t" << ENERGY   [0][b] << "\t" << ENERGY   [1][b]
                                      << "\t" << ENERGY   [2][b] << "\t" << WEIGHT[b] << "\n";
    if(SUPER[1][b]>1.e-8)
    std::cout << "Superfluid_fraction     \t" << SUPER    [0][b] << "\t" << SUPER    [1][b]
                                      << "\t" << SUPER    [2][b] << "\t" << WEIGHT[b] << "\n";
    WriteOut("coord.dat","index.dat");
    PrintG(b);
    PrintS(b);
    PrintRho(b);
    PrintGxy(b);
    PrintWorm();
    PrintAcc();
  };
  std::cout << "\nACCEPTATIONS:\n\n";
  std::cout << "Advance\t" << TRY[0] << "\t" << ACC[0]/TRY[0]*100.0 << std::endl;
  std::cout << "Recede \t" << TRY[1] << "\t" << ACC[1]/TRY[1]*100.0 << std::endl;
  std::cout << "Open   \t" << TRY[2] << "\t" << ACC[2]/TRY[2]*100.0 << std::endl;
  std::cout << "Close  \t" << TRY[3] << "\t" << ACC[3]/TRY[3]*100.0 << std::endl;
  std::cout << "Swap   \t" << TRY[4] << "\t" << ACC[4]/TRY[4]*100.0 << std::endl;
  std::cout << "Wiggle \t" << TRY[5] << "\t" << ACC[5]/TRY[5]*100.0 << std::endl;
  std::cout << "Shift  \t" << TRY[6] << "\t" << ACC[6]/TRY[6]*100.0 << std::endl;
  std::cout << "\nAVERAGE SECTOR:   " << sect/scount << "\n\n";
  std::cout << "Worm simulation finished! \n\n";
};

void System::PrintAcc() const{
  std::ofstream acfile;
  acfile.open("accept.dat");
  acfile << "\nTried & accepted moves:\n\n";
  acfile << std::setprecision(5);
  acfile << "Advance\t" << TRY[0] << "\t" << ACC[0]/TRY[0]*100.0 << " %" << std::endl;
  acfile << "Recede \t" << TRY[1] << "\t" << ACC[1]/TRY[1]*100.0 << " %" << std::endl;
  acfile << "Open   \t" << TRY[2] << "\t" << ACC[2]/TRY[2]*100.0 << " %" << std::endl;
  acfile << "Close  \t" << TRY[3] << "\t" << ACC[3]/TRY[3]*100.0 << " %" << std::endl;
  acfile << "Swap   \t" << TRY[4] << "\t" << ACC[4]/TRY[4]*100.0 << " %" << std::endl;
  acfile << "Wiggle \t" << TRY[5] << "\t" << ACC[5]/TRY[5]*100.0 << " %" << std::endl;
  acfile << "Shift  \t" << TRY[6] << "\t" << ACC[6]/TRY[6]*100.0 << " %" << std::endl;
  acfile << std::endl;
  acfile.close();
  return;
};

void System::PrintSyst() const{
  std::cout << "\nSYSTEM DATA\n";
  std::cout << "N Part      " << NUM       << std::endl;
  std::cout << "Mass        " << MASS      << std::endl;
  std::cout << "Lambda      " << LAMBDA    << std::endl;
  std::cout << "Dimensions  " << DIM       << std::endl;
  std::cout << "Density     " << NUM/Vol() << std::endl;
  std::cout << "Volume      " << Vol()     << std::endl;
  std::cout << "Box         ";
  for(int i=0;i<DIM;i++)
    std::cout << L(i) << "    ";
  std::cout << std::endl;
  std::cout << "Rmax        " << RMAX      << std::endl;
  std::cout << "Ecut        " << ECUT      << std::endl;
  std::cout << "Timeslices  " << Slices()  << std::endl;
  std::cout << "Tau         " << Tau()     << std::endl;
  std::cout << "Temperature " << T()       << std::endl;
  if(Sector()==0)
    std::cout << "Sector      Z" << std::endl;
  else
    std::cout << "Sector      G" << std::endl;
  std::cout <<"\nSIMULATION DATA\n";
  std::cout << "Equilibration blocks " << Equi()   << std::endl;
  std::cout << "Computation   blocks " << Blocks() << std::endl;
  std::cout << "Sweeps/block         " << Steps()  << std::endl;
  std::cout << "Moves/sweep          " << Sparse() << std::endl;
  std::cout <<"\nMOVE PARAMETERS\n";
  std::cout << "C0       " << C0()*(Vol()*double(Slices()*MBar())) << std::endl;
  std::cout << "M worm   " <<  MBar() << std::endl;
  std::cout << "M bridge " <<  MBri() << std::endl;
  std::cout << "Shifts   ";
  for(int i=0;i<DIM;i++)
    std::cout << DX[i] << "    ";
  std::cout << std::endl;
  return;
};

void System::TestPot() const{
  int imax = SIDE[0]/DR*0.5;
  std::ofstream potf;
  potf.open("testpot");
  for(int i=1;i<imax+1;i++)
    potf << i*DR << "\t" << U0(i*DR) << std::endl;
  potf << SIDE[0]*0.5 << "\t" << U0(SIDE[0]*0.5) << std::endl;
  return;
};

void System::NormGofR(int block){
  double a;
  if(DIM==1) a=1.0;
  if(DIM==2) a=Pi;
  if(DIM==3) a=4./3.*Pi;
  for(int i=0;i<DIM;i++)
    a = a / SIDE[i];
  double pref = 0.5 * a * pow(NLnk(0),2.0) * Slices();
  double v0=0.0;
  double v1;
  double rr=0.0; 
  for(int i=0;i<NG-1;i++){
    rr+=DR;
    v1 = pow(rr,double(DIM));
    GOFR[0][block][i] = GOFR[0][block][i]/(pref*(v1-v0));
    v0 = v1;
  };
  return;
};
