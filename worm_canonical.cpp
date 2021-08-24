#include "System.hpp"

int main(int argc, char** argv){
  System ToBeSim(argv[1]);
  ToBeSim.TestPot();
  ToBeSim.Run();
  return 0;
};
