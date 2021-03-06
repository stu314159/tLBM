#ifndef PROBLEM_H
#define PROBLEM_H

// an object whose function is to read and hold the problem
// parameters for a tLBM simulation.

#include <string>
#include "TLBM_definitions.h"

class Problem{

public:

  void load_input();
  
  // really just a structure, so make all members public

  std::string latticeType;  // [ "D3Q15" | "D3Q19" | "D3Q27" ]
  int dynamics; // [1 = LBGK | 2 = RBGK | 3 = MRT]
  int numTs;
  int tsRepFreq;
  int warmupTs; // number of timesteps to run prior to writing data to disk
  int plotFreq; // timesteps between data writes
  real cs;  // turbulence parameter
  real rhoLBM;
  real uLBM;
  real omega;
  int nx;  // lattice size in X,Y, and Z directions
  int ny;
  int nz;
  int restartFlag; // [0 = no re-start | 1 = re-start]
  int timeAvgFlag; // [0 = no time average | 1 = time average]
  real lx_p;  // physical lenght of domain in each dimension
  real ly_p;
  real lz_p;
  real tConvFact;  // conversion factors for time, length, and pressure
  real lConvFact;
  real pConvFact;
  int pRefIdx;   // node index for reference pressure
  int ssDataFlag;  // [0 = do not take subset data | 1 = take subset data ]


  private:

};
#endif
