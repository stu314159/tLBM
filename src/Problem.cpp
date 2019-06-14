#include "Problem.h"
#include <fstream>
#include <stdexcept>

void Problem::load_input(){
  std::ifstream input_params("params.lbm",std::ios::in);
  if(!input_params.is_open()){
    throw std::runtime_error("Could not open params.lbm!");
  }
  input_params >> latticeType;
  input_params >> dynamics;
  input_params >> numTs;
  input_params >> tsRepFreq;
  input_params >> warmupTs;
  input_params >> plotFreq;
  input_params >> cs;
  input_params >> rhoLBM;
  input_params >> uLBM; 
  input_params >> omega;
  input_params >> nx;
  input_params >> ny;
  input_params >> nz;
  input_params >> restartFlag;
  input_params >> timeAvgFlag;
  input_params >> lx_p;
  input_params >> ly_p;
  input_params >> lz_p;
  input_params >> tConvFact;
  input_params >> lConvFact;
  input_params >> pConvFact;
  input_params >> pRefIdx;
  input_params >> ssDataFlag;
  input_params.close();


}
