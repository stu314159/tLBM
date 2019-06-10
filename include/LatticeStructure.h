#ifndef LATTICESTRUCTURE_H
#define LATTICESTRUCTURE_H

#include <vector>

class LatticeStructure
{
public:
  LatticeStructure();
  ~LatticeStructure();
  

private:
  unsigned int numSpd;
  std::vector<int> ex;
  std::vector<int> ey;
  std::vector<int> ez;
  std::vector<int> bbSpd;
  std::vector<double> w;

};

#endif
