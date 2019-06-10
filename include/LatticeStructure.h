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
  vector<int> ex;
  vector<int> ey;
  vector<int> ez;
  vector<int> bbSpd;
  vector<double> w;

}

#endif
