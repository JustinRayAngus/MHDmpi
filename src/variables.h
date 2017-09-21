/***
 * 
 * variables function class header file
 *
 * This file should be included in each physics 
 * module file where the functions declared
 * below are defined
 *
***/

#ifndef variables_h
#define variables_h

//#include "domainGrid.h"
//#include "timeDomain.h"
//#include "HDF5dataFile.h"

using namespace std;

class variables
{
public:
   void initialize(const domainGrid&, const Json::Value&, HDF5dataFile&);
   void setdtSim(double&, const timeDomain&, const domainGrid&);
   void advance(const domainGrid&, const double&);
};


#endif
