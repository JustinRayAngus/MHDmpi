/***
 * 
 * variables function class header file
 *
 * This file should be included in each physics 
 * module file where the below functions are
 * defined in detail
 *
***/

#ifndef variables_h
#define variables_h

#include "domainGrid.h"

using namespace std;

class variables
{
public:
   void initialize(const domainGrid&, const Json::Value&, HDF5dataFile&);
   void setdtSim(double&, const timeDomain&, const domainGrid&);
   void advanceF0(const domainGrid&, const double&);
};

/*
void variables::initialize(const domainGrid& Xgrid, const Json::Value& root, 
                      HDF5dataFile& dataFile)
{
   cout << "HELLO VARIABLE CLASS" << endl;
}
*/

#endif
