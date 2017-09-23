/***
 * 
 * time domain class
 *
***/

#ifndef timeDomain_h
#define timeDomain_h

#include "json/json.h"
//#include "vectorMath.h"   // causes duplicate symbol error!!!
//#include "HDF5dataFile.h" // causes duplicate symbol error!!!
//#include "mpi.h"


using namespace std;

class timeDomain
{
 
public:
  double dtOut, tOutSteps;   // Output intervals and number of steps
  double tmax;           // Max time
  double dtSim;              // Simulation time-step
  double dtFrac;             // dtSim = dtmax/dtFrac 
  vector<double> tOutVec;        // vector of output times
  double tOut;           // current output time 
  
  void initialize(const Json::Value&);
  void updatetOut(double thistOut);

};


#endif
