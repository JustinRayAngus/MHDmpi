/***
 * 
 * time domain class
 *
***/

#ifndef timeDomain_h
#define timeDomain_h

#include "EEDF.h"
#include "domainGrid.h"

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
  void setdtSim(const EEDF&, const domainGrid&, const double&);

};

void timeDomain::initialize(const Json::Value& root)
{
   const Json::Value defValue; // used for default reference
   const Json::Value Time = root.get("Time",defValue);
   if(Time.isObject()) {
      printf("Initializing time domain ...\n");
      Json::Value dtOutVal      = Time.get("dtOut",defValue);
      Json::Value tOutStepsVal  = Time.get("tOutSteps",defValue);
      Json::Value dtFracVal     = Time.get("dtFrac",defValue);
      if(dtOutVal == defValue || tOutStepsVal == defValue || dtFracVal == defValue) {
         printf("ERROR: default 'Time' variables not declared in input file\n");
         exit (EXIT_FAILURE);
      }
      dtOut = dtOutVal.asDouble();
      tOutSteps = tOutStepsVal.asDouble();
      dtFrac = dtFracVal.asDouble();
      tmax = dtOut*tOutSteps;
      //
      cout << "tmax = " << dtOut*tOutSteps << endl;
      cout << "tOut intervals = " << dtOut << endl;
      cout << "dtFrac = " << dtFrac << endl;
   }
   else {
      cout << "value for key \"Time\" is not object type !" << endl;
   }
  
   tOutVec.assign(tOutSteps+1,0.0);
   for(auto n=0; n<tOutSteps+1; n++) {
      tOutVec[n] = n*dtOut;  // vector of output times for references
      //cout << tOutVec[n] << endl;
   }
   
   tOut = 0.0; // first output is always at t=0
   cout << endl;
}

void timeDomain::updatetOut(const double thistOut)
{
   tOut = thistOut;
}

void timeDomain::setdtSim(const EEDF& eedf, const domainGrid& Xgrid, const double& K)
{
   const double dX = Xgrid.dX;
   double dtmax = 0.5*dX*dX/K;
   dtSim = min(dtmax/dtFrac,dtOut);
   cout << "max stable time step is " << dtmax << endl;
}

#endif
