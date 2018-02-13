/***
 * 
 * time domain class
 *
***/


#include "timeDomain.h"
#include "mpi.h"
//#include "vectorMath.h"   // causes duplicate symbol error!!!
//#include "HDF5dataFile.h" // causes duplicate symbol error!!!

using namespace std;


void timeDomain::initialize(const Json::Value& root)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);

   const Json::Value defValue; // used for default reference
   const Json::Value Time = root.get("Time",defValue);
   if(Time.isObject()) {
      if(procID==0) printf("\nInitializing time domain ...\n");
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
      if(procID==0) {
         cout << "tmax = " << dtOut*tOutSteps << endl;
         cout << "tOut intervals = " << dtOut << endl;
         cout << "dtFrac = " << dtFrac << endl;
         cout << endl;
      }
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
}

void timeDomain::updatetOut(const double thistOut)
{
   tOut = thistOut;
}


