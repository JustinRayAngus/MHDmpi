#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <typeinfo>
#include <algorithm>
#include <mpi.h>
#include <assert.h>

#include "json/json.h"
#include "vectorMath.h"
#include "HDF5dataFile.h"
#include "domainGrid.h"
#include "timeDomain.h"
#include "Physics.h"
//#include "physicsMods/burgers1D.cpp" // physics module
//#include "burgers1D.cpp" // physics module


using namespace std;

//static int procID, numProcs;

#ifndef H5_NO_NAMESPACE
   using namespace H5;
#endif


// have to give initial definition to 
// static variables before main()
// (How does BOUT++ do this?)
//
//domainGrid defGrid;
//domainGrid* domainGrid::mesh = &defGrid;
domainGrid* domainGrid::mesh = NULL;
timeDomain* timeDomain::tmesh = NULL;


int main(int argc, char** argv) {   

   // start MPI stuff
   //
   int procID, numProcs;
   //double wtime;

   //MPI_Status status;
   MPI_Init (&argc, &argv);
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   double time_start = MPI_Wtime();
   if(procID==0) {
      cout << endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl;
      cout << "" << endl;
      cout << "Initiating simulation" << endl;
      cout << "mpi job with " << numProcs << " processors" << endl;
      cout << endl;
   }


   // Parse the specified json input file
   //
   const string inputFile = "./input.json";
   Json::Value inputRoot; // will contain root value after parsing
   Json::Reader reader; 
   ifstream ifile(inputFile);
   //bool isJsonOK = (ifile != NULL && reader.parse(ifile, inputRoot));
   bool isJsonOK = (reader.parse(ifile, inputRoot));
   if(isJsonOK) {
      if(procID==0) {
         cout << "Input file " << inputFile << " parsed successfully" << endl;   
      }
   }
   else {
      cout << "ERROR: json input file not found or cannot be parsed due to errors" << endl;
      exit (EXIT_FAILURE);
   }

   // set output file in HD5FdataFile class
   //
   HDF5dataFile dataFile;
   const string outputFile = "output" + to_string(procID) + ".h5";
   dataFile.setOutputFile(outputFile);
   //double procID2 = (double)procID;
   dataFile.add(procID,"procID",0);

   /*
   // the following are for testing appending matrix to data file
   //
   vector<vector<double>> Matrix;
   Matrix.resize(3,vector<double>(11,0.0));
   dataFile.add(Matrix,"Matrix",1);
   */


   // initialize spatial grid and time domain
   //
   domainGrid Xgrid;
   domainGrid::mesh = &Xgrid;
   //Xgrid.mesh = &Xgrid;
   Xgrid.initialize(inputRoot);
   dataFile.add(Xgrid.Xcc, "Xcc", 0); 
   dataFile.add(Xgrid.Xce, "Xce", 0); 
   
   //
   timeDomain tDom;
   timeDomain::tmesh = &tDom;
   tDom.initialize(inputRoot);
   dataFile.add(tDom.tOut, "tout", 1); // actual output time   

   // initialize Physics object
   //  
   Physics phys;
   phys.initialize(Xgrid, inputRoot, dataFile);   
   //initialize(Xgrid, inputRoot, dataFile);   


   // set initial time-step
   //
   double dtSim;
   phys.setdtSim(dtSim, tDom, Xgrid);
   tDom.dtSim = dtSim;
   if(procID==0) {
      cout << "Initial simulation time step: " << dtSim << endl << endl;
   }
   double thist=0;
   tDom.tSim = thist;
   //timeDomain::ptSim = &thist;

   int thistOutInt = 1;
   vector<double> F0m(Xgrid.nX,0.0), errorVec(Xgrid.nX,0.0);

  
   // advance variables in time
   //
   while(thist<tDom.tmax) {
      thist = thist + dtSim; // new time at end of this time step
      tDom.tSim = thist;

      // advance variables from n => n+1
      //
      phys.advance(Xgrid, dtSim);

      // check if thist is an output time
      //
      if(thist >= tDom.tOutVec[thistOutInt]) {
         tDom.updatetOut(thist);
         dataFile.writeAll(); // append extendable outputs
         thistOutInt = thistOutInt+1;
         phys.setdtSim(dtSim, tDom, Xgrid);
         if(procID==0) {
            cout << "Output variables dumped at t = " 
                 << thist << " units?" << endl;
            //cout << "Simulation time step = " << dtSim << endl;
         }
      }
   }

   if(procID==0) {
      double time_end = MPI_Wtime();
      cout << endl << "Final simulation time step = " << dtSim << endl;
      cout << endl << "Ending simulation: wall time = " 
           << time_end-time_start << endl;
      cout << endl << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%" << endl << endl;
   }

   MPI_Finalize(); 
   return 0;
}


