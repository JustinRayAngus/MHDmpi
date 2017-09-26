/***
 * 
 * physics module for sod shock test
 *
***/

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
#include <vector>
#include <cmath>

#include "json/json.h"
#include "vectorMath.h"
#include "domainGrid.h"
#include "timeDomain.h"
#include "HDF5dataFile.h"

#include "Physics.h"


using namespace std;

string advScheme0;    // advection differencing scheme
double gamma0;        // adiabatic coefficient
vector<double> N, M, E;   // function
vector<double> F0, Cs, V, P;      // function
vector<double> F0old, Nold, Mold, Eold;
vector<double> FluxRatio, FluxLim;
vector<double> FluxR, FluxL;  // flux at cell-edges   
vector<double> FluxN, FluxM, FluxE;

void computeFluxes(const domainGrid&);
void setXminBoundary(vector<double>&, const double&);
void setXmaxBoundary(vector<double>&, const double&);

// alternative means for getting grid, opposed
// to passing predefined instance (Xgrid) to 
// functions
//
//domainGrid* mesh = domainGrid::mesh;
//const int nXg = mesh->nXg;
//cout << "nXg = " << mesh->nXg << endl;
//cout << "mesh->nXcc = " << (*mesh).nXcc << endl;
//cout << "mesh->nXcc = " << mesh->nXcc << endl;


void Physics::initialize(const domainGrid& Xgrid, const Json::Value& root, 
                      HDF5dataFile& dataFile)
{
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   //domainGrid* mesh = domainGrid::mesh;
   
   const int nXcc = Xgrid.Xcc.size();
   F0.assign(nXcc,0.0);
   F0old.assign(nXcc,0.0);
   N.assign(nXcc,0.0);
   Nold.assign(nXcc,0.0);
   M.assign(nXcc,0.0);
   Mold.assign(nXcc,0.0);
   E.assign(nXcc,0.0);
   Eold.assign(nXcc,0.0);
   P.assign(nXcc,0.0);
   Cs.assign(nXcc,0.0);
   V.assign(nXcc,0.0);
   //
   const int nXce = Xgrid.Xce.size();
   FluxRatio.assign(nXce,0.0);
   FluxLim.assign(nXce,0.0);
   FluxN.assign(nXce,0.0);
   FluxM.assign(nXce,0.0);
   FluxE.assign(nXce,0.0);
   FluxR.assign(nXce,0.0);
   FluxL.assign(nXce,0.0);
   //
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      Json::Value advScheme = Phys.get("advScheme",defValue);
      Json::Value gammaVal = Phys.get("gammaC",defValue);
      if(advScheme == defValue || gammaVal == defValue) {
         cout << "ERROR: advScheme or gamma is " << endl;
         cout << "not declared in input file" << endl;
         exit (EXIT_FAILURE);
      } 
      
      advScheme0 = advScheme.asString();
      if(advScheme0=="C2" || advScheme0=="U1" || 
         advScheme0=="QUICK" || advScheme0=="TVD") {
         if(procID==0) {
            cout << "advection diff/interp scheme is " << advScheme0 << endl;
         }
      }
      else {
         cout << "advection scheme " << advScheme0 << " is not valid " << endl;
         cout << "valid types are C2, U1, QUICK, and TVD " << endl;
         exit (EXIT_FAILURE);
      }

      gamma0 = gammaVal.asDouble();
      if(procID==0) cout << "adiabatic coefficent = " << gamma0 << endl;
      if(gamma0 < 1.0) {
         printf("ERROR: adiabatic coefficient can't be < 1\n");
         exit (EXIT_FAILURE);
      }

   }
   else {
      cout << "value for key \"Physics\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
  

   //   get initial profiles for variables
   //
   const Json::Value Nvar = Phys.get("N",defValue);
   if(Nvar.isObject()) { 
      Xgrid.setInitialProfile(N,Nvar);
      if(procID==0) setXminBoundary(N, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(N, 0.125);   
      Xgrid.communicate(N);
      Nold  = N;
      //computeFluxes(Xgrid); // inital calculation before add to output   
   } else {
      cout << "value for Physics variable \"N\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
 
   const Json::Value Vvar = Phys.get("V",defValue);
   if(Vvar.isObject()) { 
      Xgrid.setInitialProfile(V,Vvar);
      if(procID==0) setXminBoundary(V, 0.0);   
      if(procID==numProcs-1) setXmaxBoundary(V, 0.0);   
      Xgrid.communicate(V);
      M = N*V;
      Mold  = M;
   } else {
      cout << "value for Physics variable \"V\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   const Json::Value Pvar = Phys.get("P",defValue);
   if(Pvar.isObject()) { 
      Xgrid.setInitialProfile(P,Pvar);
      if(procID==0) setXminBoundary(P, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(P, 0.1);   
      Xgrid.communicate(P);
      E = P/(gamma0-1.0) + 0.5*M*V;
      Eold  = E;
   } else {
      cout << "value for Physics variable \"P\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   Cs = sqrt(gamma0*P/N);

   computeFluxes(Xgrid); // inital calculation before add to output   
  
   dataFile.add(N, "N", 1);  // density 
   dataFile.add(M, "M", 1);  // momentum density 
   dataFile.add(E, "E", 1);  // energy density
   dataFile.add(P, "P", 1);  // pressure
   dataFile.add(V, "V", 1);  // velocity
   dataFile.add(Cs,"Cs",1);  // sound speed
   dataFile.add(gamma0,"gamma0",0); 
   //
   dataFile.add(FluxN, "FluxN", 1);  
   dataFile.add(FluxM, "FluxM", 1);  
   dataFile.add(FluxE, "FluxE", 1);  
   //
   dataFile.add(FluxRatio, "FluxRatio", 1);  
   dataFile.add(FluxLim, "FluxLim", 1);  
   dataFile.add(FluxR, "FluxR", 1);
   dataFile.add(FluxL, "FluxL", 1);

} // end Physics.initilize


void Physics::advance(const domainGrid& Xgrid, const double& dt)
{
   const int nMax = N.size();
   const int nXg = Xgrid.nXg;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   //domainGrid* mesh = domainGrid::mesh;

   // Explicit forward advance from n to n+1/2 using Flux(n)
   // (calc of Flux(t=0) is done during initilization)
   //
   for (auto i=nXg; i<nMax-nXg; i++) {
      N.at(i) = Nold.at(i) - dt/2.0*(FluxN.at(i)-FluxN.at(i-1))/Xgrid.dX;
      M.at(i) = Mold.at(i) - dt/2.0*(FluxM.at(i)-FluxM.at(i-1))/Xgrid.dX;
      E.at(i) = Eold.at(i) - dt/2.0*(FluxE.at(i)-FluxE.at(i-1))/Xgrid.dX;
   }


   // apply boundary conditions, communicate, compute Flux
   //
   if(procID==0) {
      setXminBoundary(N, 1.0);   
      setXminBoundary(M, 0.0);   
      setXminBoundary(E, 1.0/(gamma0-1.0));   
   }
   if(procID==numProcs-1) {
      setXmaxBoundary(N, 0.125);   
      setXmaxBoundary(M, 0.0);   
      setXmaxBoundary(E, 0.1/(gamma0-1.0));   
   }
   Xgrid.communicate(N);
   Xgrid.communicate(M);
   Xgrid.communicate(E);
   computeFluxes(Xgrid); // compute RHS using N(n+1/2)


   // Explicit forward advance from n to n+1 using Flux(n+1/2)
   //
   for (auto i=nXg; i<nMax-nXg; i++) {
      N.at(i) = Nold.at(i) - dt*(FluxN.at(i)-FluxN.at(i-1))/Xgrid.dX;
      M.at(i) = Mold.at(i) - dt*(FluxM.at(i)-FluxM.at(i-1))/Xgrid.dX;
      E.at(i) = Eold.at(i) - dt*(FluxE.at(i)-FluxE.at(i-1))/Xgrid.dX;
   }
   

   // apply boundary conditions, communicate, compute Flux
   //
   if(procID==0) {
      setXminBoundary(N, 1.0);   
      setXminBoundary(M, 0.0);   
      setXminBoundary(E, 1.0/(gamma0-1.0));   
   }
   if(procID==numProcs-1) {
      setXmaxBoundary(N, 0.125);   
      setXmaxBoundary(M, 0.0);   
      setXmaxBoundary(E, 0.1/(gamma0-1.0));   
   }
   Xgrid.communicate(N);
   Xgrid.communicate(M);
   Xgrid.communicate(E);
   computeFluxes(Xgrid);
   

   //  update Nold
   //
   Nold = N;
   Mold = M;
   Eold = E;

} // end Physics.advance


void computeFluxes(const domainGrid& Xgrid)
{
   // FluxN = M
   // FluxM = M^2/N + P
   //

   const int nCE = FluxN.size();
   const int nCC = N.size();
   vector<double> Cspeed, FluxNcc, FluxMcc, FluxEcc;
   vector<double> FluxP;
   FluxNcc.assign(nCC,0.0);
   FluxMcc.assign(nCC,0.0);
   FluxEcc.assign(nCC,0.0);
   FluxP.assign(nCE,0.0);
   

   //  define derived variables
   //
   V  = M/N;
   P  = (E-0.5*M*M/N)*(gamma0-1);
   Cs = sqrt(gamma0*P/N);


   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeed = V + Cs; // adv flux jacobian
   FluxNcc = M;
   FluxMcc = M*M/N + P;
   FluxEcc = V*(E+P);
   

   // compute advective flux using
   // specified scheme from input file
   //
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxNcc,Cspeed,N);
      Xgrid.computeFluxTVD(FluxM,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxMcc,Cspeed,M);
      Xgrid.computeFluxTVD(FluxE,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxEcc,Cspeed,E);
   }      
   else {
      Xgrid.InterpToCellEdges(FluxN,FluxNcc,Cspeed,advScheme0);
      Xgrid.InterpToCellEdges(FluxM,FluxMcc,Cspeed,advScheme0);
      Xgrid.InterpToCellEdges(FluxE,FluxEcc,Cspeed,advScheme0);
   } 


} // end computeFluxes


void setXminBoundary(vector<double>& var, const double& C)
{
   
   domainGrid* mesh = domainGrid::mesh;
   //F0.front() = 2.0*C-F0.at(1); 
   //F0.front() = C; 
   for (auto i=0; i<mesh->nXg; i++) {
      var.at(i) = C;
   }
   
}


void setXmaxBoundary(vector<double>& var, const double& C)
{
   
   domainGrid* mesh = domainGrid::mesh;
   //F0[nXsub+1] = 2.0*C-F0[nXsub]; 
   //F0[nXsub+1] = C;
   for (auto i=mesh->nXcc-mesh->nXg; i<mesh->nXcc; i++) {
      var.at(i) = C;
   }
   //var.back() = C;
   //cout << "var.size() = " var.size() << endl; 
      
}


void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   vector<double> Cchar;
   double Cmax;
   Cchar = abs(V)+Cs;
   Cmax = max(Cchar);
   //cout << "Cmax = " << Cmax << endl;

   const double dX = Xgrid.dX;
   double dtmax = dX/Cmax;
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);
   if(procID==0) {
      cout << endl; 
      cout << "max stable time step is " << dtmax << endl;
      cout << endl; 
   }
}

