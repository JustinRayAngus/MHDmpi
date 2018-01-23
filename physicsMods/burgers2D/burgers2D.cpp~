/***
 * 
 * physics module for 2D burgers equation
 *
 * Note: I original computed advection and diff fluxes
 *       separately when using TVD scheme, which I only 
 *       used for adv flux. However, I found that sometimes
 *       I needed a smaller time step than expected or else
 *       oscillations would occur.
 *
 *       When doing TVD, I now use the total flux and it works
 *       just fine
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

string advScheme0;   // advection differencing scheme
double K;            // diffusion coefficient
double countJRA;
int order;         // order in time (1 or 2)
int fluxDir;       // flux direction (+/-1)
vector<vector<double>> F0, F0old;    // function
vector<vector<double>> FluxLimL, FluxLimR;
vector<vector<double>> Flux, FluxR, FluxL;  // flux at cell-edges   


void computeFluxes(const domainGrid&, const int);
void setXminBoundary(vector<vector<double>>&, 
		     const double, const double);
void setXmaxBoundary(vector<vector<double>>&, 
		     const double, const double);


void Physics::initialize(const domainGrid& Xgrid, const Json::Value& root, 
                      HDF5dataFile& dataFile)
{
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   // alternative means for getting grid, opposed
   // to passing predefined instance (Xgrid) to 
   // functions
   //
   domainGrid* mesh = domainGrid::mesh;
   cout << "mesh->nXcc = " << (*mesh).nXcc << endl;
   cout << "mesh->nXcc = " << mesh->nXcc << endl;
   
   const int nXcc = Xgrid.Xcc.size();
   const int nXce = Xgrid.Xce.size();
   const int nZcc = Xgrid.Zcc.size();
   const int nZce = Xgrid.Zce.size();
   
   F0.assign(nXcc,vector<double>(nZcc));
   F0old.assign(nXcc,vector<double>(nZcc));
   FluxLimL.assign(nXce,vector<double>(nZcc));
   FluxLimR.assign(nXce,vector<double>(nZcc));
   Flux.assign(nXce,vector<double>(nZcc));
   FluxR.assign(nXce,vector<double>(nZcc));
   FluxL.assign(nXce,vector<double>(nZcc));
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      
      //   get advection scheme from input file
      //
      Json::Value advScheme  = Phys.get("advScheme",defValue);
      if(advScheme != defValue) advScheme0 = advScheme.asString();
      else advScheme0 = "U1";
      
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

      //   get diffusion coefficent from input file
      //
      Json::Value KVal = Phys.get("diffC",defValue);
      if(KVal != defValue) K = KVal.asDouble();
      else K = 0.0;

      if(procID==0) cout << "diffusion coefficent = " << K << endl;
      if(K < 0.0) {
         printf("ERROR: diffusion coefficient can't be < 0\n");
         exit (EXIT_FAILURE);
      }

      //   get order in time from input file
      //
      Json::Value orderVal = Phys.get("order",defValue);
      if(orderVal != defValue) order = orderVal.asInt();
      else order = 2;
      
      if(procID==0) cout << "order in time is " << order << endl;
      if(order != 1 && order != 2) {
         printf("ERROR: order in time can only be 1 or 2\n");
         exit (EXIT_FAILURE);
      }
      
      //   get flux direction from input file
      //
      Json::Value fluxDirVal = Phys.get("fluxDir",defValue);
      if(fluxDirVal != defValue) fluxDir = fluxDirVal.asInt();
      else fluxDir = 1;

      if(fluxDir==1) {
         if(procID==0) cout << "flux direction is in +X direction" << endl;
      } 
      else if(fluxDir==-1) {
         if(procID==0) cout << "flux direction is in -X direction" << endl;
      }
      else {
         printf("ERROR: fluxDir can only be +/-1\n");
         exit (EXIT_FAILURE);
      }

   }
   else {
      cout << "value for key \"Physics\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
  

   Xgrid.setInitialProfile(F0,Phys);
   if(procID==0) setXminBoundary(F0, 0.0, 0.0);   
   if(procID==numProcs-1) setXmaxBoundary(F0, 0.0, 0.0);   
   Xgrid.communicate(F0);
   F0old  = F0;
   computeFluxes(Xgrid,1); // inital calculation before add to output   
  

   dataFile.add(K, "K", 0);
   dataFile.add(F0, "F0", 1); 
   dataFile.add(FluxLimL, "FluxLimL", 1);  
   dataFile.add(FluxLimR, "FluxLimR", 1);  
   dataFile.add(Flux, "Flux", 1);  
   dataFile.add(FluxR, "FluxR", 1);
   dataFile.add(FluxL, "FluxL", 1);

   cout << "F0.size() = " << F0.size() << endl;
   cout << "F0[0].size() = " << F0[0].size() << endl;
} // end Physics.initilize


void Physics::advance(const domainGrid& Xgrid, const double dt)
{
   const int nMax = F0.size();
   const int nXg = Xgrid.nXg;
   const int nZcc = Xgrid.nZcc;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   countJRA = countJRA + 1;
   //cout << " countJRA = " << countJRA << endl;
   //cout << "F0.size() = " << F0.size() << endl;
   //cout << "F0[0].size() = " << F0[0].size() << endl;

   // Explicit forward advance from n to n+1/2 using Flux(n)
   // (calc of Flux(t=0) is done during initilization)
   //
   for (auto i=nXg; i<nMax-nXg; i++) {
      for (auto j=0; j<nZcc; j++) {
         F0[i][j] = F0old[i][j] - dt/2.0*(Flux[i][j]-Flux[i-1][j])/Xgrid.dX;
      }
   }

   // apply boundary conditions, communicate, compute Flux
   //
   if(procID==0) setXminBoundary(F0, 0.0, 0.0);   
   if(procID==numProcs-1) setXmaxBoundary(F0, 0.0, 0.0);   
   Xgrid.communicate(F0);
   computeFluxes(Xgrid,order); // compute RHS using F0(n+1/2)


   // Explicit forward advance from n to n+1 using Flux(n+1/2)
   //
   for (auto i=nXg; i<nMax-nXg; i++) {
      for (auto j=0; j<nZcc; j++) {
         F0[i][j] = F0old[i][j] - dt*(Flux[i][j]-Flux[i-1][j])/Xgrid.dX;
      }
   }
   

   // apply boundary conditions, communicate, compute Flux
   //
   if(procID==0) setXminBoundary(F0, 0.0, 0.0);   
   if(procID==numProcs-1) setXmaxBoundary(F0, 0.0, 0.0);   
   Xgrid.communicate(F0);
   computeFluxes(Xgrid,1);
   

   //  update F0old
   //
   F0old = F0;

} // end Physics.advance


void computeFluxes(const domainGrid& Xgrid, const int order)
{
   // compute Flux = F0^2/2 - K*dF0/dX
   //              = FluxAdv + FluxDif
   //
   // FluxAdv = (FluxR+FluxL)/2 is computed using upwinding schemes
   // FluxDif is computed using standard centered scheme


   const int nCE = Flux.size();
   const int nXcc = F0.size();
   const int nZcc = F0[0].size();
   
   vector<vector<double>> Cspeed, FluxAdvCC, FluxDifCC, FluxAdv, FluxDif;
   Cspeed.assign(nXcc,vector<double>(nZcc));
   FluxAdvCC.assign(nXcc,vector<double>(nZcc));
   FluxDifCC.assign(nXcc,vector<double>(nZcc));
   FluxAdv.assign(nCE,vector<double>(nZcc));
   FluxDif.assign(nCE,vector<double>(nZcc));
   

   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeed = abs(F0); // adv flux jacobian
   FluxAdvCC = double(fluxDir)*F0*F0*0.5;


   // compute diffusive flux using 
   // standard centered scheme
   //
   Xgrid.DDX(FluxDifCC,F0);
   FluxDifCC = -K*FluxDifCC;
   Xgrid.communicate(FluxDifCC);

   Xgrid.DDX(FluxDif,F0);
   FluxDif = -K*FluxDif;
   Xgrid.communicate(FluxDif);
   //FluxDif = DDX(F0,Xgrid.dX);


   // compute advective flux using
   // specified scheme from input file
   //
   if(advScheme0 == "TVD") {
      //Xgrid.computeFluxTVD(FluxAdv,FluxL,FluxR,FluxLimL,FluxLimR,
      //                     FluxAdvCC,Cspeed,F0,0,order);
      Xgrid.computeFluxTVD(Flux,FluxL,FluxR,FluxLimL,FluxLimR,
                           FluxAdvCC+FluxDifCC,Cspeed,F0,0,order);
   }      
   else {
      Xgrid.InterpToCellEdges(FluxAdv,FluxAdvCC,Cspeed,advScheme0,0);
      Flux = FluxAdv + FluxDif;
   } 
   Xgrid.communicate(Flux);
   Xgrid.communicate(FluxL);
   Xgrid.communicate(FluxR);


} // end computeFluxes



void setXminBoundary(vector<vector<double>>& var,
		     const double C0, const double C1)
{
   //const int thisnX = var.size();
   const int thisnZ = var[0].size();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift=mesh->nXg;
  

   for (auto i=0; i<ishift; i++) {
      for (auto j=0; j<thisnZ; j++) {
         var[ishift-i-1][j] = C0 + C1*var[ishift+i][j];
      }
   }

}

void setXmaxBoundary(vector<vector<double>>& var,
		     const double C0, const double C1)
{
   const int thisnX = var.size();
   const int thisnZ = var[0].size();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift=thisnX-mesh->nXg;

   for (auto i=ishift; i<thisnX; i++) {
      for (auto j=0; j<thisnZ; j++) {
         var[i][j] = C0 + C1*var[2*ishift-i-1][j];
      }
   }   

}


void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   double Umax;
   Umax = max(abs(F0));
   //cout << "Umax = " << Umax << endl;

   const double dX = Xgrid.dX;
   double dtmaxDif = 0.5*dX*dX/K;
   double dtmaxAdv = dX/Umax;
   double dtmax = min(dtmaxDif,dtmaxAdv);
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);
   if(procID==0) {
      cout << "max stable time step is " << dtmax << endl;
      cout << endl; 
   }
}

