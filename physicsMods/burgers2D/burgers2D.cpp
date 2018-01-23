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
int order;         // order in time (1 or 2)
int XfluxDir,ZfluxDir;       // flux direction (+/-1)
vector<vector<double>> F0, F0old;    // function
vector<vector<double>> FluxLimL_x, FluxLimR_x;
vector<vector<double>> Flux_x, FluxR_x, FluxL_x;  // flux at cell-edges   
vector<vector<double>> FluxLimL_z, FluxLimR_z;
vector<vector<double>> Flux_z, FluxR_z, FluxL_z;  // flux at cell-edges   


void computeFluxes(const domainGrid&, const int);
void setXminBoundary(vector<vector<double>>&, 
		     const double, const double);
void setXmaxBoundary(vector<vector<double>>&, 
		     const double, const double);
void setZminBoundary(vector<vector<double>>&, 
		     const double, const double);
void setZmaxBoundary(vector<vector<double>>&, 
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
   //
   FluxLimL_x.assign(nXce,vector<double>(nZcc));
   FluxLimR_x.assign(nXce,vector<double>(nZcc));
   Flux_x.assign(nXce,vector<double>(nZcc));
   FluxR_x.assign(nXce,vector<double>(nZcc));
   FluxL_x.assign(nXce,vector<double>(nZcc));
   //
   FluxLimL_z.assign(nXcc,vector<double>(nZce));
   FluxLimR_z.assign(nXcc,vector<double>(nZce));
   Flux_z.assign(nXcc,vector<double>(nZce,0.0));
   FluxR_z.assign(nXcc,vector<double>(nZce));
   FluxL_z.assign(nXcc,vector<double>(nZce));

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
      Json::Value XfluxDirVal = Phys.get("XfluxDir",defValue);
      ZfluxDir = 1;
      if(XfluxDirVal != defValue) XfluxDir = XfluxDirVal.asInt();
      else XfluxDir = 1;

      if(XfluxDir==1) {
         if(procID==0) cout << "flux direction is in +X direction" << endl;
      } 
      else if(XfluxDir==-1) {
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
   setZminBoundary(F0, 0.0, 0.0);  
   setZmaxBoundary(F0, 0.0, 0.0);  
   Xgrid.communicate(F0);
   F0old  = F0;
   computeFluxes(Xgrid,1); // inital calculation before add to output   
  

   dataFile.add(K, "K", 0);
   dataFile.add(F0, "F0", 1); 
   dataFile.add(FluxLimL_x, "FluxLimL_x", 1);  
   dataFile.add(FluxLimR_x, "FluxLimR_x", 1);  
   dataFile.add(Flux_x, "Flux_x", 1);  
   dataFile.add(FluxR_x, "FluxR_x", 1);
   dataFile.add(FluxL_x, "FluxL_x", 1);
   //
   dataFile.add(FluxLimL_z, "FluxLimL_z", 1);  
   dataFile.add(FluxLimR_z, "FluxLimR_z", 1);  
   dataFile.add(Flux_z, "Flux_z", 1);  
   dataFile.add(FluxR_z, "FluxR_z", 1);
   dataFile.add(FluxL_z, "FluxL_z", 1);


} // end Physics.initilize


void Physics::advance(const domainGrid& Xgrid, const double dt)
{
   const int nMax = F0.size();
   const int nXg = Xgrid.nXg;
   const int nZcc = Xgrid.nZcc;
   const int nZg = Xgrid.nZg;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   // Explicit forward advance from n to n+1/2 using Flux(n)
   // (calc of Flux(t=0) is done during initilization)
   //
   for (auto i=nXg; i<nMax-nXg; i++) {
      for (auto j=nZg; j<nZcc-nZg; j++) {
         F0[i][j] = F0old[i][j] - dt/2.0*(Flux_x[i][j]-Flux_x[i-1][j])/Xgrid.dX
                                - dt/2.0*(Flux_z[i][j]-Flux_z[i][j-1])/Xgrid.dZ;
      }
   }

   // apply boundary conditions, communicate, compute Flux
   //
   if(procID==0) setXminBoundary(F0, 0.0, 0.0);   
   if(procID==numProcs-1) setXmaxBoundary(F0, 0.0, 0.0);   
   setZminBoundary(F0, 0.0, 0.0);
   setZmaxBoundary(F0, 0.0, 0.0);
   Xgrid.communicate(F0);
   computeFluxes(Xgrid,order); // compute RHS using F0(n+1/2)


   // Explicit forward advance from n to n+1 using Flux(n+1/2)
   //
   for (auto i=nXg; i<nMax-nXg; i++) {
      for (auto j=nZg; j<nZcc-nZg; j++) {
         F0[i][j] = F0old[i][j] - dt*(Flux_x[i][j]-Flux_x[i-1][j])/Xgrid.dX
                                - dt*(Flux_z[i][j]-Flux_z[i][j-1])/Xgrid.dZ;
      }
   }
   

   // apply boundary conditions, communicate, compute Flux
   //
   if(procID==0) setXminBoundary(F0, 0.0, 0.0);   
   if(procID==numProcs-1) setXmaxBoundary(F0, 0.0, 0.0);   
   setZminBoundary(F0, 0.0, 0.0);
   setZmaxBoundary(F0, 0.0, 0.0);
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


   const int nXce = Flux_x.size();
   const int nZce = Flux_z[0].size();
   const int nXcc = F0.size();
   const int nZcc = F0[0].size();
   
   vector<vector<double>> Cspeed, FluxAdvCC_x, FluxDifCC_x, FluxAdv_x, FluxDif_x;
   Cspeed.assign(nXcc,vector<double>(nZcc));
   FluxAdvCC_x.assign(nXcc,vector<double>(nZcc));
   FluxDifCC_x.assign(nXcc,vector<double>(nZcc));
   FluxAdv_x.assign(nXce,vector<double>(nZcc));
   FluxDif_x.assign(nXce,vector<double>(nZcc));
   //
   vector<vector<double>> FluxAdvCC_z, FluxDifCC_z, FluxAdv_z, FluxDif_z;
   Cspeed.assign(nXcc,vector<double>(nZcc));
   FluxAdvCC_z.assign(nXcc,vector<double>(nZcc));
   FluxDifCC_z.assign(nXcc,vector<double>(nZcc));
   FluxAdv_z.assign(nXcc,vector<double>(nZce));
   FluxDif_z.assign(nXcc,vector<double>(nZce));
   

   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeed = abs(F0); // adv flux jacobian
   FluxAdvCC_x = double(XfluxDir)*F0*F0*0.5;
   FluxAdvCC_z = double(ZfluxDir)*F0*F0*0.5;


   // compute diffusive flux using 
   // standard centered scheme
   //
   Xgrid.DDX(FluxDifCC_x,F0);
   FluxDifCC_x = -K*FluxDifCC_x;
   Xgrid.communicate(FluxDifCC_x);

   Xgrid.DDX(FluxDif_x,F0);
   FluxDif_x = -K*FluxDif_x;
   Xgrid.communicate(FluxDif_x);
   //
   Xgrid.DDZ(FluxDifCC_z,F0);
   FluxDifCC_z = -K*FluxDifCC_z;
   Xgrid.communicate(FluxDifCC_z);

   Xgrid.DDZ(FluxDif_z,F0);
   FluxDif_z = -K*FluxDif_z;
   Xgrid.communicate(FluxDif_z);


   // compute advective flux using
   // specified scheme from input file
   //
   if(advScheme0 == "TVD") {
      //Xgrid.computeFluxTVD(FluxAdv,FluxL,FluxR,FluxLimL,FluxLimR,
      //                     FluxAdvCC,Cspeed,F0,0,order);
      Xgrid.computeFluxTVD(Flux_x,FluxL_x,FluxR_x,FluxLimL_x,FluxLimR_x,
                           FluxAdvCC_x+FluxDifCC_x,Cspeed,F0,0,order);
      Xgrid.computeFluxTVD(Flux_z,FluxL_z,FluxR_z,FluxLimL_z,FluxLimR_z,
                           FluxAdvCC_z+FluxDifCC_z,Cspeed,F0,1,order);
   }      
   else {
      Xgrid.InterpToCellEdges(FluxAdv_x,FluxAdvCC_x,Cspeed,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxAdv_z,FluxAdvCC_z,Cspeed,advScheme0,1);
      Flux_x = FluxAdv_x + FluxDif_x;
      Flux_z = FluxAdv_z + FluxDif_z;
   } 
   Xgrid.communicate(Flux_x);
   Xgrid.communicate(FluxL_x);
   Xgrid.communicate(FluxR_x);
   //
   Xgrid.communicate(Flux_z);
   Xgrid.communicate(FluxL_z);
   Xgrid.communicate(FluxR_z);


} // end computeFluxes



void setXminBoundary(vector<vector<double>>& var,
		     const double C0, const double C1)
{
   //const int thisnX = var.size();
   const int thisnZ = var[0].size();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift=mesh->nXg;
   //const int nZg    = mesh->nZg;

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
   const int ishift = thisnX-mesh->nXg;
   //const int nZg    = mesh->nZg;

   for (auto i=ishift; i<thisnX; i++) {
      for (auto j=0; j<thisnZ; j++) {
      //for (auto j=nZg; j<thisnZ-nZg; j++) {
         var[i][j] = C0 + C1*var[2*ishift-i-1][j];
      }
   }   

}

void setZminBoundary(vector<vector<double>>& var,
		     const double C0, const double C1)
{
   const int thisnX = var.size();
   //const int thisnZ = var[0].size();

   domainGrid* mesh = domainGrid::mesh;
   const int jshift=mesh->nZg;
   //const int nZg    = mesh->nZg;

   for (auto i=0; i<thisnX; i++) {
      for (auto j=0; j<jshift; j++) {
         var[i][jshift-j-1] = C0 + C1*var[i][jshift+j];
      }
   }

}


void setZmaxBoundary(vector<vector<double>>& var,
		     const double C0, const double C1)
{
   const int thisnX = var.size();
   const int thisnZ = var[0].size();

   domainGrid* mesh = domainGrid::mesh;
   const int jshift = thisnZ-mesh->nZg;
   //const int nZg    = mesh->nZg;

   for (auto i=0; i<thisnX; i++) {
      for (auto j=jshift; j<thisnZ; j++) {
         var[i][j] = C0 + C1*var[i][2*jshift-j-1];
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
   const double dZ = Xgrid.dX;
   double dtmaxDif = 0.5/K*(dX*dX+dZ*dZ)/(dX*dX+dZ*dZ);
   double dtmaxAdv = 0.5*dX*dZ/(dX+dZ)/Umax;
   double dtmax = min(dtmaxDif,dtmaxAdv);
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);
   if(procID==0) {
      cout << "max stable time step is " << dtmax << endl;
      cout << endl; 
   }
}

