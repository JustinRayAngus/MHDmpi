/***
 * 
 * physics module for 2D burgers equation:
 * df/dt + d(cx*f^2/2)/dx + d(cz*f^2/2)/dz = K*(d2f/dx2 + d2f/dz2)
 *
 * Note: I original computed advection and diff fluxes
 *       separately when using TVD scheme, which I only 
 *       used for adv flux. However, I found that sometimes
 *       I needed a smaller time step than expected or else
 *       oscillations would occur.
 *
 *       When doing TVD, I now use the total flux and it works
 *       just fine. Though this is not proper use...
 *
 *       Not sure if above still holds. I think I had time-step calculated wrong...
 *       Ok. Pretty sure I'm just not doing time-step correct when using upwinding
 *       Rather than C2 for advection... Need to clean this up
 *       See Section 11.8 in Pozrikidis
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
#include <cmath>

#include "json/json.h"
#include "vectorMath.h"
#include "matrix2D.h"
#include "domainGrid.h"
#include "timeDomain.h"
#include "HDF5dataFile.h"

#include "Physics.h"


using namespace std;

string advScheme0, TVDlimiter;
bool useLaxSplitting;
double K, cx, cz;
int order;         // order in time (1 or 2)
matrix2D<double> F0, F0old;    // function
matrix2D<double> FluxLimL_x, FluxLimR_x;
matrix2D<double> Flux_x, FluxR_x, FluxL_x;  // flux at cell-edges   
matrix2D<double> FluxLimL_z, FluxLimR_z;
matrix2D<double> Flux_z, FluxR_z, FluxL_z;  // flux at cell-edges   


void computeFluxes(const domainGrid&, const int);
void setXminBoundary(matrix2D<double>&, 
		     const double, const double);
void setXmaxBoundary(matrix2D<double>&, 
		     const double, const double);
void setZminBoundary(matrix2D<double>&, 
		     const double, const double);
void setZmaxBoundary(matrix2D<double>&, 
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
   
   F0.initialize(nXcc,nZcc,0.0);
   F0old.initialize(nXcc,nZcc,0.0);
   //
   FluxLimL_x.initialize(nXce,nZcc,0.0);
   FluxLimR_x.initialize(nXce,nZcc,0.0);
   Flux_x.initialize(nXce,nZcc,0.0);
   FluxR_x.initialize(nXce,nZcc,0.0);
   FluxL_x.initialize(nXce,nZcc,0.0);
   //
   FluxLimL_z.initialize(nXcc,nZce,0.0);
   FluxLimR_z.initialize(nXcc,nZce,0.0);
   Flux_z.initialize(nXcc,nZce,0.0);
   FluxR_z.initialize(nXcc,nZce,0.0);
   FluxL_z.initialize(nXcc,nZce,0.0);

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
         advScheme0=="QUICK" || advScheme0=="TVD" || "WENO5") {
         if(procID==0) {
            cout << "advection diff/interp scheme is " << advScheme0 << endl;
         }
      }
      else {
         cout << "advection scheme " << advScheme0 << " is not valid " << endl;
         cout << "valid types are C2, U1, QUICK, TVD, and WENO5 " << endl;
         exit (EXIT_FAILURE);
      }
      TVDlimiter = "superbee";
      if(advScheme0=="TVD") {
         Json::Value TVDlimiterVal = Phys.get("TVDlimiter",defValue);
         if(TVDlimiterVal != defValue) {
	    TVDlimiter = TVDlimiterVal.asString();
	 }
         if(procID==0) cout << "TVD limiter = " << TVDlimiter << endl;
      }
      useLaxSplitting = false;
      Json::Value useLaxSplittingVal = Phys.get("useLaxSplitting",defValue);
      if(useLaxSplittingVal != defValue) {
         useLaxSplitting = useLaxSplittingVal.asBool();
         if(procID==0 && useLaxSplitting) cout << "using Lax Friedrich's flux splitting " << endl;
      }

      //   get diffusion and advection coefficents from input file
      //
      Json::Value KVal = Phys.get("diffC",defValue);
      if(KVal != defValue) K = KVal.asDouble();
      else K = 0.0;
      Json::Value cxVal = Phys.get("advCx",defValue);
      if(cxVal != defValue) cx = cxVal.asDouble();
      else cx = 1.0;
      Json::Value czVal = Phys.get("advCz",defValue);
      if(czVal != defValue) cz = czVal.asDouble();
      else cz = 1.0;

      if(procID==0) {
         cout << "diffusion coefficent = " << K << endl;
         cout << "x-advection coefficent = " << cx << endl;
         cout << "z-advection coefficent = " << cz << endl;
      }
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
      
   }
   else {
      cout << "value for key \"Physics\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
  

   Xgrid.setInitialProfile(F0,Phys);
   if(procID==0) setXminBoundary(F0, 0.0, 1.0);   
   if(procID==numProcs-1) setXmaxBoundary(F0, 0.0, 1.0); 
   setZminBoundary(F0, 0.0, 1.0);  
   setZmaxBoundary(F0, 0.0, 1.0);  
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
   const int nXcc = Xgrid.nXcc;
   const int nXg = Xgrid.nXg;
   const int nZcc = Xgrid.nZcc;
   const int nZg = Xgrid.nZg;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   // Explicit forward advance from n to n+1/2 using Flux(n)
   // (calc of Flux(t=0) is done during initilization)
   //
   for (auto i=nXg; i<nXcc-nXg; i++) {
      for (auto j=nZg; j<nZcc-nZg; j++) {
         F0(i,j) = F0old(i,j) - dt/2.0*(Flux_x(i,j) - Flux_x(i-1,j))/Xgrid.dX
                              - dt/2.0*(Flux_z(i,j) - Flux_z(i,j-1))/Xgrid.dZ;
      }
   }

   // apply boundary conditions, communicate, compute Flux
   //
   if(procID==0) setXminBoundary(F0, 0.0, 1.0);   
   if(procID==numProcs-1) setXmaxBoundary(F0, 0.0, 1.0);   
   setZminBoundary(F0, 0.0, 1.0);
   setZmaxBoundary(F0, 0.0, 1.0);
   Xgrid.communicate(F0);
   computeFluxes(Xgrid,order); // compute RHS using F0(n+1/2)


   // Explicit forward advance from n to n+1 using Flux(n+1/2)
   //
   for (auto i=nXg; i<nXcc-nXg; i++) {
      for (auto j=nZg; j<nZcc-nZg; j++) {
         F0(i,j) = F0old(i,j) - dt*(Flux_x(i,j) - Flux_x(i-1,j))/Xgrid.dX
                              - dt*(Flux_z(i,j) - Flux_z(i,j-1))/Xgrid.dZ;
      }
   }
   

   // apply boundary conditions, communicate, compute Flux
   //
   if(procID==0) setXminBoundary(F0, 0.0, 1.0);   
   if(procID==numProcs-1) setXmaxBoundary(F0, 0.0, 1.0);   
   setZminBoundary(F0, 0.0, 1.0);
   setZmaxBoundary(F0, 0.0, 1.0);
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


   const int nXce = Flux_x.size0();
   const int nZce = Flux_z.size1();
   const int nXcc = F0.size0();
   const int nZcc = F0.size1();
   
   matrix2D<double> Cspeed, FluxAdvCC_x, FluxDifCC_x, FluxAdv_x, FluxDif_x;
   Cspeed.initialize(nXcc,nZcc,0.0);
   FluxAdvCC_x.initialize(nXcc,nZcc,0.0);
   FluxDifCC_x.initialize(nXcc,nZcc,0.0);
   FluxAdv_x.initialize(nXce,nZcc,0.0);
   FluxDif_x.initialize(nXce,nZcc,0.0);
   //
   matrix2D<double> FluxAdvCC_z, FluxDifCC_z, FluxAdv_z, FluxDif_z;
   FluxAdvCC_z.initialize(nXcc,nZcc,0.0);
   FluxDifCC_z.initialize(nXcc,nZcc,0.0);
   FluxAdv_z.initialize(nXcc,nZce,0.0);
   FluxDif_z.initialize(nXcc,nZce,0.0);
   

   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   FluxAdvCC_x = cx*F0*F0*0.5;
   FluxAdvCC_z = cz*F0*F0*0.5;
   //Cspeed_x = abs(cx*F0); // adv flux jacobian
   //Cspeed_z = abs(cz*F0); // adv flux jacobian
   Cspeed = sqrt(cx*cx+cz*cz)*abs(F0); // adv flux jacobian
   //Cspeed = abs(F0); // adv flux jacobian

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
      Xgrid.computeFluxTVD(FluxAdv_x,FluxL_x,FluxR_x,FluxLimL_x,FluxLimR_x,
                           FluxAdvCC_x,Cspeed,F0,TVDlimiter,0,order);
      Xgrid.computeFluxTVD(FluxAdv_z,FluxL_z,FluxR_z,FluxLimL_z,FluxLimR_z,
                           FluxAdvCC_z,Cspeed,F0,TVDlimiter,1,order);
   }      
   else {
      if(useLaxSplitting) {
         Xgrid.computeFluxTVDsimple(FluxAdv_x,FluxL_x,FluxR_x,
                                    FluxAdvCC_x,Cspeed,F0,advScheme0,0);
         Xgrid.computeFluxTVDsimple(FluxAdv_z,FluxL_z,FluxR_z,
                                    FluxAdvCC_z,Cspeed,F0,advScheme0,1);
      } else {
         Xgrid.InterpToCellEdges(FluxAdv_x,FluxAdvCC_x,Cspeed,advScheme0,0);
         Xgrid.InterpToCellEdges(FluxAdv_z,FluxAdvCC_z,Cspeed,advScheme0,1);
      }
   } 
   Flux_x = FluxAdv_x + FluxDif_x;
   Flux_z = FluxAdv_z + FluxDif_z;
   //
   Xgrid.communicate(Flux_x);
   Xgrid.communicate(FluxL_x);
   Xgrid.communicate(FluxR_x);
   //
   Xgrid.communicate(Flux_z);
   Xgrid.communicate(FluxL_z);
   Xgrid.communicate(FluxR_z);


} // end computeFluxes



void setXminBoundary(matrix2D<double>& var,
		     const double C0, const double C1)
{
   //const int thisnX = var.size0();
   const int thisnZ = var.size1();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift=mesh->nXg;
   //const int nZg    = mesh->nZg;

   for (auto i=0; i<ishift; i++) {
      for (auto j=0; j<thisnZ; j++) {
         var(ishift-i-1,j) = C0 + C1*var(ishift+i,j);
      }
   }

}

void setXmaxBoundary(matrix2D<double>& var,
		     const double C0, const double C1)
{
   const int thisnX = var.size0();
   const int thisnZ = var.size1();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift = thisnX-mesh->nXg;
   //const int nZg    = mesh->nZg;

   for (auto i=ishift; i<thisnX; i++) {
      for (auto j=0; j<thisnZ; j++) {
      //for (auto j=nZg; j<thisnZ-nZg; j++) {
         var(i,j) = C0 + C1*var(2*ishift-i-1,j);
      }
   }   

}

void setZminBoundary(matrix2D<double>& var,
		     const double C0, const double C1)
{
   const int thisnX = var.size0();
   //const int thisnZ = var.size1();

   domainGrid* mesh = domainGrid::mesh;
   const int jshift=mesh->nZg;
   //const int nZg    = mesh->nZg;

   for (auto i=0; i<thisnX; i++) {
      for (auto j=0; j<jshift; j++) {
         var(i,jshift-j-1) = C0 + C1*var(i,jshift+j);
      }
   }

}


void setZmaxBoundary(matrix2D<double>& var,
		     const double C0, const double C1)
{
   const int thisnX = var.size0();
   const int thisnZ = var.size1();

   domainGrid* mesh = domainGrid::mesh;
   const int jshift = thisnZ-mesh->nZg;
   //const int nZg    = mesh->nZg;

   for (auto i=0; i<thisnX; i++) {
      for (auto j=jshift; j<thisnZ; j++) {
         var(i,j) = C0 + C1*var(i,2*jshift-j-1);
      }
   }   

}


void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid, const int verbose)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   double dtmax, Factor, Umax_x, Umax_z;
   Umax_x = max(abs(cx*F0));
   Umax_z = max(abs(cz*F0));
   //Umax = max(sqrt(cx*cx + cz*cz)*abs(F0));
   //cout << "Umax = " << Umax << endl;

   const double dX = Xgrid.dX;
   const double dZ = Xgrid.dX;
   dtmax = 0.5*dX*dZ/(Umax_x*dZ+Umax_z*dX);

   if(advScheme0=="C2") {
      if(K>0.0) {
         double dtmaxDif = 0.5/K*dX*dX*dZ*dZ/(dX*dX + dZ*dZ);
         double dtmaxAD = 2.0*K/(Umax_x*Umax_x + Umax_z+Umax_z);
         dtmax = min(dtmax,dtmaxAD);
         dtmax = min(dtmax,dtmaxDif);
      }
   } else { // This is for U1 adv + C2 diff
      Factor = 2.0*K/dX/dX + 2.0*K/dZ/dZ + Umax_x/dX + Umax_z/dZ;
      dtmax = min(dtmax,1.0/Factor);
   }
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);


   if(procID==0 && verbose) {
      //if(dtmax==dtmaxDif) {
      //   cout << "max stable time step is set by diffusion" << endl;
      //} else {
      //   cout << "max stable time step is set by advection" << endl;
      //}
      cout << "dtSim = " << dtSim << endl;
      cout << endl; 
   }
}


