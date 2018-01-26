/***
 * 
 * 2D physics module for testing shock capturing schemes
 *
 * see G.A. Sod, J. of Comp. Phys. 27, 1-31 (1978)
 *
 * upwding (U1) scheme seems to fail at really high
 * grid resolution (fine at nx=3000, fails at nx=4000)
 * I think this is result of how pressure term in 
 * momentum equation is handled
 *
 * QUICK doesn't work well at all
 *
 * TVD works very well, and no problem at high resolution
 * (nx=6000 is fine)
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
#include "matrix2D.h"
#include "domainGrid.h"
#include "timeDomain.h"
#include "HDF5dataFile.h"

#include "Physics.h"


using namespace std;

string advScheme0;    // advection differencing scheme
double gamma0;        // adiabatic coefficient
int Nsub;             // time-solver subcycle steps
vector<vector<double>> N, Mx, Mz, E;   // function
vector<vector<double>> F0, Cs, Vx, Vz, P;      // function
vector<vector<double>> F0old, Nold, Mxold, Mzold, Eold;
vector<vector<double>> FluxLimLx, FluxLimRx, FluxRx, FluxLx;  // flux at cell-edges   
vector<vector<double>> FluxLimLz, FluxLimRz, FluxRz, FluxLz;  // flux at cell-edges   
vector<vector<double>> FluxN_x, FluxMx_x, FluxMz_x, FluxE_x;
vector<vector<double>> FluxN_z, FluxMx_z, FluxMz_z, FluxE_z;


matrix2D<double> matA(2,2,1);
matrix2D<double> matB;
matrix2D<double> matC;


void computeFluxes(const domainGrid&, const int);
void setXminBoundary(vector<vector<double>>&, const double, const double);
void setXmaxBoundary(vector<vector<double>>&, const double, const double);
void setZminBoundary(vector<vector<double>>&, const double, const double);
void setZmaxBoundary(vector<vector<double>>&, const double, const double);

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
   //cout << "mesh->nXcc = " << (*mesh).nXcc << endl;
   //cout << "mesh->nXcc = " << mesh->nXcc << endl;

   const int nXcc = Xgrid.Xcc.size();
   const int nXce = Xgrid.Xce.size();
   const int nZcc = Xgrid.Zcc.size();
   const int nZce = Xgrid.Zce.size();
   

   //matA(nXcc,nZcc,1.0);
   //matB(nXcc,nZcc,1.0);
   matA.initialize(nXcc,nZcc,2.0);
   matB.initialize(nXcc,nZcc,2.0);
   matC.initialize(nXcc,nZcc,1.0);
   matC = matC/matB;
   matC = matC*10.0;
   matC *= matC;
   matC = 10.0*matC;
   matC *= 2.0;
   matC += matC;
   matC -= 0.5*matC;
   matC /= matB;
   matC += 1.0;
   matC = matC+1.0;
   matC = 1.0+matC;
   matC /= matB;
   matC = matC/10.0;
   matC = matC-1.0;
   matC = -10.0/matC;
   matC = abs(matC);
   matC = pow(matC,3.1);
   matC = cos(matC);
   matC = sin(matC);
   matC = exp(matC);
   matC = sqrt(matC);
   matC = log(matC);
   matC = tanh(matC);
   matC(1,1) = -0.1;

   cout << "matB.size0() = " << matB.size0() << endl;
   cout << "matC.size0() = " << matC.size0() << endl;
   //cout << "matA.nX      = " << matA.nX << endl;
   cout << "matB.size1() = " << matB.size1() << endl;
   cout << "matC.size1() = " << matC.size1() << endl;
   //cout << "matA.nZ      = " << matA.nZ << endl;
   cout << "matB(2,2)    = " << matB(2,2) << endl;
   cout << "matC(2,2)    = " << matC(2,2) << endl;
   cout << "min(matC)    = " << min(matC) << endl;
   cout << "max(matC)    = " << max(matC) << endl;


   F0.assign(nXcc,vector<double>(nZcc,0.0));
   F0old.assign(nXcc,vector<double>(nZcc,0.0));
   N.assign(nXcc,vector<double>(nZcc,0.0));
   Nold.assign(nXcc,vector<double>(nZcc,0.0));
   Mx.assign(nXcc,vector<double>(nZcc,0.0));
   Mxold.assign(nXcc,vector<double>(nZcc,0.0));
   Mz.assign(nXcc,vector<double>(nZcc,0.0));
   Mzold.assign(nXcc,vector<double>(nZcc,0.0));
   E.assign(nXcc,vector<double>(nZcc,0.0));
   Eold.assign(nXcc,vector<double>(nZcc,0.0));
   P.assign(nXcc,vector<double>(nZcc,0.0));
   Cs.assign(nXcc,vector<double>(nZcc,0.0));
   Vx.assign(nXcc,vector<double>(nZcc,0.0));
   Vz.assign(nXcc,vector<double>(nZcc,0.0));
   //
   FluxRx.assign(nXce,vector<double>(nZcc,0.0));
   FluxLx.assign(nXce,vector<double>(nZcc,0.0));
   FluxLimLx.assign(nXce,vector<double>(nZcc,0.0));
   FluxLimRx.assign(nXce,vector<double>(nZcc,0.0));
   //
   FluxRz.assign(nXcc,vector<double>(nZce,0.0));
   FluxLz.assign(nXcc,vector<double>(nZce,0.0));
   FluxLimLz.assign(nXcc,vector<double>(nZce,0.0));
   FluxLimRz.assign(nXcc,vector<double>(nZce,0.0));
   //
   FluxN_x.assign(nXce,vector<double>(nZcc,0.0));
   FluxMx_x.assign(nXce,vector<double>(nZcc,0.0));
   FluxMz_x.assign(nXce,vector<double>(nZcc,0.0));
   FluxE_x.assign(nXce,vector<double>(nZcc,0.0));
   //
   FluxN_z.assign(nXcc,vector<double>(nZce,0.0));
   FluxMx_z.assign(nXcc,vector<double>(nZce,0.0));
   FluxMz_z.assign(nXcc,vector<double>(nZce,0.0));
   FluxE_z.assign(nXcc,vector<double>(nZce,0.0));
   //
   //
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      Json::Value advScheme = Phys.get("advScheme",defValue);
      Json::Value gammaVal  = Phys.get("gammaC",defValue);
      Json::Value NsubVal   = Phys.get("Nsub",defValue);
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
      
      Nsub = NsubVal.asInt();
      if(procID==0) cout << "Nsub = " << Nsub << endl;
      if(gamma0 < 1.0) {
         printf("ERROR: Nsub must be int >= 1\n");
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
      Xgrid.setInitialProfile(matA,Nvar);
      if(procID==0) setXminBoundary(N, 0.0, 1.0);   
      //if(procID==numProcs-1) setXmaxBoundary(N, 0.125);   
      if(procID==numProcs-1) setXmaxBoundary(N, 0.0, 1.0);   
      setZminBoundary(N,0.0,1.0);
      setZmaxBoundary(N,0.0,1.0);
      Xgrid.communicate(N);
      Nold  = N;
   } else {
      cout << "value for Physics variable \"N\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
 
   const Json::Value Vxvar = Phys.get("Vx",defValue);
   if(Vxvar.isObject()) { 
      Xgrid.setInitialProfile(Vx,Vxvar);
      if(procID==0) setXminBoundary(Vx, 0.0, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(Vx, 0.0, 1.0);   
      Xgrid.communicate(Vx);
      Mx = N*Vx;
      Mxold  = Mx;
   } else {
      cout << "value for Physics variable \"Vx\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   
   const Json::Value Vzvar = Phys.get("Vz",defValue);
   if(Vzvar.isObject()) { 
      Xgrid.setInitialProfile(Vz,Vzvar);
      if(procID==0) setXminBoundary(Vz, 0.0, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(Vz, 0.0, 1.0);   
      Xgrid.communicate(Vz);
      Mz = N*Vz;
      Mzold  = Mz;
   } else {
      cout << "value for Physics variable \"Vz\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   const Json::Value Pvar = Phys.get("P",defValue);
   if(Pvar.isObject()) { 
      Xgrid.setInitialProfile(P,Pvar);
      if(procID==0) setXminBoundary(P, 0.0, 1.0);   
      //if(procID==numProcs-1) setXmaxBoundary(P, 0.1);   
      if(procID==numProcs-1) setXmaxBoundary(P, 0.0, 1.0);   
      Xgrid.communicate(P);
      setZminBoundary(P,0.0,1.0);
      setZmaxBoundary(P,0.0,1.0);
      E = P/(gamma0-1.0) + 0.5*(Mx*Vx+Mz*Vz);
      Eold  = E;
   } else {
      cout << "value for Physics variable \"P\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   Cs = sqrt(gamma0*P/N);

   // need inital flux calculation before add to output
   // and before first call to Physics.advance   
   computeFluxes(Xgrid,1);   
  
   cout << "matA(2,10) = " << matA(2,10) << endl;
   cout << "matA(2,150) = " << matA(2,150) << endl;

  
   dataFile.add(matA, "matA", 1);  // density 
   dataFile.add(N, "N", 1);  // density 
   dataFile.add(Mx, "Mx", 1);  // momentum density 
   dataFile.add(Mz, "Mz", 1);  // momentum density 
   dataFile.add(E, "E", 1);  // energy density
   dataFile.add(P, "P", 1);  // pressure
   dataFile.add(Vx, "Vx", 1);  // velocity
   dataFile.add(Vz, "Vz", 1);  // velocity
   dataFile.add(Cs,"Cs",1);  // sound speed
   dataFile.add(gamma0,"gamma0",0); 
   //
   dataFile.add(FluxN_x,  "FluxN_x",  1);  
   dataFile.add(FluxMx_x, "FluxMx_x", 1);  
   dataFile.add(FluxMz_x, "FluxMz_x", 1);  
   dataFile.add(FluxE_x,  "FluxE_x",  1);  
   //
   dataFile.add(FluxN_z,  "FluxN_z",  1);  
   dataFile.add(FluxMx_z, "FluxMx_z", 1);  
   dataFile.add(FluxMz_z, "FluxMz_z", 1);  
   dataFile.add(FluxE_z,  "FluxE_z",  1);  
   //
   //dataFile.add(FluxLimLx, "FluxLimLx", 1);  
   //dataFile.add(FluxLimRx, "FluxLimRx", 1);  
   //dataFile.add(FluxRx, "FluxRx", 1);
   //dataFile.add(FluxLx, "FluxLx", 1);

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

   //domainGrid* mesh = domainGrid::mesh;


   // Explicit forward advance from n to n+1 using subcycling in time 
   //  
   double thisdt;
   for (auto n=1; n<Nsub+1; n++) {
      thisdt = dt*n/Nsub;

      for (auto i=nXg; i<nXcc-nXg; i++) {
         for (auto j=nZg; j<nZcc-nZg; j++) {
            N[i][j] = Nold[i][j] - thisdt*(FluxN_x[i][j] - FluxN_x[i-1][j])/Xgrid.dX
                                 - thisdt*(FluxN_z[i][j] - FluxN_z[i][j-1])/Xgrid.dZ;
            Mx[i][j]= Mxold[i][j]- thisdt*(FluxMx_x[i][j]-FluxMx_x[i-1][j])/Xgrid.dX
                                 - thisdt*(FluxMx_z[i][j]-FluxMx_z[i][j-1])/Xgrid.dZ;
            Mz[i][j]= Mzold[i][j]- thisdt*(FluxMz_x[i][j]-FluxMz_x[i-1][j])/Xgrid.dX
                                 - thisdt*(FluxMz_z[i][j]-FluxMz_z[i][j-1])/Xgrid.dZ;
            E[i][j] = Eold[i][j] - thisdt*(FluxE_x[i][j] - FluxE_x[i-1][j])/Xgrid.dX
                                 - thisdt*(FluxE_z[i][j] - FluxE_z[i][j-1])/Xgrid.dZ;
         }
      }

      if(procID==0) {
         setXminBoundary(N, 0.0, 1.0);   
         setXminBoundary(Mx, 0.0, 1.0);   
         setXminBoundary(Mz, 0.0, 1.0);   
         //setXminBoundary(E, 1.0/(gamma0-1.0));   
         setXminBoundary(E, 0.0, 1.0);   
      }
      if(procID==numProcs-1) {
         //setXmaxBoundary(N, 0.125);   
         //setXmaxBoundary(Mx, 0.0);   
         //setXmaxBoundary(E, 0.1/(gamma0-1.0));   
         setXmaxBoundary(N, 0.0, 1.0);   
         setXmaxBoundary(Mz, 0.0, 1.0);   
         setXmaxBoundary(Mz, 0.0, 1.0);   
         setXmaxBoundary(E, 0.0, 1.0);   
      }
      setZmaxBoundary(N, 0.0, 1.0);   
      setZmaxBoundary(Mx, 0.0, 1.0);   
      setZmaxBoundary(Mz, 0.0, 1.0);   
      setZmaxBoundary(E, 0.0, 1.0);   
      setZminBoundary(N, 0.0, 1.0);   
      setZminBoundary(Mx, 0.0, 1.0);   
      setZminBoundary(Mz, 0.0, 1.0);   
      setZminBoundary(E, 0.0, 1.0);   
      Xgrid.communicate(N);
      Xgrid.communicate(Mx);
      Xgrid.communicate(Mz);
      Xgrid.communicate(E);
      
      if(n==1 && Nsub==2) computeFluxes(Xgrid,2);
      else computeFluxes(Xgrid,1);
      //computeFluxes(Xgrid,2);

   }

   Nold = N;
   Mxold = Mx;
   Mzold = Mz;
   Eold = E;

} // end Physics.advance


void computeFluxes(const domainGrid& Xgrid, const int order)
{

   //const int nXce = FluxN_x.size();
   //const int nZce = FluxN_z[0].size();
   const int nXcc = N.size();
   const int nZcc = N[0].size();


   vector<vector<double>> Cspeed, Cspeed2;
   vector<vector<double>> FluxNcc_x, FluxMxcc_x, FluxMzcc_x, FluxEcc_x;
   vector<vector<double>> FluxNcc_z, FluxMxcc_z, FluxMzcc_z, FluxEcc_z;
   
   FluxNcc_x.assign(nXcc,vector<double>(nZcc,0.0));
   FluxMxcc_x.assign(nXcc,vector<double>(nZcc,0.0));
   FluxMzcc_x.assign(nXcc,vector<double>(nZcc,0.0));
   FluxEcc_x.assign(nXcc,vector<double>(nZcc,0.0));
   //
   FluxNcc_z.assign(nXcc,vector<double>(nZcc,0.0));
   FluxMxcc_z.assign(nXcc,vector<double>(nZcc,0.0));
   FluxMzcc_z.assign(nXcc,vector<double>(nZcc,0.0));
   FluxEcc_z.assign(nXcc,vector<double>(nZcc,0.0));
   //
   Cspeed.assign(nXcc,vector<double>(nZcc,0.0));
   Cspeed2.assign(nXcc,vector<double>(nZcc,0.0));

   //  define derived variables
   //
   Vx  = Mx/N;
   Vz  = Mz/N;
   P  = (E-0.5*(Mx*Mx+Mz*Mz)/N)*(gamma0-1);
   Cs = sqrt(gamma0*P/N);

   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   
   Cspeed = sqrt(Vx*Vx+Vz*Vz)+Cs; // adv flux jacobian
   /*
   Cspeed  = abs(Vx + Cs); // adv flux jacobian
   Cspeed2 = abs(Vx - Cs); // adv flux jacobian
   for (auto i=0; i<nXcc; i++) {
      for (auto j=0; j<nZcc; j++) {
         Cspeed[i][j] = max(Cspeed[i][j],Cspeed2[i][j]);
      }
   }
   */

   FluxNcc_x = Mx;
   FluxMxcc_x = Mx*Mx/N + P;
   FluxMzcc_x = Mz*Mx/N;
   FluxEcc_x = Vx*(E + P);
   //
   FluxNcc_z = Mz;
   FluxMxcc_z = Mx*Mz/N;
   FluxMzcc_z = Mz*Mz/N + P;
   FluxEcc_z = Vz*(E + P);
   

   // compute advective flux using
   // specified scheme from input file
   //
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN_x, FluxLx,FluxRx,FluxLimLx,FluxLimRx,
                           FluxNcc_x,Cspeed,N,0,order);
      Xgrid.computeFluxTVD(FluxMx_x,FluxLx,FluxRx,FluxLimLx,FluxLimRx,
                           FluxMxcc_x,Cspeed,Mx,0,order);
      Xgrid.computeFluxTVD(FluxMz_x,FluxLx,FluxRx,FluxLimLx,FluxLimRx,
                           FluxMzcc_x,Cspeed,Mz,0,order);
      Xgrid.computeFluxTVD(FluxE_x, FluxLx,FluxRx,FluxLimLx,FluxLimRx,
                           FluxEcc_x,Cspeed,E,0,order);
      //
      Xgrid.computeFluxTVD(FluxN_z, FluxLz,FluxRz,FluxLimLz,FluxLimRz,
                           FluxNcc_z,Cspeed,N,1,order);
      Xgrid.computeFluxTVD(FluxMx_z,FluxLz,FluxRz,FluxLimLz,FluxLimRz,
                           FluxMxcc_z,Cspeed,Mx,1,order);
      Xgrid.computeFluxTVD(FluxMz_z,FluxLz,FluxRz,FluxLimLz,FluxLimRz,
                           FluxMzcc_z,Cspeed,Mz,1,order);
      Xgrid.computeFluxTVD(FluxE_z, FluxLz,FluxRz,FluxLimLz,FluxLimRz,
                           FluxEcc_z,Cspeed,E,1,order);
   }      
   else {
      Xgrid.InterpToCellEdges(FluxN_x, FluxNcc_x,N,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxMx_x,FluxMxcc_x,Mx,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxMz_x,FluxMzcc_x,Mz,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxE_x, FluxEcc_x,E,advScheme0,0);
      //
      Xgrid.InterpToCellEdges(FluxN_z, FluxNcc_z,N,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxMx_z,FluxMxcc_z,Mx,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxMz_z,FluxMzcc_z,Mz,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxE_z, FluxEcc_z,E,advScheme0,0);
   } 


} // end computeFluxes


void setXminBoundary(vector<vector<double>>& var,
		     const double C0, const double C1)
{
   const int thisnZ = var[0].size();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift = mesh->nXg;

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

   for (auto i=ishift; i<thisnX; i++) {  
      for (auto j=0; j<thisnZ; j++) {
	 var[i][j] = C0 + C1*var[2*ishift-i-1][j];
      }
   }

}

void setZminBoundary(vector<vector<double>>& var,
		     const double C0, const double C1)
{
   const int thisnX = var.size();

   domainGrid* mesh = domainGrid::mesh;
   const int jshift = mesh->nZg;

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
   
   vector<vector<double>> Cchar;
   double Cmax;
   //Cchar = abs(Vx)+Cs;
   Cchar = sqrt(Vx*Vx+Vz*Vz)+Cs;
   Cmax = max(Cchar);
   //cout << "Cmax = " << Cmax << endl;

   const double dX = Xgrid.dX;
   const double dZ = Xgrid.dZ;
   
   double dtmax = 0.5*dX*dZ/(dX+dZ)/Cmax;
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);
   if(procID==0) {
      cout << "max stable time step is " << dtmax << endl;
   }
}

