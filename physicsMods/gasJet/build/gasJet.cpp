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
string geometry0;
double gamma0;        // adiabatic coefficient
int Nsub;             // time-solver subcycle steps
matrix2D<double> N, Mx, Mz, E;   // function
matrix2D<double> F0, Cs, Vx, Vz, P;      // function
matrix2D<double> F0old, Nold, Mxold, Mzold, Eold;
matrix2D<double> FluxLimLx, FluxLimRx, FluxRx, FluxLx;  // flux at cell-edges   
matrix2D<double> FluxLimLz, FluxLimRz, FluxRz, FluxLz;  // flux at cell-edges   
matrix2D<double> FluxN_x, FluxMx_x, FluxMz_x, FluxE_x;
matrix2D<double> FluxN_z, FluxMx_z, FluxMz_z, FluxE_z;

matrix2D<double> hy_cc, hy_ce, hy_z, hy_xz; // y-dimension lame coefficient

double Rinlet, Pinlet, Minlet, Tplenum;
double Pplenum, Nplenum;
double Ninlet, Tinlet, Einlet, Vzinlet, Mzinlet;

double Tthresh = 1.0e-2, Nthresh = 1.0e-4;
double Ethresh, Pthresh = Nthresh*Tthresh;

void computeFluxes(const domainGrid&, const int);

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

   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXce = Xgrid.nXce2;
   const int nZce = Xgrid.nZce2;
   

   F0.initialize(nXcc,nZcc,0.0);
   F0old.initialize(nXcc,nZcc,0.0);
   N.initialize(nXcc,nZcc,0.0);
   Nold.initialize(nXcc,nZcc,0.0);
   Mx.initialize(nXcc,nZcc,0.0);
   Mxold.initialize(nXcc,nZcc,0.0);
   Mz.initialize(nXcc,nZcc,0.0);
   Mzold.initialize(nXcc,nZcc,0.0);
   E.initialize(nXcc,nZcc,0.0);
   Eold.initialize(nXcc,nZcc,0.0);
   P.initialize(nXcc,nZcc,0.0);
   Cs.initialize(nXcc,nZcc,0.0);
   Vx.initialize(nXcc,nZcc,0.0);
   Vz.initialize(nXcc,nZcc,0.0);
   //
   FluxRx.initialize(nXce,nZcc,0.0);
   FluxLx.initialize(nXce,nZcc,0.0);
   FluxLimLx.initialize(nXce,nZcc,0.0);
   FluxLimRx.initialize(nXce,nZcc,0.0);
   //
   FluxRz.initialize(nXcc,nZce,0.0);
   FluxLz.initialize(nXcc,nZce,0.0);
   FluxLimLz.initialize(nXcc,nZce,0.0);
   FluxLimRz.initialize(nXcc,nZce,0.0);
   //
   FluxN_x.initialize(nXce,nZcc,0.0);
   FluxMx_x.initialize(nXce,nZcc,0.0);
   FluxMz_x.initialize(nXce,nZcc,0.0);
   FluxE_x.initialize(nXce,nZcc,0.0);
   //
   FluxN_z.initialize(nXcc,nZce,0.0);
   FluxMx_z.initialize(nXcc,nZce,0.0);
   FluxMz_z.initialize(nXcc,nZce,0.0);
   FluxE_z.initialize(nXcc,nZce,0.0);
   //
   //
   hy_cc.initialize(nXcc,nZcc,1.0);
   hy_ce.initialize(nXce,nZcc,1.0);
   hy_z.initialize(nXcc,nZce,1.0);
   hy_xz.initialize(nXce,nZce,1.0);
   //
   //
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      Json::Value geometry  = Phys.get("geometry",defValue);
      Json::Value advScheme = Phys.get("advScheme",defValue);
      Json::Value gammaVal  = Phys.get("gammaC",defValue);
      Json::Value NsubVal   = Phys.get("Nsub",defValue);
      if(advScheme == defValue || gammaVal == defValue) {
         cout << "ERROR: advScheme or gamma is " << endl;
         cout << "not declared in input file" << endl;
         exit (EXIT_FAILURE);
      } 
      
      geometry0 = geometry.asString();
      if(geometry0=="CAR" || geometry0=="CYL") {
         if(procID==0) {
            cout << "geometry is " << geometry0 << endl;
         }
         if(geometry0=="CYL") {
            for (auto i=0; i<nXce; i++) {
               for (auto j=0; j<nZce; j++) {
                  if(i<nXcc && j<nZcc) hy_cc(i,j) = Xgrid.Xcc.at(i);
                  if(j<nZcc) hy_ce(i,j)  = Xgrid.Xce2.at(i);
                  if(i<nXcc) hy_z(i,j) = Xgrid.Xcc.at(i);
                  hy_xz(i,j) = Xgrid.Xce2.at(i);
               }
            }
         }
      }
      else {
         cout << "valid geometry types are CAR and CYL" << endl;
         exit (EXIT_FAILURE);
      }

      
      advScheme0 = advScheme.asString();
      if(procID==0) {
         cout << "advection diff/interp scheme is " << advScheme0 << endl;
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
  
   //   get inlet parameters
   //
   const Json::Value Inlet = Phys.get("Inlet",defValue);
   if(Inlet.isObject()) { 
      Json::Value RinletVal = Inlet.get("Radius",defValue);
      Json::Value PinletVal = Inlet.get("Pressure",defValue);
      Json::Value MinletVal = Inlet.get("Mach",defValue);
      Json::Value TplenumVal = Inlet.get("Tplenum",defValue);
      Rinlet = RinletVal.asDouble(); 
      Pinlet = PinletVal.asDouble(); 
      Minlet = MinletVal.asDouble(); 
      Tplenum = TplenumVal.asDouble(); 
      if(procID==0) {
         cout << "Specified inlet parameters:" << endl;
         cout << "==========================" << endl;
         cout << "inlet R/R0   = " << Rinlet << endl;
         cout << "inlet P/P0  = " << Pinlet << endl;
         cout << "inlet Mach   = " << Minlet << endl;
         cout << "plenum T/T0  = " << Tplenum << endl;
      }
   
      // calculate plenum density and pressure
      //
      double exponent = gamma0/(gamma0-1.0);
      Pplenum = Pinlet*pow(1.0 + (gamma0-1.0)/2.0*Minlet*Minlet,exponent);
      Nplenum = Pplenum/Tplenum;
     
      // calculate all inlet exit values
      //
      Tinlet = Tplenum/(1.0 + (gamma0-1.0)/2.0*Minlet*Minlet);
      Ninlet = Pinlet/Tinlet;
      Vzinlet = Minlet*sqrt(gamma0*Pinlet/Ninlet);
      Mzinlet = Vzinlet*Ninlet;
      Einlet = Pinlet/(gamma0-1.0) + 0.5*Mzinlet*Mzinlet/Ninlet;
      
      if(procID==0) {
         cout << "Derived inlet parameters:" << endl;
         cout << "==========================" << endl;
         cout << "plenum P/P0  = " << Pplenum << endl;
         cout << "plenum N/N0  = " << Nplenum << endl;
         cout << "Tinlet = " << Tinlet << endl;
         cout << "Ninlet = " << Ninlet << endl;
      }

   } else {
      cout << "value for Physics variable \"Inlet\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   //   get initial profiles for variables
   //
   const Json::Value Nvar = Phys.get("N",defValue);
   if(Nvar.isObject()) { 
      Xgrid.setInitialProfile(N,Nvar);
      if(procID==0) Xgrid.setXminBoundary(N, 0.0, 1.0);   
      //if(procID==numProcs-1) setXmaxBoundary(N, 0.125);   
      if(procID==numProcs-1) Xgrid.setXmaxBoundary(N, 0.0, 1.0);   
      Xgrid.setZminBoundary(N,0.0,1.0);
      Xgrid.setZmaxBoundary(N,0.0,1.0);
      Xgrid.communicate(N);
      Nold  = N;
   } else {
      cout << "value for Physics variable \"N\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
 
   const Json::Value Vxvar = Phys.get("Vx",defValue);
   if(Vxvar.isObject()) { 
      Xgrid.setInitialProfile(Vx,Vxvar);
      if(procID==0) Xgrid.setXminBoundary(Vx, 0.0, 1.0);   
      if(procID==numProcs-1) Xgrid.setXmaxBoundary(Vx, 0.0, 1.0);   
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
      if(procID==0) Xgrid.setXminBoundary(Vz, 0.0, 1.0);   
      if(procID==numProcs-1) Xgrid.setXmaxBoundary(Vz, 0.0, 1.0);   
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
      if(procID==0) Xgrid.setXminBoundary(P, 0.0, 1.0);   
      //if(procID==numProcs-1) setXmaxBoundary(P, 0.1);   
      if(procID==numProcs-1) Xgrid.setXmaxBoundary(P, 0.0, 1.0);   
      Xgrid.communicate(P);
      Xgrid.setZminBoundary(P,0.0,1.0);
      Xgrid.setZmaxBoundary(P,0.0,1.0);
      E = P/(gamma0-1.0) + 0.5*(Mx*Vx+Mz*Vz);
      Eold  = E;
   } else {
      cout << "value for Physics variable \"P\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   Cs = sqrt(gamma0*P/N);
      

   // set inlet boundary conditions
   //
   Xgrid.setZboundaryInlet(N, 0,Rinlet,Ninlet); 
   Xgrid.setZboundaryInlet(Mz,0,Rinlet,Mzinlet); 
   Xgrid.setZboundaryInlet(E, 0,Rinlet,Einlet); 


   // need inital flux calculation before add to output
   // and before first call to Physics.advance   
   computeFluxes(Xgrid,2);   
  
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
   //
   dataFile.add(hy_cc,  "hy_cc",  0);  
   dataFile.add(hy_ce,  "hy_ce",  0);  

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
   double dX = Xgrid.dX;
   double dZ = Xgrid.dZ;
   double thisdt, thist;
   timeDomain* tmesh = timeDomain::tmesh;
   for (auto n=1; n<Nsub+1; n++) {
      thisdt = dt*n/Nsub;
      thist = tmesh->tSim + thisdt; // time at end of step
      for (auto i=nXg; i<nXcc-nXg; i++) {
         for (auto j=nZg; j<nZcc-nZg; j++) {
            N(i,j) = Nold(i,j) - thisdt*(FluxN_x(i+1,j) - FluxN_x(i,j))/hy_cc(i,j)/dX
                               - thisdt*(FluxN_z(i,j+1) - FluxN_z(i,j))/dZ;
            Mx(i,j)= Mxold(i,j)- thisdt*(FluxMx_x(i+1,j)-FluxMx_x(i,j))/hy_cc(i,j)/dX
                               - thisdt*(FluxMx_z(i,j+1)-FluxMx_z(i,j))/dZ;
            if(geometry0=="CYL") Mx(i,j) = Mx(i,j) + thisdt*P(i,j)/hy_cc(i,j);
            Mz(i,j)= Mzold(i,j)- thisdt*(FluxMz_x(i+1,j)-FluxMz_x(i,j))/hy_cc(i,j)/dX
                               - thisdt*(FluxMz_z(i,j+1)-FluxMz_z(i,j))/dZ;
            E(i,j) = Eold(i,j) - thisdt*(FluxE_x(i+1,j) - FluxE_x(i,j))/hy_cc(i,j)/dX
                               - thisdt*(FluxE_z(i,j+1) - FluxE_z(i,j))/dZ;
            
            if(N(i,j)<=Nthresh) N(i,j) = Nthresh;
            if(N(i,j)!=N(i,j)) {
               cout << "bout to go bad: N(i,j) = " << N(i,j) << endl;
               cout << "thist = " << thist << endl;
               cout << "X(i) = " << Xgrid.Xcc.at(i) << endl;
               cout << "Z(j) = " << Xgrid.Zcc.at(j) << endl;
               exit (EXIT_FAILURE);
            }
            Ethresh = N(i,j)/Nthresh*Pthresh/(gamma0-1.0)
                    + 0.5*(Mx(i,j)*Mx(i,j)+Mz(i,j)*Mz(i,j))/N(i,j);
            if(E(i,j)<=Ethresh) E(i,j) = Ethresh;
 
         }
      }

      if(procID==0) {
         Xgrid.setXminBoundary(N, 0.0, 1.0);   
         Xgrid.setXminBoundary(Mx, 0.0, -1.0);   
         Xgrid.setXminBoundary(Mz, 0.0, 1.0);   
         Xgrid.setXminBoundary(E, 0.0, 1.0);   
      }
      if(procID==numProcs-1) {
         Xgrid.setXmaxBoundary(N, 0.0, 1.0);   
         Xgrid.setXmaxBoundary(Mx, 0.0, -1.0);   
         Xgrid.setXmaxBoundary(Mz, 0.0, 1.0);   
         Xgrid.setXmaxBoundary(E, 0.0, 1.0);   
      }
      Xgrid.setZmaxBoundary(N, 0.0, 1.0);   
      Xgrid.setZmaxBoundary(Mx, 0.0, 1.0);   
      Xgrid.setZmaxBoundary(Mz, 0.0, -1.0);   
      Xgrid.setZmaxBoundary(E, 0.0, 1.0);   
      Xgrid.setZminBoundary(N, 0.0, 1.0);   
      Xgrid.setZminBoundary(Mx, 0.0, 1.0);   
      Xgrid.setZminBoundary(Mz, 0.0, -1.0);   
      Xgrid.setZminBoundary(E, 0.0, 1.0);   
      
      // set inlet boundary conditions
      //
      Xgrid.setZboundaryInlet(N, 0,Rinlet,Ninlet); 
      Xgrid.setZboundaryInlet(Mz,0,Rinlet,Mzinlet); 
      Xgrid.setZboundaryInlet(E, 0,Rinlet,Einlet); 
      
      Xgrid.communicate(N);
      Xgrid.communicate(Mx);
      Xgrid.communicate(Mz);
      Xgrid.communicate(E);


      //if(n==1 && Nsub==2) computeFluxes(Xgrid,2);
      //else computeFluxes(Xgrid,1);
      computeFluxes(Xgrid,2);

   }
 
   Nold  = N;
   Mxold = Mx;
   Mzold = Mz;
   Eold  = E;

} // end Physics.advance


void computeFluxes(const domainGrid& Xgrid, const int order)
{

   const int nXg = Xgrid.nXg;
   const int nZg = Xgrid.nZg;
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXce = Xgrid.nXce2;
   const int nZce = Xgrid.nZce2;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);


   matrix2D<double> Cspeed;
   matrix2D<double> FluxNcc_x, FluxMxcc_x, FluxMzcc_x, FluxEcc_x;
   matrix2D<double> FluxNcc_z, FluxMxcc_z, FluxMzcc_z, FluxEcc_z;
   matrix2D<double> Pce_z, Vz_z, FluxMzce_z;   

   FluxNcc_x.initialize(nXcc,nZcc,0.0);
   FluxMxcc_x.initialize(nXcc,nZcc,0.0);
   FluxMzcc_x.initialize(nXcc,nZcc,0.0);
   FluxEcc_x.initialize(nXcc,nZcc,0.0);
   //
   FluxNcc_z.initialize(nXcc,nZcc,0.0);
   FluxMxcc_z.initialize(nXcc,nZcc,0.0);
   FluxMzcc_z.initialize(nXcc,nZcc,0.0);
   FluxEcc_z.initialize(nXcc,nZcc,0.0);
   //
   Cspeed.initialize(nXcc,nZcc,0.0);
   //
   Pce_z.initialize(nXcc,nZce,0.0);
   Vz_z.initialize(nXcc,nZce,0.0);
   FluxMzce_z.initialize(nXcc,nZce,0.0);
   

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


   FluxNcc_x  = hy_cc*Mx;
   FluxMxcc_x = hy_cc*(Mx*Mx/N + P);
   FluxMzcc_x = hy_cc*Mz*Mx/N;
   FluxEcc_x  = hy_cc*Vx*(E + P);
   //
   FluxNcc_z  = Mz;
   FluxMxcc_z = Mx*Vz;
   FluxMzcc_z = Mz*Vz + P;
   FluxEcc_z  = Vz*(E + P);
   
   
   //  define correct FluxMz_z on zmin boundary
   //
   Xgrid.InterpToCellEdges(Pce_z,P,FluxMzce_z,"C2",1);
   Xgrid.InterpToCellEdges(Vz_z,Vz,FluxMzce_z,"C2",1);
   Xgrid.InterpToCellEdges(FluxMzce_z,Mz,FluxMzce_z,"C2",1);
   FluxMzce_z *=Vz_z;
   FluxMzce_z += Pce_z;

   // compute advective flux using
   // specified scheme from input file
   //
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN_x, FluxLx,FluxRx,FluxLimLx,FluxLimRx,
                           FluxNcc_x, Cspeed,hy_cc*N, "vanleer",0,order);
      Xgrid.computeFluxTVD(FluxMx_x,FluxLx,FluxRx,FluxLimLx,FluxLimRx,
                           FluxMxcc_x,Cspeed,hy_cc*Mx,"vanleer",0,order);
      Xgrid.computeFluxTVD(FluxMz_x,FluxLx,FluxRx,FluxLimLx,FluxLimRx,
                           FluxMzcc_x,Cspeed,hy_cc*Mz,"vanleer",0,order);
      Xgrid.computeFluxTVD(FluxE_x, FluxLx,FluxRx,FluxLimLx,FluxLimRx,
                           FluxEcc_x, Cspeed,hy_cc*E, "vanleer",0,order);
      //
      Xgrid.computeFluxTVD(FluxN_z, FluxLz,FluxRz,FluxLimLz,FluxLimRz,
                           FluxNcc_z, Cspeed,N, "vanleer",1,order);
      Xgrid.computeFluxTVD(FluxMx_z,FluxLz,FluxRz,FluxLimLz,FluxLimRz,
                           FluxMxcc_z,Cspeed,Mx,"vanleer",1,order);
      Xgrid.computeFluxTVD(FluxMz_z,FluxLz,FluxRz,FluxLimLz,FluxLimRz,
                           FluxMzcc_z,Cspeed,Mz,"vanleer",1,order);
      Xgrid.computeFluxTVD(FluxE_z, FluxLz,FluxRz,FluxLimLz,FluxLimRz,
                           FluxEcc_z, Cspeed,E, "vanleer",1,order);
   }      
   else {
      Xgrid.computeFluxTVDsimple(FluxN_x, FluxLx,FluxRx,FluxNcc_x, 
		                 Cspeed,hy_cc*N, advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxMx_x,FluxLx,FluxRx,FluxMxcc_x,
		                 Cspeed,hy_cc*Mx,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxMz_x,FluxLx,FluxRx,FluxMzcc_x,
		                 Cspeed,hy_cc*Mz,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxE_x, FluxLx,FluxRx,FluxEcc_x, 
		                 Cspeed,hy_cc*E, advScheme0,0);
      //
      Xgrid.computeFluxTVDsimple(FluxN_z, FluxLz,FluxRz,FluxNcc_z, 
		                 Cspeed,N, advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxMx_z,FluxLz,FluxRz,FluxMxcc_z,
		                 Cspeed,Mx,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxMz_z,FluxLz,FluxRz,FluxMzcc_z,
		                 Cspeed,Mz,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxE_z, FluxLz,FluxRz,FluxEcc_z, 
		                 Cspeed,E, advScheme0,1);
   } 
  
   
   // set BCs for fluxes
   //
   vector<double> P0;
   P0.assign(nZcc,0.0);
   if(procID==0) {
      Xgrid.setXminFluxBC(FluxN_x, 0.0, 0.0);
      for (auto j=0; j<nZcc; j++) {
         P0.at(j) = hy_ce(nXg,j)*(P(nXg,j)+P(nXg-1,j))/2.0;
      }
      Xgrid.setXminFluxBC(FluxMx_x, P0);
      Xgrid.setXminFluxBC(FluxMz_x, 0.0, 0.0);
      Xgrid.setXminFluxBC(FluxE_x, 0.0, 0.0);
   }
   if(procID==numProcs-1) {
      Xgrid.setXmaxFluxBC(FluxN_x, 0.0, 0.0);
      const int thisnX = P.size0();
      for (auto j=0; j<nZcc; j++) {
         P0.at(j) = hy_ce(thisnX-nXg-1,j)*(P(thisnX-nXg-1,j)+P(thisnX-nXg-2,j))/2.0;
      }
      Xgrid.setXmaxFluxBC(FluxMx_x, P0);
      Xgrid.setXmaxFluxBC(FluxMz_x, 0.0, 0.0);
      Xgrid.setXmaxFluxBC(FluxE_x, 0.0, 0.0);
   }

   // set BCs for fluxes on Z-boundaries
   //
   Xgrid.setZminFluxBC(FluxN_z,  FluxNcc_z);
   Xgrid.setZminFluxBC(FluxMx_z, FluxMxcc_z);
   Xgrid.setZminFluxBC(FluxMz_z, FluxMzce_z);
   Xgrid.setZminFluxBC(FluxE_z,  FluxEcc_z);
   //
   Xgrid.setZmaxFluxBC(FluxN_z,  0.0);
   Xgrid.setZmaxFluxBC(FluxMx_z, 0.0);
   Xgrid.setZmaxFluxBC(FluxMz_z, P);
   Xgrid.setZmaxFluxBC(FluxE_z,  0.0);
   

} // end computeFluxes


void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid, const int verbose)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   matrix2D<double> Cchar(Cs);
   double Cmax;
   //Cchar = abs(Vx)+Cs;
   Cchar += sqrt(Vx*Vx+Vz*Vz);
   Cmax = max(Cchar);
   //cout << "Cmax = " << Cmax << endl;

   const double dX = Xgrid.dX;
   const double dZ = Xgrid.dZ;
   
   double dtmax = 0.5*dX*dZ/(dX+dZ)/Cmax;
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);
   if(procID==0 && verbose) {
      cout << "max stable time step is " << dtmax << endl;
   }
}

