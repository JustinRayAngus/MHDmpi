/***
 * 
 * physics module for 2D zpinch with m=0 using Ideal MHD
 * using energy conservative methods
 *
 * dN/dt  + d(r*Mr)/r/dr + d(Mz)/dz = 0
 * dMr/dt + d(r*(Mr*Vr + P + B^2/2))/r/dr -(P-B^2/2)/r + d(Mr*Vz)/dz = 0
 * dMz/dt + d(Mz*Vz + P + B^2/2)/dz + d(r*Mz*Vr)/r/dr = 0
 * dE/dt  + d(rVr*(E+P))/r/dr +d(Vz*(E+P))/dz = Jz*Ez + Jr*Er  
 * dBy/dt + d(Vr*By)/dr + d(Vz*By)/dr = 0
 *
 * Ez = -Vr*By
 * Er = Vz*By
 * Jz = d(r*By)/r/dr
 * Jr = -d(By)/dz
 *
 * Vr = Mr/N
 * Vz = Mz/N
 * P = (E-0.5*N*(Vr^2 + Vz^2))*(gamma0-1)
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
double Nthresh;        // density threshold
int Nsub;             // time-solver subcycle steps
matrix2D<double> N, Mx, Mz, E, By, Ez, Ex;  // time-evolving variables
matrix2D<double> P0, deltaP;                    // initial perturbation
matrix2D<double> eta, Cs, Vx, Vz, P, T, S;     // derived variables
matrix2D<double> Jz, Jz0, Jx, Jx0;          // derived variables
matrix2D<double> eta_x, eta_z, Jzcc, Ezcc, Jxcc, VxBy_x, VzBy_z;
matrix2D<double> Nold, Mxold, Mzold, Eold, Byold;
matrix2D<double> Exold, Exhalf, Ezold, Ezhalf;
matrix2D<double> FluxR_x, FluxL_x, FluxRatio_x, FluxLim_x;
matrix2D<double> FluxR_z, FluxL_z, FluxRatio_z, FluxLim_z;  
matrix2D<double> FluxN_x, FluxE_x, FluxMx_x, FluxMz_x, FluxBy_x, FluxEz_x;
matrix2D<double> FluxN_z, FluxE_z, FluxMx_z, FluxMz_z, FluxBy_z, FluxEx_z;

matrix2D<double> Fx, rcc, rce_x;

// Set lower threshold values for N, T, and S
double Tthresh=2.5e-2, Pthresh, Ethresh;

void computeFluxes(const domainGrid&, const int);
void setXminBoundary(matrix2D<double>&, const double, const double);
void setXminExtrap(matrix2D<double>&, const int);
void setXmaxExtrap(matrix2D<double>&, const int);
void setXmaxBy(matrix2D<double>&, const double);
void setXminBoundaryEz(matrix2D<double>&);
void setXmaxBoundary(matrix2D<double>&, const double, const double);
void setZboundaryPeriodic(matrix2D<double>&);

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
   const int nXce = Xgrid.Xce.size();
   const int nZcc = Xgrid.Zcc.size();
   const int nZce = Xgrid.Zce.size();


   N.initialize(nXcc,nZcc,0.0);
   Nold.initialize(nXcc,nZcc,0.0);
   deltaP.initialize(nXcc,nZcc,0.0);
   P0.initialize(nXcc,nZcc,0.0);
   Mx.initialize(nXcc,nZcc,0.0);
   Mxold.initialize(nXcc,nZcc,0.0);
   Mz.initialize(nXcc,nZcc,0.0);
   Mzold.initialize(nXcc,nZcc,0.0);
   E.initialize(nXcc,nZcc,0.0);
   Eold.initialize(nXcc,nZcc,0.0);
   By.initialize(nXcc,nZcc,0.0);
   Byold.initialize(nXcc,nZcc,0.0);
   P.initialize(nXcc,nZcc,0.0);
   T.initialize(nXcc,nZcc,0.0);
   Cs.initialize(nXcc,nZcc,0.0);
   Vx.initialize(nXcc,nZcc,0.0);
   Vz.initialize(nXcc,nZcc,0.0);
   
   Jz.initialize(nXcc,nZcc,0.0);
   Ez.initialize(nXcc,nZcc,0.0);
   //
   Jx.initialize(nXcc,nZcc,0.0); 
   Ex.initialize(nXcc,nZcc,0.0);
   //
   S.initialize(nXcc,nZcc,0.0);
   //
   FluxN_x.initialize(nXce, nZcc, 0.0);
   FluxMx_x.initialize(nXce,nZcc, 0.0);
   FluxMz_x.initialize(nXce,nZcc, 0.0);
   FluxE_x.initialize(nXce, nZcc, 0.0);
   FluxBy_x.initialize(nXce,nZcc, 0.0);
   
   FluxN_z.initialize(nXcc, nZce, 0.0);
   FluxMx_z.initialize(nXcc,nZce, 0.0);
   FluxMz_z.initialize(nXcc,nZce, 0.0);
   FluxE_z.initialize(nXcc, nZce, 0.0);
   FluxBy_z.initialize(nXcc,nZce, 0.0);
 
   FluxRatio_x.initialize(nXce,nZcc,0.0);
   FluxLim_x.initialize(nXce,nZcc,0.0);
   FluxR_x.initialize(nXce,nZcc,0.0);
   FluxL_x.initialize(nXce,nZcc,0.0);
   //
   FluxRatio_z.initialize(nXcc,nZce,0.0);
   FluxLim_z.initialize(nXcc,nZce,0.0);
   FluxR_z.initialize(nXcc,nZce,0.0);
   FluxL_z.initialize(nXcc,nZce,0.0);
   
   
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      Json::Value advScheme = Phys.get("advScheme",defValue);
      Json::Value gammaVal  = Phys.get("gammaC",defValue);
      Json::Value NsubVal   = Phys.get("Nsub",defValue);
      Json::Value NthreshVal= Phys.get("Nthresh",defValue);
      if(advScheme == defValue || gammaVal == defValue ||
	 NsubVal == defValue || NthreshVal == defValue) {
         cout << "ERROR: advScheme or gamma " << endl;
         cout << "or Nsub or Nthresh is " << endl;
         cout << "not declared in input file" << endl;
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

      Nthresh = NthreshVal.asDouble();
      if(procID==0) cout << "Nthresh = " << Nthresh << endl;
      Pthresh = Nthresh*Tthresh;
      //Sthresh = Pthresh/pow(Nthresh,gamma0-1);
      Ethresh = Pthresh/(gamma0-1.0);
      

   }
   else {
      cout << "value for key \"Physics\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   //   get initial profiles for variables
   //
   const Json::Value Pvar = Phys.get("P",defValue);
   if(Pvar.isObject()) { 
      Xgrid.setInitialProfile(P,Pvar);
      if(procID==0) setXminExtrap(P, 0);   
      if(procID==numProcs-1) setXmaxExtrap(P,0); 
      setZboundaryPeriodic(P);  
      Xgrid.communicate(P);
   } else {
      cout << "value for Physics variable \"P\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
 
   const Json::Value Tvar = Phys.get("T",defValue);
   if(Tvar.isObject()) { 
      Xgrid.setInitialProfile(T,Tvar);
      //if(procID==0) setXminBoundary(T, 0.0, 1.0);   
      //if(procID==numProcs-1) setXmaxBoundary(T, 0.0, 1.0); 
      setZboundaryPeriodic(T);  
      Xgrid.communicate(T);
   } else {
      cout << "value for Physics variable \"T\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   
   const Json::Value Vzvar = Phys.get("Vz",defValue);
   if(Vzvar.isObject()) { 
      Xgrid.setInitialProfile(Vz,Vzvar);
      if(procID==0) setXminExtrap(Vz, 0);   
      if(procID==numProcs-1) setXmaxExtrap(Vz, 0); 
      setZboundaryPeriodic(Vz);  
      Xgrid.communicate(Vz);
   } else {
      cout << "value for Physics variable \"Vz\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }


   const Json::Value Byvar = Phys.get("By",defValue);
   if(Byvar.isObject()) { 
      Xgrid.setInitialProfile(By,Byvar);
      if(procID==0) setXminBoundary(By,0.0,-1.0);   
      if(procID==numProcs-1) setXmaxBy(By,-1);
      //if(procID==numProcs-1) setXmaxExtrap(By,2);
      setZboundaryPeriodic(By);  
      Xgrid.communicate(By);
      Byold = By;
   } else {
      cout << "value for Physics variable \"By\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   
   const Json::Value deltaPvar = Phys.get("deltaP",defValue);
   if(deltaPvar.isObject()) { 
      Xgrid.setInitialProfile(deltaP,deltaPvar);
      if(procID==0) setXminExtrap(deltaP, 0);   
      if(procID==numProcs-1) setXmaxExtrap(deltaP,0); 
      setZboundaryPeriodic(deltaP);  
      Xgrid.communicate(deltaP);
   } else {
      cout << "value for Physics variable \"deltaP\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   // create a 2D matrix for r at cell center and at cell edges
   //
   Fx.initialize(nXcc,nZcc,0.0);
   rcc.initialize(nXcc,nZcc,0.0);
   rce_x.initialize(nXce,nZcc,0.0);
   for (auto i=0; i<nXcc; i++) {
      for (auto j=0; j<nZcc; j++) {
	 rcc(i,j) = Xgrid.Xcc.at(i);
	 if(i<nXce) rce_x(i,j) = Xgrid.Xce.at(i);
      }
   }
   
   // calculate current density 
   //
   Xgrid.DDX(Jz,rcc*By); 
   Jz /= rcc;
   setXmaxBoundary(Jz,0.0,0.0);
   //setXminBoundary(Jz,sqrt(2.0)*2.0/0.3333,0.0);


   // add perturbation to pressure profile
   //
   P0 = P;
   P *= 1.0+deltaP;
   double minP = min(P);
   if(minP<Nthresh) P += Nthresh; // assumed total T=Ti+Te=1
   N = P/T;
   Nold = N;
 
   S = P/pow(N,gamma0-1.0);
   Cs = pow(gamma0*P/N + By*By/N,0.5);

   Mz = N*Vz;
   Mzold = N*Vz;
   Ex = Vz*By;

   E = 0.5*(Mx*Mx+Mz*Mz)/N + P/(gamma0-1.0);   
   Eold = E;

   computeFluxes(Xgrid, 2); // 1 for first order calculation   


   // add stuff to output files
   //
   dataFile.add(rcc, "rcc", 0);  // force density x-direction 
   dataFile.add(Fx, "Fx", 1);  // force density x-direction 
   dataFile.add(N,  "N",  1);  // density 
   dataFile.add(deltaP,  "deltaP",  0);  // density perturbation 
   dataFile.add(P0,  "P0",  0);  // initial unperturbed pressure profile 
   dataFile.add(Mx, "Mx", 1);  // momentum density 
   dataFile.add(Mz, "Mz", 1);  // momentum density 
   dataFile.add(E,  "E",  1);  // total fluid energy
   dataFile.add(S,  "S",  1);  // entropy density
   dataFile.add(By, "By", 1);  // magnetic field
   dataFile.add(P, "P", 1);  // pressure
   dataFile.add(T,  "T",  1);  // temperature
   dataFile.add(Vx, "Vx", 1);  // x-velocity
   dataFile.add(Vz, "Vz", 1);  // z-velocity
   dataFile.add(Jz, "Jz", 1);  // current density
   dataFile.add(Ez, "Ez", 1);  // z-electric field
   dataFile.add(Jx, "Jx", 1);  // current density
   dataFile.add(Ex, "Ex", 1);  // x-electric field
   dataFile.add(Cs, "Cs", 1);  // sound speed
   dataFile.add(gamma0,"gamma0",0); 
   //
   dataFile.add(FluxN_x,  "FluxN_x",  1);  
   dataFile.add(FluxMx_x, "FluxMx_x", 1);  
   dataFile.add(FluxMz_x, "FluxMz_x", 1);  
   dataFile.add(FluxE_x,  "FluxE_x",  1);  
   dataFile.add(FluxBy_x, "FluxBy_x", 1);  
   //
   dataFile.add(FluxRatio_x, "FluxRatio_x", 1);  
   dataFile.add(FluxLim_x, "FluxLim_x", 1);  
   dataFile.add(FluxR_x, "FluxR_x", 1);
   dataFile.add(FluxL_x, "FluxL_x", 1);

} // end Physics.initilize


void Physics::advance(const domainGrid& Xgrid, const double dt)
{
   const int nXcc = Xgrid.nXcc;
   const int nXg = Xgrid.nXg;
   const int nZcc = Xgrid.nZcc;
   const int nZg = Xgrid.nZg;
   
   vector<double> Xcc, Xce;
   Xcc = Xgrid.Xcc;
   Xce = Xgrid.Xce;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   //domainGrid* mesh = domainGrid::mesh;

   // Explicit forward advance from n to n+1 using subcycling in time 
   //  
   double thisdt;
   for (auto n=1; n<Nsub+1; n++) {
      thisdt = dt*n/Nsub;

      // Update N, M, S, and B terms from n to n+1
      for (auto i=nXg; i<nXcc-nXg; i++) {
         for (auto j=nZg; j<nZcc-nZg; j++) {
	 
         Fx(i,j) = - (FluxMx_x(i,j)-FluxMx_x(i-1,j))/Xcc.at(i)/Xgrid.dX
		   + (P(i,j)-By(i,j)*By(i,j)/2.0)/Xcc.at(i);
	 if(procID==numProcs-1 && i==nXcc-nXg-1) Fx(i,j) = Fx(i-1,j)/3.0;
	 //if(i==nXcc-nXg-1) Fx(i,j) = 2.0*Fx(i-1,j)-Fx(i-2,j); // doesn't work well

	 N(i,j)  = Nold(i,j)  - thisdt*(FluxN_x(i,j)  - FluxN_x(i-1,j))/Xcc.at(i)/Xgrid.dX
	                      - thisdt*(FluxN_z(i,j)  - FluxN_z(i,j-1))/Xgrid.dZ;
         Mx(i,j) = Mxold(i,j) + thisdt*Fx(i,j)
                              - thisdt*(FluxMx_z(i,j) - FluxMx_z(i,j-1))/Xgrid.dZ;
         Mz(i,j) = Mzold(i,j) - thisdt*(FluxMz_x(i,j) - FluxMz_x(i-1,j))/Xcc.at(i)/Xgrid.dX
                              - thisdt*(FluxMz_z(i,j) - FluxMz_z(i,j-1))/Xgrid.dZ;
         E(i,j)  = Eold(i,j)  - thisdt*(FluxE_x(i,j)  - FluxE_x(i-1,j))/Xcc.at(i)/Xgrid.dX
                              - thisdt*(FluxE_z(i,j)  - FluxE_z(i,j-1))/Xgrid.dZ
                              + thisdt*(Jx(i,j)*Ex(i,j) + Jz(i,j)*Ez(i,j));
	 By(i,j) = Byold(i,j) - thisdt*(FluxBy_x(i,j) - FluxBy_x(i-1,j))/Xgrid.dX
	                      - thisdt*(FluxBy_z(i,j) - FluxBy_z(i,j-1))/Xgrid.dZ;
	 
	 if(N(i,j)<=Nthresh) N(i,j) = Nthresh;
	 //if(M(i,j)<0.0) M(i,j) = 0.0;
	 Ethresh = Pthresh/(gamma0-1.0) + 0.5*(Mx(i,j)*Mx(i,j)+Mz(i,j)*Mz(i,j))/N(i,j);
	 if(E(i,j)<=Ethresh) E(i,j) = Ethresh;

         }
      }

      timeDomain* tmesh = timeDomain::tmesh;
      //cout << "thist = " << tmesh->tSim << endl;
      //double thist = tmesh->tSim;
      
      if(procID==0) {
         setXminExtrap(N, 0);   
         setXminBoundary(Mx,0.0,-1.0);   
         setXminExtrap(Mz,0);   
         setXminExtrap(E,0);
         setXminBoundary(By,0.0,-1.0);   
      }
     
      if(procID==numProcs-1) {
         setXmaxBoundary(Mx, 0.0, -1.0);   
         setXmaxExtrap(N, 0);   
         setXmaxExtrap(Mz, 0);   
         setXmaxExtrap(E, 0);   
         setXmaxBy(By,-1);
      }
      setZboundaryPeriodic(N);
      setZboundaryPeriodic(Mx);
      setZboundaryPeriodic(Mz);
      setZboundaryPeriodic(E);
      setZboundaryPeriodic(By);
      Xgrid.communicate(Fx);
      Xgrid.communicate(N);
      Xgrid.communicate(Mx);
      Xgrid.communicate(Mz);
      Xgrid.communicate(E);
      Xgrid.communicate(By);

      // compute fluxes using partialy updated fields at n+1/2 for 2nd order
      // or n+1 for 1st order
      computeFluxes(Xgrid, Nsub); // compute second order fluxes at n+1/2

   } // finish subcycle steps for N, M, S, and B

   // update old fields
   Nold = N;
   Mxold = Mx;
   Mzold = Mz;
   Eold = E;
   Byold = By;
   
   // compute fluxes using fully updated fields at n+1
   computeFluxes(Xgrid, 2);

} // end Physics.advance


void computeFluxes(const domainGrid& Xgrid, const int order)
{

   const int nXcc = Xgrid.nXcc;
   //const int nXce = Xgrid.nXce;
   const int nZcc = Xgrid.nZcc;
   //const int nZce = Xgrid.nZce;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   matrix2D<double> FluxNcc_x, FluxMxcc_x, FluxMzcc_x, FluxEcc_x, FluxBycc_x;
   matrix2D<double> FluxNcc_z, FluxMxcc_z, FluxMzcc_z, FluxEcc_z, FluxBycc_z;
   matrix2D<double> CspeedBx, CspeedBz, Cspeedx, Cspeedz; 
 
   CspeedBx.initialize(nXcc,nZcc,0.0);
   CspeedBz.initialize(nXcc,nZcc,0.0);
   Cspeedx.initialize(nXcc,nZcc,0.0);
   Cspeedz.initialize(nXcc,nZcc,0.0);
   FluxNcc_x.initialize(nXcc,nZcc,0.0);
   FluxNcc_z.initialize(nXcc,nZcc,0.0);
   FluxMxcc_x.initialize(nXcc,nZcc,0.0);
   FluxMxcc_z.initialize(nXcc,nZcc,0.0);
   FluxMzcc_x.initialize(nXcc,nZcc,0.0);
   FluxMzcc_z.initialize(nXcc,nZcc,0.0);
   FluxEcc_x.initialize(nXcc,nZcc,0.0);
   FluxEcc_z.initialize(nXcc,nZcc,0.0);
   FluxBycc_x.initialize(nXcc,nZcc,0.0);
   FluxBycc_z.initialize(nXcc,nZcc,0.0);

   //  define derived variables
   //
   Vx = Mx/N;
   Vz = Mz/N;
   P  = (E - 0.5*(Vx*Mx+Vz*Mz))*(gamma0-1.0);
   T  = P/N;
   S = P/pow(N,gamma0-1.0);
   Cs = sqrt(gamma0*P/N + By*By/N);

   Xgrid.DDX(Jz,rcc*By); 
   Jz /= rcc;
   setXmaxBoundary(Jz,0.0,0.0);
   Xgrid.communicate(Jz);
   
   Xgrid.DDZ(Jx,By);
   Xgrid.communicate(Jx);
   Jx *= -1.0; 
   
   Ez = -Vx*By;
   Ex = Vz*By;

   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeedx = sqrt(Vx*Vx);
   Cspeedz = sqrt(Vz*Vz);
   Cspeedx = sqrt(Vx*Vx) + Cs;
   Cspeedz = sqrt(Vz*Vz) + Cs; 
   
   FluxNcc_x = rcc*Mx;
   FluxNcc_z = Mz;
   FluxMxcc_x = rcc*(Mx*Vx + P + By*By/2.0);
   FluxMxcc_z = Mx*Vz;
   FluxMzcc_x = rcc*Mz*Vx;
   FluxMzcc_z = Mz*Vz + P + By*By/2.0;
   FluxEcc_x = rcc*Vx*(E+P);
   FluxEcc_z = Vz*(E+P);
   FluxBycc_x = Vx*By;
   FluxBycc_z = Vz*By;

   
   // compute advective flux using
   // specified scheme from input file
   //
   Nsub = 2; // hard code this for force balance at initiation
   if(advScheme0 == "TVD") {
      //Xgrid.computeFluxTVDnew(FluxN_x, FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
      //                     FluxNcc_x, Cspeedx,rcc*N, 0,Nsub);
      Xgrid.computeFluxTVD(FluxN_x, FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxNcc_x, Cspeedx,rcc*N, "vanleer",0,Nsub);
      Xgrid.computeFluxTVD(FluxMx_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMxcc_x,Cspeedx,rcc*Mx,"vanleer",0,Nsub);
      Xgrid.computeFluxTVD(FluxMz_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMzcc_x,Cspeedx,rcc*Mz,"vanleer",0,Nsub);
      Xgrid.computeFluxTVD(FluxE_x, FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxEcc_x, Cspeedx,rcc*E, "vanleer",0,Nsub);
      Xgrid.computeFluxTVD(FluxBy_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxBycc_x,Cspeedx,By,    "vanleer",0,Nsub);
      //
      //
      Xgrid.computeFluxTVD(FluxN_z, FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxNcc_z, Cspeedz,N, "vanleer",1,Nsub);
      Xgrid.computeFluxTVD(FluxMx_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxMxcc_z,Cspeedz,Mx,"vanleer",1,Nsub);
      Xgrid.computeFluxTVD(FluxMz_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxMzcc_z,Cspeedz,Mz,"vanleer",1,Nsub);
      Xgrid.computeFluxTVD(FluxE_z, FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxEcc_z, Cspeedz,E, "vanleer",1,Nsub);
      Xgrid.computeFluxTVD(FluxBy_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxBycc_z,Cspeedz,By,"vanleer",1,Nsub);
   }
   else {
      Xgrid.computeFluxTVDsimple(FluxN_x, FluxL_x,FluxR_x,FluxNcc_x,
                                 Cspeedx,rcc*N, advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxMx_x,FluxL_x,FluxR_x,FluxMxcc_x,
                                 Cspeedx,rcc*Mx,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxMz_x,FluxL_x,FluxR_x,FluxMzcc_x,
                                 Cspeedx,rcc*Mz,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxE_x, FluxL_x,FluxR_x,FluxEcc_x,
                                 Cspeedx,rcc*E, advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxBy_x,FluxL_x,FluxR_x,FluxBycc_x,
                                 Cspeedx,By,    advScheme0,0);
      //
      Xgrid.computeFluxTVDsimple(FluxN_z, FluxL_z,FluxR_z,FluxNcc_z,
                                 Cspeedz,N, advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxMx_z,FluxL_z,FluxR_z,FluxMxcc_z,
                                 Cspeedz,Mx,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxMz_z,FluxL_z,FluxR_z,FluxMzcc_z,
                                 Cspeedz,Mz,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxE_z, FluxL_z,FluxR_z,FluxEcc_z,
                                 Cspeedz,E, advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxBy_z,FluxL_z,FluxR_z,FluxBycc_z,
                                 Cspeedz,By,advScheme0,1);
   } 

   if(procID==0)  {
      setXminBoundary(FluxN_x,  0.0, 0.0);   
      setXminBoundary(FluxMx_x, 0.0, 0.0);   
      setXminBoundary(FluxMz_x, 0.0, 0.0);   
      setXminBoundary(FluxE_x,  0.0, 0.0);
      setXminBoundary(FluxBy_x, 0.0, 0.0);
      
      //setXminBoundary(FluxMx_z, 0.0, 0.0);   
   }   
   if(procID==numProcs-1)  {
      setXmaxBoundary(FluxN_x,  0.0, 0.0);   
      //setXmaxBoundary(FluxMx, (P.at(2)+P.at(1))/2.0, 0.0);   
      setXmaxBoundary(FluxMz_x, 0.0, 0.0);   
      setXmaxBoundary(FluxE_x,  0.0, 0.0);
      setXmaxBoundary(FluxBy_x, 0.0, 0.0);
      
      //setXmaxBoundary(FluxMx_z, 0.0, 0.0);   
   }   

   
   Xgrid.communicate(FluxL_x);   
   Xgrid.communicate(FluxR_x);   
   Xgrid.communicate(FluxRatio_x);   
   Xgrid.communicate(FluxLim_x);   
   
    
} // end computeFluxes


void setXminBoundary(matrix2D<double>& var, 
		     const double C0, const double C1)
{
   const int thisnZ = var.size1();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift = mesh->nXg;

   for (auto i=0; i<ishift; i++) {
      for (auto j=0; j<thisnZ; j++) {
         var(ishift-i-1,j) = C0 + C1*var(ishift+i,j);
      }
   }

}

void setXminExtrap(matrix2D<double>& var, const int order)
{
   
   const int thisnZ = var.size1();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift = mesh->nXg;

   for (auto j=0; j<thisnZ; j++) {
      if(order==0) {
         var(1,j) = var(2,j);
         var(0,j) = var(2,j);
      }
      else {
         //var(1,j) = 2.0*var(2,j) - var(3,j);
         //var(0,j) = 2.0*var(1,j) - var(2,j);
         var(1,j) = 3.0*(var(2,j) - var(3,j)) + var(4,j);
         var(0,j) = 3.0*(var(1,j) - var(2,j)) + var(3,j);
      }
   }

}

void setXmaxExtrap(matrix2D<double>& var, const int order)
{
   
   const int thisnX = var.size0();
   const int thisnZ = var.size1();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift = thisnX-mesh->nXg;

   for (auto j=0; j<thisnZ; j++) {
      if(order==0) {
         var(ishift,j) = var(ishift-1,j);
         var(ishift+1,j) = var(ishift-1,j);
      }
      else {
         //var(ishift,j)   = 2.0*var(ishift-1,j) - var(ishift-2,j);
         //var(ishift+1,j) = 2.0*var(ishift,j)   - var(ishift-1,j);
         var(ishift,j)   = 3.0*(var(ishift-1,j) - var(ishift-2,j)) + var(ishift-3,j);
         var(ishift+1,j) = 3.0*(var(ishift,j) - var(ishift-1,j)) + var(ishift-2,j);
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
   
   for (auto i=ishift; i<thisnX; i++) {
      for (auto j=0; j<thisnZ; j++) {
         var(i,j) = C0 + C1*var(2*ishift-i-1,j);
      }
   }
      
}

void setXmaxBy(matrix2D<double>& var, const double exponent)
{
   const int thisnX = var.size0(); 
   const int thisnZ = var.size1(); 

   domainGrid* mesh = domainGrid::mesh;
   const int ishift = thisnX-mesh->nXg;
   vector<double> Xcc = mesh->Xcc;
   
   for (auto i=ishift; i<thisnX; i++) {
      for (auto j=0; j<thisnZ; j++) {
         var(i,j) = var(i-1,j)*pow(Xcc.at(i)/Xcc.at(i-1),exponent);
      }
   }
      
}

void setXminBoundaryEz(matrix2D<double>& var)
{
   //const int thisnX = var.size0(); 
   const int thisnZ = var.size1(); 

   domainGrid* mesh = domainGrid::mesh;
   const int ishift = mesh->nXg;
   
   for (auto j=0; j<thisnZ; j++) {
      var(ishift-1,j) = 4.0/3.0*var(ishift,j) - 1.0/3.0*var(ishift+1,j);
   }
      

}

void setZboundaryPeriodic(matrix2D<double>& var)
{
   
   const int thisnX = var.size0();
   const int thisnZ = var.size1();

   domainGrid* mesh = domainGrid::mesh;
   const int nZg = mesh->nZg;
 
   //assert(jshift==2 || jshift==1);

   for (auto i=0; i<thisnX; i++) {
      for (auto j=0; j<nZg; j++) {
         var(i,j) = var(i,thisnZ-2*nZg+j);
         var(i,thisnZ-nZg+j) = var(i,nZg+j);
            
         //var(i,1) = var(i,thisnZ-3);
         //var(i,0) = var(i,thisnZ-4);

         //var(i,thisnZ-2) = var(i,2);
         //var(i,thisnZ-1) = var(i,3);
      }
   }

}


void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid, const int verbose)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   matrix2D<double> Cchar_x(Cs), Cchar_z(Cs);
   double Cmax_x, Cmax_z, Cmax;
   Cchar_x += abs(Vx);
   Cchar_z += abs(Vz);
   Cmax_x = max(Cchar_x);
   Cmax_z = max(Cchar_z);
   Cmax = max(Cmax_x,Cmax_z);
   //cout << "Cmax = " << Cmax << endl;

   const double dX = Xgrid.dX;
   const double dZ = Xgrid.dZ;
   
   double dtCFL_sound = 0.5*dX*dZ/(dX+dZ)/Cmax;
   double dtmax = dtCFL_sound;
   //double dtmax = dtCFL_sound;
   
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);


   //dtSim = 1.0e-4;
   if(procID==0 && verbose) {
      cout << "dtSim = " << dtSim << endl;
   }
}

