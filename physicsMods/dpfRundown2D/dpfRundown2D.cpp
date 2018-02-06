/***
 * 
 * physics module for 2D resistive MHD dpf rundown
 * (1D railgun)
 *
 * dN/dt  + d(Mx)/dx + d(Mz)/dz = 0
 * dMx/dt + d(Mx*Ux + P)/dx + d(Mx*Uz)/dz = -Jz*By
 * dMz/dt + d(Mz*Uz + P)/dz + d(Mz*Ux)/dx = Jx*By
 * dS/dt  + d(S*Ux)/dx +d(S*Uz)/dz = eta*(gamma-1)/N^(gamma-1)*J^2  
 * dBy/dt + d(-Ex)/dx + d(Ez)/dx = 0
 * dEz/dt + d(-By)/dx = -Jz
 * dEx/dt + d(By)/dz  = -Jx
 *
 * Jx = [Ex - Uz*By]/eta
 * Jz = [Ez + Ux*By]/eta
 * S  = P/N^(gamma-1)
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
double eta0=0.08;     // resistivity (may need to decrease dt and or dx more if eta0 really low)
                      // I think the issue is boundary related? May need vacuum resistivity
double delta0=1.0e-4; // relaxation const (v/c)^2
double B0 = 0.0;      // boundary value of magnetic field
int Nsub;             // time-solver subcycle steps
matrix2D<double> N, Mx, Mz, S, By, Ez, Ex;  // time-evolving variables
matrix2D<double> eta, Cs, Vx, Vz, P, T;     // derived variables
matrix2D<double> Jz, Jz0, Jx, Jx0;          // derived variables
matrix2D<double> eta_x, eta_z, Jzcc, Jxcc, VxBy_x, VzBy_z;
matrix2D<double> Nold, Mxold, Mzold, Sold, Byold;
matrix2D<double> Exold, Exhalf, Ezold, Ezhalf;
matrix2D<double> FluxR_x, FluxL_x, FluxRatio_x, FluxLim_x;
matrix2D<double> FluxR_z, FluxL_z, FluxRatio_z, FluxLim_z;  
matrix2D<double> FluxN_x, FluxS_x, FluxMx_x, FluxMz_x, FluxBy_x, FluxEz_x;
matrix2D<double> FluxN_z, FluxS_z, FluxMx_z, FluxMz_z, FluxBy_z, FluxEx_z;

// Set lower threshold values for N, T, and S
double Nthresh=1.0e-4, Tthresh=2.5e-2, Sthresh;

void computeFluxes(const domainGrid&, const int);
void setXminBoundary(matrix2D<double>&, const double, const double);
void setXminExtrap(matrix2D<double>&, const int);
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
   Mx.initialize(nXcc,nZcc,0.0);
   Mxold.initialize(nXcc,nZcc,0.0);
   Mz.initialize(nXcc,nZcc,0.0);
   Mzold.initialize(nXcc,nZcc,0.0);
   S.initialize(nXcc,nZcc,0.0);
   Sold.initialize(nXcc,nZcc,0.0);
   By.initialize(nXcc,nZcc,0.0);
   Byold.initialize(nXcc,nZcc,0.0);
   P.initialize(nXcc,nZcc,0.0);
   T.initialize(nXcc,nZcc,0.0);
   eta.initialize(nXcc,nZcc,0.0);
   Cs.initialize(nXcc,nZcc,0.0);
   Vx.initialize(nXcc,nZcc,0.0);
   Vz.initialize(nXcc,nZcc,0.0);
   
   Jz.initialize(nXce,nZcc,0.0); // J defined at cell edges
   Jzcc.initialize(nXcc,nZcc,0.0);
   eta_x.initialize(nXce,nZcc,0.0);
   VxBy_x.initialize(nXce,nZcc,0.0);
   Jz0.initialize(nXce,nZcc,0.0);
   //
   Jx.initialize(nXcc,nZce,0.0); // J defined at cell edges
   Jxcc.initialize(nXcc,nZcc,0.0);
   eta_z.initialize(nXcc,nZce,0.0);
   VzBy_z.initialize(nXcc,nZce,0.0);
   Jx0.initialize(nXcc,nZce,0.0);
   // Ez is defined on cell edges in x
   Ez.initialize(nXce,nZcc,0.0);
   Ezold.initialize(nXce,nZcc,0.0);
   Ezhalf.initialize(nXce,nZcc,0.0);
   // Ex is defined on cell edges in z
   Ex.initialize(nXcc,nZce,0.0);
   Exold.initialize(nXcc,nZce,0.0);
   Exhalf.initialize(nXcc,nZce,0.0);

   //
   FluxN_x.initialize(nXce,nZcc,0.0);
   FluxMx_x.initialize(nXce,nZcc,0.0);
   FluxMz_x.initialize(nXce,nZcc,0.0);
   FluxS_x.initialize(nXce,nZcc,0.0);
   FluxBy_x.initialize(nXce,nZcc,0.0);
   FluxEz_x.initialize(nXcc,nZcc,0.0); // Flux for Ez on cell-center
   
   FluxN_z.initialize(nXcc,nZce,0.0);
   FluxMx_z.initialize(nXcc,nZce,0.0);
   FluxMz_z.initialize(nXcc,nZce,0.0);
   FluxS_z.initialize(nXcc,nZce,0.0);
   FluxBy_z.initialize(nXcc,nZce,0.0);
   FluxEx_z.initialize(nXcc,nZcc,0.0); // Flux for Ex on cell-center
 
 
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
      Json::Value deltaVal  = Phys.get("delta0",defValue);
      if(advScheme == defValue || gammaVal == defValue ||
	 NsubVal == defValue || deltaVal == defValue) {
         cout << "ERROR: advScheme or gamma " << endl;
         cout << "or Nsub or delta0 is " << endl;
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
      Sthresh = Nthresh*Tthresh/pow(Nthresh,gamma0-1);
      
      delta0 = deltaVal.asDouble();
      if(procID==0) cout << "relaxation constant = " << delta0 << endl;
      if(delta0 >= 1.0) {
         printf("ERROR: delta0>=1 ==> cvac<=V \n");
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
      if(procID==0) setXminBoundary(N, 0.0, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(N, 0.0, 1.0); 
      setZboundaryPeriodic(N);  
      Xgrid.communicate(N);
      Nold  = N;
   } else {
      cout << "value for Physics variable \"N\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
 
   const Json::Value Vxvar = Phys.get("Vx",defValue);
   if(Vxvar.isObject()) { 
      Xgrid.setInitialProfile(Vx,Vxvar);
      if(procID==0) setXminBoundary(Vx, 0.0, -1.0);   
      if(procID==numProcs-1) setXmaxBoundary(Vx, 0.0, -1.0);   
      setZboundaryPeriodic(Vx);  
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
      if(procID==0) setXminBoundary(Vz, 0.0, -1.0);   
      if(procID==numProcs-1) setXmaxBoundary(Vz, 0.0, -1.0);   
      setZboundaryPeriodic(Vz);  
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
      if(procID==numProcs-1) setXmaxBoundary(P, 0.0, 1.0);   
      setZboundaryPeriodic(P);  
      Xgrid.communicate(P);
      S = P/pow(N,gamma0-1.0);
      Sold  = S;
   } else {
      cout << "value for Physics variable \"P\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   const Json::Value Byvar = Phys.get("By",defValue);
   if(Byvar.isObject()) { 
      Xgrid.setInitialProfile(By,Byvar);
      if(procID==0) setXminBoundary(By, B0, 0.0);   
      if(procID==numProcs-1) setXmaxBoundary(By, 0.0, 0.0);
      setZboundaryPeriodic(By);  
      Xgrid.communicate(By);
      Byold = By;
   } else {
      cout << "value for Physics variable \"By\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   T = P/N;
   if(min(P)<0.0) cout << " T IS LESS THAN ZERO " << endl;
   //Cs = sqrt(gamma0*P/N);
   Cs = pow(gamma0*P/N,0.5);
   Xgrid.DDX(Jz,By); 
   Xgrid.DDX(Jzcc,By); 
   Xgrid.communicate(Jz);
   Xgrid.communicate(Jzcc);
   Jz0 = Jz;
   if(min(T)<0.0) cout << " T IS LESS THAN ZERO " << endl;
   //eta = eta0/T/sqrt(T);
   eta = eta0/T/pow(T,0.5);
   Xgrid.InterpToCellEdges(VxBy_x,Vx*By,By,"C2",0);
   Xgrid.InterpToCellEdges(eta_x,eta,eta,"C2",0);
   Ez = eta_x*Jz - VxBy_x; 
   Ezhalf = Ez;
   //
   Xgrid.InterpToCellEdges(VzBy_z,Vz*By,By,"C2",1);
   Xgrid.InterpToCellEdges(eta_z,eta,eta,"C2",1);
   Ex = eta_z*Jx + VzBy_z; 
   Exhalf = Ex;


   // need inital flux calculation before add to output
   // and before first call to Physics.advance   
   computeFluxes(Xgrid, 1); // 1 for first order calculation   
  

   // add stuff to output files
   //
   dataFile.add(N, "N", 1);  // density 
   dataFile.add(Mx, "Mx", 1);  // momentum density 
   dataFile.add(Mz, "Mz", 1);  // momentum density 
   dataFile.add(S, "S", 1);  // entropy density
   dataFile.add(By, "By", 1);  // magnetic field
   dataFile.add(P, "P", 1);  // pressure
   dataFile.add(T, "T", 1);  // temperature
   dataFile.add(eta, "eta", 1);  // resistivity
   dataFile.add(Vx, "Vx", 1);  // x-velocity
   dataFile.add(Vz, "Vz", 1);  // z-velocity
   dataFile.add(Jz, "Jz", 1);  // current density
   dataFile.add(Jz0, "Jz0", 1);  // curl of B
   dataFile.add(Ez, "Ez", 1);  // z-electric field
   dataFile.add(Jx, "Jx", 1);  // current density
   dataFile.add(Jx0, "Jx0", 1);  // curl of B
   dataFile.add(Ex, "Ex", 1);  // x-electric field
   dataFile.add(Cs,"Cs",1);  // sound speed
   dataFile.add(gamma0,"gamma0",0); 
   //
   dataFile.add(FluxN_x, "FluxN_x", 1);  
   dataFile.add(FluxMx_x, "FluxMx_x", 1);  
   dataFile.add(FluxMz_x, "FluxMz_x", 1);  
   dataFile.add(FluxS_x, "FluxS_x", 1);  
   dataFile.add(FluxBy_x, "FluxBy_x", 1);  
   dataFile.add(FluxEz_x, "FluxEz_x", 1);  
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
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   //domainGrid* mesh = domainGrid::mesh;

   // Explicit forward advance from n to n+1 using subcycling in time 
   //  
   double thisdt, expFact, Jsq;
   for (auto n=1; n<Nsub+1; n++) {
      thisdt = dt*n/Nsub;

      // Update N, M, S, and B terms from n to n+1
      for (auto i=nXg; i<nXcc-nXg; i++) {
         for (auto j=nZg; j<nZcc-nZg; j++) {
	 
         Jsq = Jzcc(i,j)*Jzcc(i,j)+Jxcc(i,j)*Jxcc(i,j);

	 N(i,j)  = Nold(i,j)  - thisdt*(FluxN_x(i,j)  - FluxN_x(i-1,j))/Xgrid.dX
	                      - thisdt*(FluxN_z(i,j)  - FluxN_z(i,j-1))/Xgrid.dZ;
         Mx(i,j) = Mxold(i,j) - thisdt*(FluxMx_x(i,j) - FluxMx_x(i-1,j))/Xgrid.dX
                              - thisdt*(FluxMx_z(i,j) - FluxMx_z(i,j-1))/Xgrid.dZ
		              - thisdt*Jzcc(i,j)*Byold(i,j);
         Mz(i,j) = Mzold(i,j) - thisdt*(FluxMz_x(i,j) - FluxMz_x(i-1,j))/Xgrid.dX
                              - thisdt*(FluxMz_z(i,j) - FluxMz_z(i,j-1))/Xgrid.dZ
			      + thisdt*Jxcc(i,j)*Byold(i,j);
         S(i,j)  = Sold(i,j)  - thisdt*(FluxS_x(i,j)  - FluxS_x(i-1,j))/Xgrid.dX
                              - thisdt*(FluxS_z(i,j)  - FluxS_z(i,j-1))/Xgrid.dZ
		 + thisdt*(gamma0-1.0)*eta(i,j)*pow(Nold(i,j),1.0-gamma0)*Jsq;
         
	 By(i,j) = Byold(i,j) - thisdt*(FluxBy_x(i,j)-FluxBy_x(i-1,j))/Xgrid.dX
	                      - thisdt*(FluxBy_z(i,j)-FluxBy_z(i,j-1))/Xgrid.dZ;
	 //By(i,j) = Byold(i,j) + thisdt*(Ezold(i,j)-Ezold(i-1,j))/Xgrid.dX;
	 
	 if(N(i,j)<=Nthresh) N(i,j) = Nthresh;
	 //if(M(i,j)<0.0) M(i,j) = 0.0;
	 if(S(i,j)<=Sthresh) S(i,j) = Sthresh;

         }
      }

      timeDomain* tmesh = timeDomain::tmesh;
      //cout << "thist = " << tmesh->tSim << endl;
      double thist = tmesh->tSim;
      B0 = thist*40.0;
      if(B0>4.0) B0 = 4.0;
      
      if(procID==0) {
         setXminExtrap(N,0);   
         setXminExtrap(Mx,0);   
         setXminExtrap(Mz,0);   
         setXminExtrap(S,0);
         //setXminExtrap(N,2);
         //setXminExtrap(Mx,2);
         //setXminExtrap(S,2);
         setXminBoundary(By, B0, 0.0);   

	 /*
	 if(N.at(1)<Nthresh) {
            N.at(0) = Nthresh;
            N.at(1) = Nthresh;
	 }
	 if(N.at(0)<Nthresh) {
            N.at(0) = Nthresh;
	 }
	 if(S.at(1)<Sthresh) {
            S.at(0) = Sthresh;
            S.at(1) = Sthresh;
	 }
	 if(S.at(0)<Sthresh) {
            S.at(0) = Sthresh;
	 }
	 
	 if(M.at(1)<0.0) {
            M.at(0) = M.at(2);
            M.at(1) = M.at(2);
	 }
	 if(M.at(0)<0.0) {
            M.at(0) = M.at(1);
	 }
	 */
	 
      }
     
      if(procID==numProcs-1) {
         setXmaxBoundary(N, 0.0, 1.0);   
         setXmaxBoundary(Mx, 0.0, 1.0);   
         setXmaxBoundary(Mz, 0.0, 1.0);   
         setXmaxBoundary(S, 0.0, 1.0);   
         setXmaxBoundary(By, 0.0, 0.0);   
      }
      setZboundaryPeriodic(N);
      setZboundaryPeriodic(Mx);
      setZboundaryPeriodic(Mz);
      setZboundaryPeriodic(S);
      setZboundaryPeriodic(By);
      Xgrid.communicate(N);
      Xgrid.communicate(Mx);
      Xgrid.communicate(Mz);
      Xgrid.communicate(S);
      Xgrid.communicate(By);

      // compute fluxes using partialy updated fields at n+1/2 for 2nd order
      // or n+1 for 1st order
      computeFluxes(Xgrid, Nsub); // compute second order fluxes at n+1/2

   } // finish subcycle steps for N, M, S, and B


   // Now update electric field (stag in time wrt others)
   //
   //Xgrid.InterpToCellEdges(VxBy_x,M/N*By,By,"C2");
   for (auto i=nXg-1; i<nXcc-nXg; i++) {
      for (auto j=nZg; j<nZcc-nZg; j++) {

	 // solve for z-electric field
	 //
         if(thisdt/delta0/eta_x(i,j)<=1.0e-3) {
            expFact = 1.0-thisdt/delta0/eta_x(i,j)+0.5*pow(thisdt/delta0/eta_x(i,j),2);
         }
         else {
            expFact = exp(-thisdt/delta0/eta_x(i,j));
         }
         Ez(i,j) = Ezold(i,j)*expFact 
    	         - ( VxBy_x(i,j) - eta_x(i,j)*(By(i+1,j)-By(i,j))/Xgrid.dX 
	           )*(1.0-expFact);

	 // solve for x-electric field
	 //
	 //Ex(i,j) = Exold(i,j)*0.0;
         if(thisdt/delta0/eta_z(i,j)<=1.0e-3) {
            expFact = 1.0-thisdt/delta0/eta_z(i,j)+0.5*pow(thisdt/delta0/eta_z(i,j),2);
         }
         else {
            expFact = exp(-thisdt/delta0/eta_z(i,j));
         }
         Ex(i,j) = Exold(i,j)*expFact 
    	         + ( VzBy_z(i,j) - eta_z(i,j)*(By(i,j+1)-By(i,j))/Xgrid.dZ 
	           )*(1.0-expFact);


      }
   }
   setZboundaryPeriodic(Ez);
   setZboundaryPeriodic(Ex);
   Xgrid.communicate(Ez);
   Xgrid.communicate(Ex);

   // update old fields
   Nold = N;
   Mxold = Mx;
   Mzold = Mz;
   Sold = S;
   Byold = By;
   Ezhalf = (Ez+Ezold)/2.0;
   Exhalf = (Ex+Exold)/2.0;
   Ezold = Ez;
   Exold = Ex;
   
   // compute fluxes using fully updated fields at n+1
   computeFluxes(Xgrid, 1);

} // end Physics.advance


void computeFluxes(const domainGrid& Xgrid, const int order)
{

   const int nXcc = Xgrid.nXcc;
   const int nXce = Xgrid.nXce;
   const int nZcc = Xgrid.nZcc;
   const int nZce = Xgrid.nZce;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   matrix2D<double> FluxNcc_x, FluxMxcc_x, FluxMzcc_x, FluxScc_x;
   matrix2D<double> FluxNcc_z, FluxMxcc_z, FluxMzcc_z, FluxScc_z;
   matrix2D<double> Cspeed, Cspeed2, Cvaceff, FluxB0cc, Exprime, Ezprime; 
 
   Cspeed.initialize(nXcc,nZcc,0.0);
   Cspeed2.initialize(nXcc,nZcc,0.0);
   Cvaceff.initialize(nXcc,nZcc,0.0);
   FluxNcc_x.initialize(nXcc,nZcc,0.0);
   FluxMxcc_x.initialize(nXcc,nZcc,0.0);
   FluxMzcc_x.initialize(nXcc,nZcc,0.0);
   FluxScc_x.initialize(nXcc,nZcc,0.0);
   FluxNcc_z.initialize(nXcc,nZcc,0.0);
   FluxMxcc_z.initialize(nXcc,nZcc,0.0);
   FluxMzcc_z.initialize(nXcc,nZcc,0.0);
   FluxScc_z.initialize(nXcc,nZcc,0.0);
   FluxB0cc.initialize(nXcc,nZcc,0.0);
   Exprime.initialize(nXcc,nZce,0.0);
   Ezprime.initialize(nXce,nZcc,0.0);

   //  define derived variables
   //
   Vx = Mx/N;
   Vz = Mz/N;
   P  = S*pow(N,gamma0-1.0);
   T  = P/N;
   Cs = sqrt(gamma0*P/N);
   eta = eta0/T/sqrt(T); //*(1.0+10*Nthresh/N);
   Xgrid.InterpToCellEdges(eta_x,eta,eta,"C2",0);
   Xgrid.communicate(eta_x);
   Xgrid.InterpToCellEdges(eta_z,eta,eta,"C2",1);
   Xgrid.communicate(eta_z);

   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeed  = abs(Vx) + Cs; // adv flux jacobian (may want to include B^2/N in pressure)
   Cspeed2 = abs(Vz) + Cs; // adv flux jacobian
   for (auto i=0; i<nXcc; i++) {
      for (auto j=0; j<nZcc; j++) {
         if(Cspeed2(i,j)>Cspeed(i,j)) Cspeed(i,j) = Cspeed2(i,j);
      }
   }
   Cvaceff.initialize(nXcc,nZcc,1.0/sqrt(delta0));
   FluxNcc_x = Mx;
   FluxNcc_z = Mz;
   FluxMxcc_x = Mx*Vx + P;
   FluxMxcc_z = Mx*Vz;
   FluxMzcc_x = Mz*Vx;
   FluxMzcc_z = Mz*Vz + P;
   FluxScc_x = Vx*S;
   FluxScc_z = Vz*S;
   FluxB0cc = Vx*By;
   FluxEz_x = -By;
   FluxEx_z = By;

   // compute advective flux using
   // specified scheme from input file
   //
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxNcc_x,Cspeed,N,0,Nsub);
      Xgrid.computeFluxTVD(FluxMx_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMxcc_x,Cspeed,Mx,0,Nsub);
      Xgrid.computeFluxTVD(FluxMz_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMzcc_x,Cspeed,Mz,0,Nsub);
      Xgrid.computeFluxTVD(FluxS_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxScc_x,Cspeed,S,0,Nsub);
      //Xgrid.computeFluxTVD(VxBy_x,FluxL,FluxR,FluxRatio,FluxLim,
      //                     FluxB0cc,Cspeed,B,0,Nsub);
      
      Xgrid.computeFluxTVD(FluxN_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxNcc_z,Cspeed,N,1,Nsub);
      Xgrid.computeFluxTVD(FluxMx_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxMxcc_z,Cspeed,Mx,1,Nsub);
      Xgrid.computeFluxTVD(FluxMz_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxMzcc_z,Cspeed,Mz,1,Nsub);
      Xgrid.computeFluxTVD(FluxS_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxScc_z,Cspeed,S,1,Nsub);
   }
   else {
      Xgrid.InterpToCellEdges(FluxN_x,FluxNcc_x,Vx,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxMx_x,FluxMxcc_x,Vx,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxMz_x,FluxMzcc_x,Vz,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxS_x,FluxScc_x,Vx,advScheme0,0);
      //Xgrid.InterpToCellEdges(FluxB,FluxBcc,Vx,advScheme0,0);
      //Xgrid.InterpToCellEdges(FluxEz_x,FluxEzcc,Ez,advScheme0,0);
      //Xgrid.InterpToCellEdges(FluxEz,FluxEzcc,Ez,"C2",0);
   } 
   Xgrid.InterpToCellEdges(VxBy_x,FluxB0cc,By,"C2",0);
   Xgrid.InterpToCellEdges(VzBy_z,Vz*By,By,"C2",1);

   if(procID==0) {
      setXminBoundary(FluxN_x, 0.0, 0.0);   
      //setXminBoundary(FluxMx, (P.at(2)+P.at(1))/2.0, 0.0);   
      //setXminBoundary(FluxMx, (FluxMxcc.at(2)+FluxMxcc.at(1))/2.0, 0.0);   
      setXminBoundary(FluxS_x, 0.0, 0.0);
      //setXminBoundary(VxBy_x, (FluxB0cc.at(2)+FluxB0cc.at(1))/2.0, 0.0);
   }   
   //FluxBy_x = FluxBy-Ezprime;
   FluxBy_x = -Ez;
   FluxBy_z = Ex;

   Xgrid.communicate(VxBy_x);
   Ezprime = Ez + VxBy_x;
   Jz = (Ezhalf + VxBy_x)/eta_x;
   Xgrid.communicate(Jz);
   Xgrid.InterpToCellCenter(Jzcc,Jz);
   Xgrid.communicate(Jzcc);
   Xgrid.DDX(Jz0,By); 
   Xgrid.communicate(Jz0);

   Xgrid.communicate(VzBy_z);
   Exprime = Ex - VzBy_z;
   Jx = (Exhalf - VzBy_z)/eta_z;
   Xgrid.communicate(Jz);
   Xgrid.InterpToCellCenter(Jxcc,Jx);
   Xgrid.communicate(Jxcc);
   Xgrid.DDZ(Jx0,By); 
   Jx0 *= -1.0; 
   Xgrid.communicate(Jz0);
   

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

void setZboundaryPeriodic(matrix2D<double>& var)
{
   
   const int thisnX = var.size0();
   const int thisnZ = var.size1();

   domainGrid* mesh = domainGrid::mesh;
   const int jshift = mesh->nZg;
 
   assert(jshift==2 || jshift==1);

   if(jshift==2) {
      for (auto i=0; i<thisnX; i++) {
         var(i,1) = var(i,thisnZ-3);
         var(i,0) = var(i,thisnZ-4);

         var(i,thisnZ-2) = var(i,2);
         var(i,thisnZ-1) = var(i,3);
      }
   } 
   else { 
      for (auto i=0; i<thisnX; i++) {
         var(i,0) = var(i,thisnZ-2);
         var(i,thisnZ-1) = var(i,1);
      }
   }



}


void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   matrix2D<double> Cchar(Cs);
   double Cmax;
   Cchar += sqrt(Vx*Vx);
   Cmax = max(Cchar);
   //cout << "Cmax = " << Cmax << endl;

   const double dX = Xgrid.dX;
   double dtCFL_sound = dX/Cmax;
   double dtCFL_light = dX*sqrt(delta0);
   double dtmax = min(dtCFL_sound,dtCFL_light);
   //double dtmax = dtCFL_sound;
   
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);


   //dtSim = 1.0e-4;
   if(procID==0) {
      cout << "sigma_0*dt/delta = " << dtSim/delta0/eta0 << endl;
      cout << "dtSim = " << dtSim << endl;
   }
}

