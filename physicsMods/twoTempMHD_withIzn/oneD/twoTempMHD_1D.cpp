/***
 * 
 * physics module for 1D two Temp MHD dpf rundown
 * (1D railgun) with ionization physics included
 *
 * Includes Hall field (Ex) in direction of flow
 * assuming force balance in that direction.
 * Because 1D, there is no dEx/dz and so Hall field
 * only enters in how the ions gain energy
 *
 *
 * Input Scales:
 * Xscale (length [m])
 * Nscale (density [1/m^3])
 * Amass  (atomic mass)
 * Iscale (current [Amps])
 *
 * Derived Scales:
 * Pscale   = mu0*I^2/(8*pi^2*Xscale^2) [J/m^3]
 * rhoScale = Mi*Nscale [kg/m^3], Mi=Amass*amu [kg]
 * Vscale   = sqrt(Pscale/rhoScale) [m/s]
 * tscale   = Xscale/Vscale [s]
 * Bscale   = sqrt(mu0*Pscale) [tesla]
 * Jscale   = Bscale/Xscale/mu0 [Amps/m^2]
 * Escale   = Vscale*Bscale [V/m]
 *
 * Governing Equations:
 * dN/dt  + d(Mx)/dx = 0
 * dMx/dt + d(Mx*Vx + P)/dx = -Jz*By
 * dEi/dt + d((Ei+Pi)*Vx)/dx = NUdotE + Qie 
 * dEe/dt + d((Ee+Pe)*Vx)/dx = Jz*Ez  - Qie - NUdotE + SEe
 * dNe/dt + d(Zbar*Mx)/dx = Se
 * dBy/dt + d(-Ez)/dx = 0
 * delta0*dEz/dt   = Jz0 - Jz
 * Le0sq*dJz/dt = Ne*([Ez + Vx*By] - eta*Jz)
 *
 * Derived Parameters: 
 * Zbar = Ne/N
 * Nn = N-Ne
 * Vx = Mx/N;
 * Jz0 = curl(By)
 * Pi = (Ei - 0.5*N*Vx^2)*(gamma0-1)
 * Pe = Ee*(gamma0-1)
 * P = Pi + Pe
 * Ti = Pi/N
 * Te = Pe/Ne
 * eta = eta0/Te^1.5
 * NUdotE = -Ux*(Jz*By + dPe/dx) (Ex from force balance)
 * Qie = 3.0*me/Mi*Ne/taue*(Te-Ti)
 * taue = taue0*Te^1.5/Ne
 * Se = Ne*Nn*kiz
 * SEe = -Uizn*Se
 *
 * Dimensionless parameters:
 * delta0 = (V0/cvac)^2
 * Le0sq = (Le0/r0)^2, Le0 = cvac/wpe0
 * gamma0 = 2/degFreedom + 1
 *
 * Dimensional Parameters:
 * eta0 = 1.03/10/Te0^1.5/(r0^2*mu0/t0)
 * taue0 = 3.44e5/10*Te0^1.5/N0/t0
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
#include <cmath>

#include "json/json.h"
#include "vectorMath.h"
#include "domainGrid.h"
#include "timeDomain.h"
#include "HDF5dataFile.h"

#include "Physics.h"


using namespace std;

string advScheme0;  // advection differencing scheme
string XlowBC, XhiBC;  // boundary condition strings
string geometry0;   // CAR or CYL
double gamma0;      // adiabatic coefficient
double eta0;        // resistivity coefficient
double taue0;       // ele collision time coefficient
double taui0;       // ion collision time coefficient
double mM;          // me/Mi
double etaVis0;     // numerical viscosity coefficient
double delta0;      // relaxation const (V0/cvac)^2
double Le0sq;    // normalized electron skin depth squared
double B0, B00;     // boundary value of magnetic field
int Nsub;           // time-solver subcycle steps
vector<double> N, Mx, Ei, Ee, By, Ez, Jz;   // time-evolving variables
vector<double> eta, Cs, Va, Cspeed, Vx, P, Pe, Pi, Te, Ti, Jz0, Jz0stagTime, Qvisc; // derived variables
vector<double> etace, Jzcc, Ezcc, Eidealz, Nece;
vector<double> Nold, Mxold, Eiold, Eeold, Byold, Ezold, Jzold;
vector<double> FluxRatio, FluxLim;
vector<double> FluxR, FluxL;  // flux at cell-edges   
vector<double> FluxN, FluxMx, FluxEi, FluxEe;
vector<double> Qie, NUdotE, JdotE, taue, nue_spi, nue_vac;

// new stuff for ionization
//
const double Uizn = 13.6;
const double g0 = 2, g1=1;
double Zmin = 1.0;
vector<double> Zbar, Neold, Se, SEe, Ne, Nn, nue_neu, nue_izn;
vector<double> FluxNe;

double Nscale, Xscale, Amass, Iscale, dtIscale;
double Pscale, Tscale, Bscale, Jscale, Ezscale, Vscale, tscale, Mi;
double epsilonRel, meRel;

vector<double> hy_cc, hy_ce; // y-dimension lame coefficient

// Set lower threshold values for N, T, and S
double Nthresh=1.0e-4, Tthresh=2.5e-2, Ethresh, Pthresh;

// Set conditions for vaccum resistivity model
double NvacC=0.01, NvacP = 4;

void initializeMemberVars(const domainGrid& );
void computeFluxes(const domainGrid&, const int);
void computeJ0(vector<double>&, const vector<double>&, const domainGrid& );
void updatePhysicalVars(const domainGrid& );
void updateCollisionTerms(const domainGrid& );
void computeIdealEatEdges( const domainGrid&, const vector<double>& );
void computeJouleHeatingTerms(const domainGrid& );
void setNeutralInteractionRates(const domainGrid&,    const vector<double>&,
                                const vector<double>&, const vector<double>& );
void advanceSemiImplicitVars(const domainGrid&, const int, const double, const double);

void parseInputFile(const domainGrid&, const Json::Value& );
void addMembersToDataFile( HDF5dataFile& );

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
   const int nXg  = Xgrid.nXg;
   const int nXcc = Xgrid.Xcc.size();
   const int nXce = Xgrid.Xce2.size();
   
   initializeMemberVars(Xgrid);

   parseInputFile(Xgrid,root);

   // define remaining vars from those specified in input file
   //
   Ne = Zbar*N;
   Nn = N-Ne; 
   Mx = N*Vx;
   Pe = Ne*Te;
   Pi = N*Ti;
   Ei = 0.5*Mx*Mx/N + Pi/(gamma0-1.0);
   Ee = Pe/(gamma0-1.0);

   // initialize current density from magnetic field
   //
   computeJ0(Jz0,By,Xgrid);
   Jz = Jz0;
   Jz0stagTime = Jz;

   Xgrid.DDX(Jzcc,hy_cc*By);
   Jzcc = Jzcc/hy_cc; 
   Xgrid.communicate(Jzcc);

   // need a complete RHS eval before writing at t=0
   // and before first call to Physics.advance   

   updatePhysicalVars(Xgrid);
   
   updateCollisionTerms(Xgrid);  
      
   computeIdealEatEdges(Xgrid,Cspeed);
   
   Ez = etace*Jz + Eidealz; // initialize electric field from Ohm's law 

   computeJouleHeatingTerms(Xgrid);
   
   computeFluxes(Xgrid, 1); // 1 for first order calculation   

   // dont forget to initialize all the old values
   //
   Nold  = N;
   Mxold = Mx;  
   Eiold = Ei;
   Eeold = Ee;
   Byold = By;
   Neold = Ne;
   Ezold = Ez;
   Jzold = Jz;

   // add selected members to output files
   //
   addMembersToDataFile( dataFile );

} // end Physics.initilize


void Physics::advance(const domainGrid& Xgrid, const double a_dt)
{
   const int nXcc = Xgrid.nXcc;
   const int nXg = Xgrid.nXg;
   const double dX = Xgrid.dX;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   //domainGrid* mesh = domainGrid::mesh;

   // Explicit forward advance from n to n+1 using subcycling in time 
   //  
   double thisdt, thist;
   double dt_stag;
   timeDomain* tmesh = timeDomain::tmesh;
   for (auto n=1; n<Nsub+1; n++) {
      
      thisdt = a_dt*n/Nsub;
      //thist = tmesh->tSim - thisdt*(2-n);
      //dt_stag = a_dt; 
      //if(thist==thisdt && n==1) dt_stag = thisdt;

      for (auto i=nXg; i<nXcc-nXg; i++) {
	 N.at(i)  = Nold.at(i)  - thisdt*(FluxN.at(i+1) -FluxN.at(i) )/hy_cc.at(i)/dX;
         Mx.at(i) = Mxold.at(i) - thisdt*(FluxMx.at(i+1)-FluxMx.at(i))/hy_cc.at(i)/dX
		  - thisdt*Jzcc.at(i)*By.at(i);
	 if(geometry0=="CYL") Mx.at(i) = Mx.at(i) + thisdt*P.at(i)/hy_cc.at(i);
         Ei.at(i) = Eiold.at(i) - thisdt*(FluxEi.at(i+1)-FluxEi.at(i))/hy_cc.at(i)/dX
      		  + thisdt*(NUdotE.at(i) + Qie.at(i));
         Ee.at(i) = Eeold.at(i) - thisdt*(FluxEe.at(i+1)-FluxEe.at(i))/hy_cc.at(i)/dX
      		  + thisdt*(JdotE.at(i) - Qie.at(i) - NUdotE.at(i) + SEe.at(i));
	 By.at(i) = Byold.at(i) + thisdt*(Ez.at(i+1)-Ez.at(i))/dX;
	 
	 if(N.at(i)<=Nthresh) N.at(i) = Nthresh;
	 if(N.at(i)!=N.at(i)) cout << "bout to go bad: N.at(i) = " << N.at(i) << endl;
	 Ethresh = Pthresh/(gamma0-1.0) + 0.5*Mx.at(i)*Mx.at(i)/N.at(i);
	 if(Ei.at(i)<=Ethresh) Ei.at(i) =Ethresh;
	 Ethresh = Pthresh/(gamma0-1.0);
	 if(Ee.at(i)<=Ethresh) Ee.at(i) =Ethresh;
	 
         Ne.at(i) = Neold.at(i) - thisdt*(FluxNe.at(i+1)-FluxNe.at(i))/hy_cc.at(i)/dX
                  + thisdt*Se.at(i);
	 if(Ne.at(i)<N.at(i)*Zmin) Ne.at(i) = N.at(i)*Zmin;
	 if(Ne.at(i)>N.at(i))      Ne.at(i) = N.at(i);

      }
      thist = tmesh->tSim + thisdt;
      dt_stag = a_dt; 
      if(thist==thisdt && n==1) dt_stag = thisdt;

      //  apply boundary conditions and communicate
      //
      B0 = thist/dtIscale*B00;
      if(B0>B00) B0 = B00;
      //if(thist>=0.03) B0 = B00*pow(0.03/thist,2);

      if(procID==0) {
         Xgrid.setXminBoundary(N,  0.0, 1.0);   
         Xgrid.setXminBoundary(Ne, 0.0, 1.0);   
         Xgrid.setXminBoundary(Mx, 0.0, -1.0);   
         Xgrid.setXminBoundary(Ei, 0.0, 1.0);
         Xgrid.setXminBoundary(Ee, 0.0, 1.0);
	 if(N.at(0)<Nthresh || N.at(1)<Nthresh) {
            N.at(0) = Nthresh;
            N.at(1) = Nthresh;
	 }
         //Xgrid.setXminBoundary(By, 0.0, -1.0);   
         
         if(XlowBC.compare("axis") == 0 || XlowBC.compare("symmetry") == 0) {
            Xgrid.setXminBoundary(By, 0.0, -1.0);   
         } 
         else if(XlowBC.compare("insulator") == 0) {
            Xgrid.setXminBoundary_J(By,B0,0.0);   
         } 
         else {
            cout << "low MagFieldBC not defined !!!! " << endl;
         }
         
      }
     
      if(procID==numProcs-1) {
         Xgrid.setXmaxBoundary(N,  0.0, 1.0);   
         Xgrid.setXmaxBoundary(Ne, 0.0, 1.0);   
         Xgrid.setXmaxBoundary(Mx, 0.0, -1.0);   
         Xgrid.setXmaxBoundary(Ei, 0.0, 1.0);   
         Xgrid.setXmaxBoundary(Ee, 0.0, 1.0); 
         
         if(XhiBC.compare("conductor") ==0 ) {
            Xgrid.setXmaxBoundary_J(By,0,1);   
         } 
         else if(XhiBC.compare("insulator") == 0) {
            //cout << "hy_cc.at(nXcc-nXg) = " << hy_cc.at(nXcc-nXg) << endl;	 
            Xgrid.setXmaxBoundary_J(By,B0,0); 
         } 
         else {
            cout << "Xhi MagFieldBC not defined !!!! " << endl;
         }
         
      }

      Xgrid.communicate(N);
      Xgrid.communicate(Ne);
      Xgrid.communicate(Mx);
      Xgrid.communicate(Ei);
      Xgrid.communicate(Ee);
      Xgrid.communicate(By);

      // update values needed for SemiImplicit advance
      // 
      updatePhysicalVars(Xgrid);
      
      updateCollisionTerms(Xgrid);  
   
      //computeIdealEatEdges(Xgrid,abs(Vx));
      computeIdealEatEdges(Xgrid,Cspeed);
      
      advanceSemiImplicitVars(Xgrid,n,dt_stag,thist);
   
      // update Joule heating source terms for energy equations
      //
      computeJouleHeatingTerms(Xgrid);
      
      // update fluxes
      //
      computeFluxes(Xgrid, Nsub);

   } // finish subcycle steps

   // update old values for explicit vars
   //
   Nold  = N;
   Mxold = Mx;
   Eiold = Ei;
   Eeold = Ee;
   Byold = By;
   Neold = Ne;

} // end Physics.advance

void advanceSemiImplicitVars(const domainGrid& Xgrid, const int stage, 
                             const double thisdt, const double thist)
{
   const int nXce = Xgrid.nXce2;
   const int nXg = Xgrid.nXg;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   //if(!procID) cout << "dt_stag = " << thisdt << endl;
   //if(!procID) cout << "thist   = " << thist << endl;
   //if(!procID) cout << "stage   = " << stage << endl;
   
   // only compute Jz0 at second stage for electric field advance 
   // (leap-frog scheme for waves => du-fort and frenkel scheme for diffusion)
   //
   if(stage==2) {
      computeJ0(Jz0stagTime,By,Xgrid);
   }

   // advance electric field and current density
   //
   double d0, e0;
   for (auto i=nXg; i<nXce-nXg; i++) {
      d0 = delta0/thisdt;	   
      e0 = meRel*Le0sq/thisdt/Nece.at(i);	   
      //
      Ez.at(i) = Jz0stagTime.at(i) + d0*Ezold.at(i)
               - (e0*Jzold.at(i) - Eidealz.at(i))/(e0 + etace.at(i));
      Ez.at(i) /= d0 + 1.0/(e0 + etace.at(i)); 
      //
      Jz.at(i) = e0*Jzold.at(i) + Ez.at(i) - Eidealz.at(i);
      Jz.at(i) /= e0 + etace.at(i);
    
   }
   //cout << "Xce(nXg-1) =" << Xgrid.Xce.at(nXg-1) << endl;
   //cout << "Xce(nXcc-Xg-1) =" << Xgrid.Xce.at(nXcc-nXg-1) << endl;
   //if(procID==numProcs-1) {
   //   Ez.at(nXcc-nXg) = Ez.at(nXcc-nXg-1);
   //}   
   Xgrid.communicate(Ez);
   Xgrid.communicate(Jz);

   if(stage==1) {
      Ezold = Ez;
      Jzold = Jz;
   }

}


void computeFluxes(const domainGrid& Xgrid, const int order)
{ // order not being used right now

   const int nXcc = Xgrid.Xcc.size();
   const int nXg = Xgrid.nXg;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   vector<double> FluxNcc, FluxMxcc, FluxEicc, FluxEecc;
   vector<double> FluxNecc, FluxVisc, dVdx; 
   FluxNcc.assign(nXcc,0.0);
   FluxMxcc.assign(nXcc,0.0);
   FluxEicc.assign(nXcc,0.0);
   FluxEecc.assign(nXcc,0.0);
   FluxNecc.assign(nXcc,0.0);
   FluxVisc.assign(nXcc,0.0);
   dVdx.assign(nXcc,0.0);


   // compute fluxes at cell center
   //
   FluxNcc = hy_cc*Mx;
   FluxMxcc = hy_cc*(Mx*Vx + P);
   FluxEicc = hy_cc*(0.5*Vx*Mx + Pi*gamma0/(gamma0-1.0) )*Vx;
   FluxEecc = hy_cc*Pe*gamma0/(gamma0-1.0)*Vx;
   FluxNecc = Zbar*FluxNcc;
   
   // compute viscous terms
   //
   vector<double> etaVisc;
   etaVisc.assign(nXcc,0.0); 
   etaVisc = N*etaVis0;
   //Xgrid.DDX(FluxVisc,-4.0*etaVisc/3.0*Vx);
   Xgrid.DDX(dVdx,Vx);
   //Xgrid.communicate(FluxVisc);
   Xgrid.communicate(dVdx);
   Qvisc = 4.0/3.0*etaVisc*dVdx*dVdx;
   FluxVisc = -4.0/3.0*etaVisc*dVdx;
   FluxMxcc = FluxMxcc + FluxVisc;

   // compute advective flux using
   // specified scheme from input file
   //
   //Nsub = 2;
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxNcc,Cspeed,hy_cc*N,"minmod",Nsub);
      Xgrid.computeFluxTVD(FluxMx,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxMxcc,Cspeed,hy_cc*Mx,"minmod",Nsub);
      Xgrid.computeFluxTVD(FluxEi,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxEicc,Cspeed,hy_cc*Ei,"minmod",Nsub);
      Xgrid.computeFluxTVD(FluxEe,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxEecc,Cspeed,hy_cc*Ee,"minmod",Nsub);
      Xgrid.computeFluxTVD(FluxNe,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxNecc,Cspeed,hy_cc*Ne,"minmod",Nsub);
   }
   else {
      Xgrid.computeFluxTVDsimple(FluxN,FluxL,FluxR,FluxNcc,Cspeed,hy_cc*N, advScheme0);
      Xgrid.computeFluxTVDsimple(FluxMx,FluxL,FluxR,FluxMxcc,Cspeed,hy_cc*Mx,advScheme0);
      Xgrid.computeFluxTVDsimple(FluxEi,FluxL,FluxR,FluxEicc,Cspeed,hy_cc*Ei,advScheme0);
      Xgrid.computeFluxTVDsimple(FluxEe,FluxL,FluxR,FluxEecc,Cspeed,hy_cc*Ee,advScheme0);
      Xgrid.computeFluxTVDsimple(FluxNe,FluxL,FluxR,FluxNecc,Cspeed,hy_cc*Ne,advScheme0);
   } 
   //FluxMx = FluxMx + FluxVisc;

   if(procID==0) {
      Xgrid.setXminFluxBC(FluxN, 0.0, 0.0);
      //cout << "hy_ce.at(1)" << hy_ce.at(1) << endl;   
      Xgrid.setXminFluxBC(FluxMx, hy_ce.at(nXg)*(P.at(nXg)+P.at(nXg-1))/2.0, 0.0);   
      Xgrid.setXminFluxBC(FluxEi, 0.0, 0.0);
      Xgrid.setXminFluxBC(FluxEe, 0.0, 0.0);
      Xgrid.setXminFluxBC(FluxNe, 0.0, 0.0);
   }
   if(procID==numProcs-1) {
      Xgrid.setXmaxFluxBC(FluxN, 0.0, 0.0);
      Xgrid.setXmaxFluxBC(FluxEi, 0.0, 0.0);
      Xgrid.setXmaxFluxBC(FluxEe, 0.0, 0.0);
      const int thisnX = P.size();
      //cout << "hy_ce.at(thisnX-3)" << hy_ce.at(thisnX-3) << endl;   
      double P0 = hy_ce.at(thisnX-nXg)*(P.at(thisnX-nXg-1)+P.at(thisnX-nXg))/2.0;
      Xgrid.setXmaxFluxBC(FluxMx, P0, 0.0);   
      Xgrid.setXmaxFluxBC(FluxNe, 0.0, 0.0);
   }   
   Xgrid.communicate(FluxN);   
   Xgrid.communicate(FluxMx);   
   Xgrid.communicate(FluxEi);   
   Xgrid.communicate(FluxEe);   
   Xgrid.communicate(FluxNe);   

   
   /*
   cout << "JRA: max(FluxN) = " << max(FluxN) << endl;
   cout << "JRA: min(FluxN) = " << min(FluxN) << endl;
   */

}

void updatePhysicalVars( const domainGrid&  Xgrid )
{
   const int nXcc = Xgrid.Xcc.size();
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   double maxZ, minZ;
   
   Vx  = Mx/N;
   Zbar = Ne/N;
   Nn = N-Ne;
   if(min(N)<0.0)    cout << " N IS LESS THAN ZERO " << endl;
   if(min(Zbar)<Zmin) {
      cout << " Zbar IS LESS THAN Zmin " << endl;
      cout << " min(Zbar) = " << min(Zbar) << endl;
   }
   if(max(Zbar)>1.0)  {
      cout << " Zbar IS LARGER THAN ONE " << endl;
      cout << " max(Zbar) = " << max(Zbar) << endl;
   }

   Pi  = (Ei - 0.5*Vx*Mx)*(gamma0-1.0);
   Pe  = Ee*(gamma0-1.0);
   P = Pe + Pi; 
   Ti  = Pi/N;
   Te  = Pe/Ne;
   if(min(Te)<0.0) cout << " Te IS LESS THAN ZERO " << endl;
   if(min(Ti)<0.0) cout << " Ti IS LESS THAN ZERO " << endl;
   Cs = sqrt(gamma0*P/N);
   Va = sqrt(By*By/N);
   Cspeed  = abs(Vx) + sqrt(Cs*Cs + Va*Va);  // adv flux jacobian
   double Clight = 1.0/sqrt(delta0);
   Cspeed  = min(Cspeed,Clight);

}

void updateCollisionTerms( const domainGrid&  Xgrid )
{
   const int nXcc = Xgrid.Xcc.size();
   const int nXg  = Xgrid.nXg;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   vector<double> nue, taue_spi, Te_eV;
   nue.assign(nXcc,0.0);
   taue_spi.assign(nXcc,0.0);
   Te_eV.assign(nXcc,0.0);
   

   //  set collision times and resistivity
   //
   Te_eV = Te*Tscale; 
   setNeutralInteractionRates(Xgrid,Te_eV,Ne,Nn);

   nue_spi = pow(Tscale,1.5)/taue0*Ne*( 1.0/pow(Te_eV,1.5) );
   nue_vac = pow(Tscale,1.5)/taue0*Ne*( 0.01*pow(NvacC/N,NvacP) );
   nue = nue_spi + nue_neu + nue_vac;
   double nueR_min = pow(Tscale,1.5)/taue0*(1.0/pow(0.1,1.5));
   const int thisnX = nue.size();
   if(procID==0 && XlowBC.compare("insulator") == 0) { // set resistivity at insulator boundary
      if(nue.at(1) < nueR_min) Xgrid.setXminBoundary(nue, nueR_min, 0.0);
   }
   if(procID==numProcs-1 && XhiBC.compare("insulator") == 0) { // set resistivity at insulator boundary
      if(nue.at(thisnX-nXg+1) < nueR_min) Xgrid.setXmaxBoundary(nue, nueR_min, 0.0);
   }
   /*  // trying different vac res method
   nue = nue_spi;
   for (auto i=0; i<thisnX; i++) {
      if(N.at(i)<=0.01) {
         nue.at(i) = nueR_min;
      }
   }
   */
   //eta = eta0*taue0*nue/Ne;   
   eta = Le0sq*nue/Ne;   
   taue = 1.0/nue; // dont use vac taue for thermalization...too small time step needed
   taue_spi = taue0*pow(Te,1.5)/Ne;
   Qie   = 2.0/(gamma0-1.0)*mM*(1.0/taue_spi + nue_neu)*Ne*(Te-Ti); // use taue_spi here for time-step

   Xgrid.InterpToCellEdges(etace,eta,eta,"C2");
   Xgrid.InterpToCellEdges(Nece,Ne,Ne,"C2");
   Xgrid.communicate(etace);
   Xgrid.communicate(eta);

}

void computeJ0( vector<double>&  a_J0_ce,
          const vector<double>&  a_By_cc,
          const domainGrid&      Xgrid )
{
   const int nXg =  Xgrid.nXg;
   const int nXcc = Xgrid.Xcc.size();
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
      
   Xgrid.DDX(a_J0_ce,hy_cc*a_By_cc);
   a_J0_ce = a_J0_ce/hy_ce;
   if(procID==0 && geometry0=="CYL" && XlowBC.compare("axis") == 0) {
      Xgrid.setXminFluxBC(a_J0_ce,2.0*a_By_cc.at(nXg)/hy_cc.at(nXg),0.0);
   } 
   Xgrid.communicate(a_J0_ce);

}

void computeIdealEatEdges( const domainGrid&      Xgrid, 
                           const vector<double>&  a_Cspeed )
{
   const int nXcc = Xgrid.Xcc.size();
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   vector<double> FluxB0cc;
   FluxB0cc.assign(nXcc,0.0);
   FluxB0cc = Vx*By;
   if(advScheme0=="TVD") {
      Xgrid.computeFluxTVD(Eidealz,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxB0cc,a_Cspeed,By,"minmod",Nsub);
   }
   else {
      Xgrid.computeFluxTVDsimple(Eidealz,FluxL,FluxR,FluxB0cc,a_Cspeed,By,advScheme0);
   }
   if(procID==0) Xgrid.setXminFluxBC(Eidealz, 0.0, 0.0);
   if(procID==numProcs-1) Xgrid.setXmaxFluxBC(Eidealz, 0.0, 0.0);
   Xgrid.communicate(Eidealz);
   Eidealz = -1.0*Eidealz;

}

void computeJouleHeatingTerms( const domainGrid&  Xgrid )
{
   const int nXcc = Xgrid.Xcc.size();
   const int nXg  = Xgrid.nXg;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   vector<double> dPedx;
   dPedx.assign(nXcc,0.0);

   // compute current density at cell center
   //
   computeJ0(Jz0,By,Xgrid);

   bool useJ0forJcc = false;
   if(useJ0forJcc) {   
      Xgrid.InterpToCellCenter(Jzcc,Jz0);
      Xgrid.communicate(Jzcc);
   } else {
      Xgrid.InterpToCellCenter(Jzcc,Jz);
      Xgrid.communicate(Jzcc);
   }

   // compute electric field at cell center
   //
   Xgrid.InterpToCellCenter(Ezcc,Ez);
   Xgrid.communicate(Ezcc);
  
   //  calculate work source terms for energy equations
   //
   Xgrid.DDX(dPedx,Pe);
   Xgrid.communicate(dPedx);   
   NUdotE = -Vx*(Jzcc*By + dPedx);
   JdotE = Jzcc*Ezcc;

}

void setNeutralInteractionRates( const domainGrid&     Xgrid, 
                                 const vector<double>&  a_Te_eV, 
                                 const vector<double>&  a_Ne,
                                 const vector<double>&  a_Nn )
{
   const int nXcc = Xgrid.Xcc.size();
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   vector<double> kizn, K3br, kmom;
   kizn.assign(nXcc,0.0);
   K3br.assign(nXcc,0.0);
   kmom.assign(nXcc,1.0e-7); // [cm^3/s]

   //if(!procID) cout << "Te_eV = " << a_Te_eV.at(10) << endl;
   //cout << "Ne    = " << a_Ne*Nscale/1.0e6 << endl;
   //cout << "Nn    = " << a_Nn*Nscale/1.0e6 << endl;

   // compute ionization and 3-body recombination rate constant
   //
   kizn = 1.0e-5*sqrt(a_Te_eV/Uizn)/Uizn/sqrt(Uizn)/(6.0+a_Te_eV/Uizn)*exp(-Uizn/a_Te_eV); // [cm^3/s]
   K3br = kizn*1.66e-22*g1/g0/pow(a_Te_eV,1.5)*exp(Uizn/a_Te_eV); // [cm^6/s]
 
   // set charge density source/sink rate
   //
   nue_izn = a_Nn*Nscale/1.0e6*kizn; // - pow(a_Ne*Nscale/1.0e6,2)*K3br; // [1/s]
   nue_izn = nue_izn*tscale; // normalized

   // compute momentum-exchange rate constant and frequency
   //
   //kmom = 1.0e-7; // [cm^3/s]
   nue_neu = a_Nn*Nscale/1.0e6*kmom; // [1/s] 
   nue_neu = nue_neu*tscale;           // normalized
   nue_neu = nue_neu + nue_izn;         // add growth term

   // compute normalized source terms for charge density and electron energy density
   //
   Se = nue_izn*a_Ne;          // electron density source term
   SEe = -1.0*Uizn/Tscale*Se; // Electron energy density source term

}

void initializeMemberVars( const domainGrid& Xgrid )
{ 
   const int nXce = Xgrid.Xce2.size();
   const int nXcc = Xgrid.Xcc.size();
   const int nXg = Xgrid.nXg;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   N.assign(nXcc,0.0);
   Mx.assign(nXcc,0.0);
   Ei.assign(nXcc,0.0);
   Ee.assign(nXcc,0.0);
   By.assign(nXcc,0.0);
   //
   Nold.assign(nXcc,0.0);
   Mxold.assign(nXcc,0.0);
   Eiold.assign(nXcc,0.0);
   Eeold.assign(nXcc,0.0);
   Byold.assign(nXcc,0.0);
   //
   P.assign(nXcc,0.0);
   Pi.assign(nXcc,0.0);
   Pe.assign(nXcc,0.0);
   Ti.assign(nXcc,0.0);
   Te.assign(nXcc,0.0);
   eta.assign(nXcc,0.0);
   Va.assign(nXcc,0.0);
   Cs.assign(nXcc,0.0);
   Cspeed.assign(nXcc,0.0);
   Vx.assign(nXcc,0.0);
   Jzcc.assign(nXcc,0.0);
   Ezcc.assign(nXcc,0.0);
   Qvisc.assign(nXcc,0.0);
   //
   Qie.assign(nXcc,0.0);
   NUdotE.assign(nXcc,0.0);
   JdotE.assign(nXcc,0.0);
   taue.assign(nXcc,0.0);
   nue_spi.assign(nXcc,0.0);
   nue_vac.assign(nXcc,0.0);
   
   // Ez and Jz are defined on cell edges
   //
   Ez.assign(nXce,0.0);
   Ezold.assign(nXce,0.0);
   Jz.assign(nXce,0.0);
   Jzold.assign(nXce,0.0);
   Jz0.assign(nXce,0.0);
   Jz0stagTime.assign(nXce,0.0);
   etace.assign(nXce,0.0);
   Nece.assign(nXce,0.0);
   Eidealz.assign(nXce,0.0);
   //
   FluxRatio.assign(nXce,0.0);
   FluxLim.assign(nXce,0.0);
   FluxN.assign(nXce,0.0);
   FluxMx.assign(nXce,0.0);
   FluxEe.assign(nXce,0.0);
   FluxEi.assign(nXce,0.0);
   FluxNe.assign(nXce,0.0);
   FluxR.assign(nXce,0.0);
   FluxL.assign(nXce,0.0);
   //
   hy_cc.assign(nXcc,1.0);
   hy_ce.assign(nXce,1.0);

   // ionization terms
   //
   //Zbar.assign(nXcc,Zmin);
   Neold.assign(nXcc,0.0);
   Se.assign(nXcc,0.0);
   SEe.assign(nXcc,0.0);
   Ne.assign(nXcc,0.0);
   Nn.assign(nXcc,0.0);
   nue_neu.assign(nXcc,0.0);
   nue_izn.assign(nXcc,0.0);

}

void addMembersToDataFile( HDF5dataFile&  dataFile )
{
   dataFile.add(N, "N", 1);      // density 
   dataFile.add(Mx, "Mx", 1);    // momentum density 
   dataFile.add(By, "By", 1);    // magnetic field
   dataFile.add(Ei, "Ei", 1);    // total ion energy
   dataFile.add(Ee, "Ee", 1);    // total ele energy
   dataFile.add(P, "P", 1);      // total pressure
   dataFile.add(Pi, "Pi", 1);    // ion pressure
   dataFile.add(Pe, "Pe", 1);    // ele pressure
   dataFile.add(Ti, "Ti", 1);    // ion temperature
   dataFile.add(Te, "Te", 1);    // ele temperature
   dataFile.add(eta, "eta", 1);  // resistivity
   dataFile.add(taue, "taue", 1);  // collision time
   dataFile.add(nue_spi, "nue_spi", 1);  // spitzer coll freq
   dataFile.add(nue_vac, "nue_vac", 1);  // vac coll freq
   
   dataFile.add(Ne,  "Ne", 1);       // Ne = Zbar*N   [cm^3/s]
   dataFile.add(Zbar,"Zbar", 1);     // Ne = Zbar*N   [cm^3/s]
   dataFile.add(nue_neu, "nue_neu", 1);   // e-n coll freq [Hz]
   dataFile.add(nue_izn, "nue_izn", 1);  // e-n izn - recom [Hz]
   
   dataFile.add(Vx, "Vx", 1);     // velocity
   dataFile.add(Jz, "Jz", 1);     // current density
   dataFile.add(Jzcc, "Jzcc", 1); // current density at cell-center
   dataFile.add(Jz0, "Jz0", 1);   // curl of By
   dataFile.add(Ez, "Ez", 1);     // z-electric field
   dataFile.add(Eidealz, "Eidealz", 1);    // z-electric field
   dataFile.add(Va,"Va",1);  // alfven speed
   dataFile.add(Cs,"Cs",1);  // sound speed
   dataFile.add(Cspeed,"Cspeed",1);  // max char speed
   dataFile.add(gamma0,"gamma0",0); 
   dataFile.add(delta0,"delta0",0); 
   dataFile.add(Le0sq,"Le0sq",0); 
   //
   dataFile.add(FluxN,  "FluxN", 1);  
   dataFile.add(FluxMx, "FluxMx", 1);  
   dataFile.add(FluxEi, "FluxEi", 1);  
   dataFile.add(FluxEe, "FluxEe", 1);  
   dataFile.add(FluxNe, "FluxNe", 1);  
   //
   dataFile.add(FluxRatio, "FluxRatio", 1);  
   dataFile.add(FluxLim, "FluxLim", 1);  
   dataFile.add(FluxR, "FluxR", 1);
   dataFile.add(FluxL, "FluxL", 1);
   //
   dataFile.add(Iscale,"Iscale",0);
   dataFile.add(Nscale,"Nscale",0);
   dataFile.add(Tscale,"Tscale",0); 
   dataFile.add(Xscale,"Xscale",0); 
   dataFile.add(Bscale,"Bscale",0);
   dataFile.add(Ezscale,"Ezscale",0);
   dataFile.add(Jscale,"Jscale",0);
   dataFile.add(Pscale,"Pscale",0);
   dataFile.add(Vscale,"Vscale",0);
   dataFile.add(tscale,"tscale",0);
   dataFile.add(Mi,"Mi",0);
   //
   dataFile.add(hy_cc,"hy_cc",0);
   dataFile.add(hy_ce,"hy_ce",0);

}

void parseInputFile( const domainGrid& Xgrid, const Json::Value& a_root )
{ 
   const int nXce = Xgrid.Xce2.size();
   const int nXcc = Xgrid.Xcc.size();
   const int nXg = Xgrid.nXg;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = a_root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      Json::Value advScheme = Phys.get("advScheme",defValue);
      Json::Value geometry  = Phys.get("geometry",defValue);
      Json::Value gammaVal  = Phys.get("gammaC",defValue);
      Json::Value ZminVal   = Phys.get("Zmin",defValue);
      Json::Value NsubVal   = Phys.get("Nsub",defValue);
      //Json::Value etaVal    = Phys.get("eta0",defValue);
      //Json::Value etaVal    = Phys.get("eta0",defValue);
      Json::Value etaVisVal   = Phys.get("etaVis0",defValue);
     
  
      // parse stuff for boundary conditions
      //
      Json::Value XlowBoundary = Phys.get("XlowBoundary",defValue);
      Json::Value XhiBoundary  = Phys.get("XhiBoundary",defValue);
      XlowBC = XlowBoundary.asString();
      XhiBC  = XhiBoundary.asString();
      if(XlowBoundary == defValue || XhiBoundary == defValue) {
         cout << "input ERROR: XlowBoundary or XhiBoundary set incorrectly" << endl;
	 exit (EXIT_FAILURE);
      } else { 
         if(procID==0) {
            cout << "X low boundary is " << XlowBC << endl;
            cout << "X hi boundary is " << XhiBC << endl;
         }
      }
      
 
      //   get characteristic scales from input file
      // 
      Json::Value NscaleVal   = Phys.get("DensScale_invmc",defValue);
      Json::Value XscaleVal   = Phys.get("SpatScale_m",    defValue);
      Json::Value IscaleVal   = Phys.get("CurrScale_Amps", defValue);
      Json::Value dtIscaleVal = Phys.get("dtCurrScale",    defValue);
      Json::Value riseTime_unitsVal = Phys.get("riseTime_units", defValue);
      Json::Value AmassVal    = Phys.get("Amass",          defValue);
      Json::Value NthreshVal  = Phys.get("Nthresh",        defValue);
      Json::Value epsilonRelVal = Phys.get("epsilonRel",defValue);
      Json::Value meRelVal    = Phys.get("meRel",defValue);
      Json::Value NvacCVal    = Phys.get("NvacC",        defValue);
      Json::Value NvacPVal    = Phys.get("NvacP",        defValue);
      
      if(NscaleVal == defValue) {
         cout << "input ERROR: did not set DensScale_invmc correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         Nscale = NscaleVal.asDouble();
         if(procID==0) { 
            cout << endl;
	    cout << "input values:" << endl;
	    cout << "density scale [1/m^3] = " << Nscale << endl;
         }
      }
      //
      if(XscaleVal == defValue) {
         cout << "input ERROR: did not set SpatScale_m correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         Xscale = XscaleVal.asDouble();
         if(procID==0) cout << "spatial scale [m] = " << Xscale << endl;
      }
      //
      if(IscaleVal == defValue) {
         cout << "input ERROR: did not set CurrScale_Amps correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         Iscale = IscaleVal.asDouble();
         if(procID==0) cout << "current scale [Amps] = " << Iscale << endl;
      }
      //
      if(dtIscaleVal == defValue) {
         cout << "input ERROR: did not set dtCurrScale_norm correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         dtIscale = dtIscaleVal.asDouble();
         if(procID==0) cout << "current rise time = " << dtIscale << endl;
      }
      //
      string riseTime_units;
      if(riseTime_unitsVal == defValue) {
         cout << "input ERROR: did not set riseTime_units correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         riseTime_units = riseTime_unitsVal.asString();
         if(procID==0) cout << "current rise time units = " << riseTime_units << endl;
      }

      //
      if(AmassVal == defValue) {
         cout << "input ERROR: did not set Amass correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         Amass = AmassVal.asDouble();
       	 if(procID==0) cout << "atomic mass = " << Amass << endl;
      }
      //
      if(NthreshVal == defValue) {
         cout << "input ERROR: did not set Nthresh correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         Nthresh = NthreshVal.asDouble();
         if(procID==0) cout << "Nthresh = " << Nthresh << endl;
      }
      //
      if(epsilonRelVal == defValue) {
         cout << "input ERROR: did not set epsilonRel correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         epsilonRel = epsilonRelVal.asDouble();
         if(procID==0) cout << "epsilon/epsilon0 = " << epsilonRel << endl;
      }
      //
      if(meRelVal == defValue) {
         cout << "input ERROR: did not set meRel correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         meRel = meRelVal.asDouble();
         if(procID==0) cout << "me/me0 = " << meRel << endl;
      }
      if(NvacCVal != defValue) {
         NvacC = NvacCVal.asDouble();
      }
      if(procID==0) cout << "NvacC = " << NvacC << endl;
      
      if(NvacPVal != defValue) {
         NvacP = NvacPVal.asDouble();
      }
      if(procID==0) cout << "NvacP = " << NvacP << endl;
      
      //   set fundamental constants
      //
      double pi  = 3.14159265; // pi
      double kB  = 1.3807e-23; // Boltzmann constant [J/K]
      double qe  = 1.6022e-19; // electron charge [C]       
      double amu = 1.6605e-27; // atomic mass unit [kg]
      double Mp  = 1.6726e-27; // proton mass [kg]
      double me  = 9.1094e-31; // electron mass [kg]
      double mu0 = pi*4.0e-7;  // permeability of free space [H/m]
      double ep0 = 8.8542e-12; // permittivity of free space [F/m]
      double cvac= 2.9979e8;   // speed of light [m/s]

      //   calculate derived parameter scales
      //
      Pscale  = mu0/2.0*pow(Iscale/2.0/pi/Xscale,2);  // mag pressure at x=R [J/m^3]
      Mi      = Amass*amu;                  // ion mass [kg]
      Vscale  = pow(Pscale/Mi/Nscale,0.5);  // velocity scale [m/s]
      Bscale  = pow(mu0*Pscale,0.5);        // magnetic field scale [T]
      Jscale  = Bscale/Xscale/mu0;          // current density scale [A/m^2]
      Tscale  = Pscale/Nscale/qe;           // temperature scale [eV]
      Ezscale = Vscale*Bscale;              // electric field scale [V/m]
      tscale  = Xscale/Vscale;              // time scale [s]
      double etascale  = Xscale*Xscale*mu0/tscale; // resistivity scale [Ohm-m]
      mM = me/Mi;
      double wpescale = 5.64e4*pow(Nscale/1.0e6,0.5); // ele plasma freq [rad/s]
      double wpiscale = wpescale*pow(me/Mi,0.5);    // ion plasma freq [rad/s]
      double wcescale = qe*Bscale/me;   // ele cyclotron freq [rad/s]
      double wciscale = qe*Bscale/Mi;   // ion cyclotron freq [rad/s]
      double tauescale = 3.44e5/10.0*pow(Tscale,1.5)/(Nscale/1.0e6); // collision time [s]
      double tauiscale = 2.09e7/10.0*pow(Tscale,1.5)/(Nscale/1.0e6)*sqrt(Mi/Mp); // collision time [s]
      //
      double Lescale  = cvac/wpescale;    // ele inertial scale [m]
      double Liscale  = cvac/wpiscale;         // ion inertial scale [m]
      
      if(riseTime_units.compare("ns") == 0 ) dtIscale = dtIscale*1.0e-9/tscale;
      
      if(procID==0) {
         cout << endl;
         cout << "derived scales:" << endl;
         cout << "pressure scale [J/m^3]     = " << Pscale << endl; 
         cout << "velocity scale [m/s]       = " << Vscale << endl; 
         cout << "electric field scale [V/m] = " << Ezscale << endl; 
         cout << "magnetic field scale [T]   = " << Bscale << endl; 
         cout << "temperature scale [eV]     = " << Tscale << endl; 
         cout << "time scale [s]             = " << tscale << endl; 
         cout << "resistivity scale [Ohm-m]  = " << etascale << endl; 
         cout << "ele plasma freq [rad/s]    = " << wpescale << endl; 
         cout << "ion plasma freq [rad/s]    = " << wpiscale << endl; 
         cout << "ele cyclotron freq [rad/s] = " << wcescale << endl; 
         cout << "ion cyclotron freq [rad/s] = " << wciscale << endl; 
         cout << "ele collision time [s]     = " << tauescale << endl; 
         cout << "ion collision time [s]     = " << tauiscale << endl; 
         cout << "ele inertial length [m]    = " << Lescale << endl; 
         cout << "ion inertial length [m]    = " << Liscale << endl; 
      }
      
      //   calculate dimensionless parameters
      //
      Json::Value etaVal = Phys.get("eta0",defValue);
      if(etaVal == defValue) {
         eta0   = 1.03e-4*10.0/pow(Tscale,1.5)/etascale; // norm res
      } else {
	 eta0 = etaVal.asDouble();
      } 
      // 
      taue0 = tauescale/tscale;
      taui0 = tauiscale/tscale;
      delta0 = pow(Vscale/cvac,2.0)*epsilonRel;
      Le0sq = pow(Lescale/Xscale,2.0);
      //B00      = mu0*Iscale/2/pi/R/Bscale;  // 
      B00      = sqrt(2.0);  // 
      if(procID==0) {
	 cout << endl;
         cout << "dimensionless parameters:" << endl;
         cout << "normalized resistivity = " << eta0 << endl;
	 if(etaVal != defValue) cout << "WARNING: USING eta0 FROM INPUT FILE !!!" << endl;
         cout << "taue/tscale = " << taue0 << endl;
         cout << "taui/tscale = " << taui0 << endl;
         cout << "wce*taue = " << wcescale*tauescale << endl;
         cout << "wci*taui = " << wciscale*tauiscale << endl;
         cout << "(Le0/r0)^2 = " << Le0sq << " (Ez relaxation const)" << endl;      
         cout << "(V0/c)^2   = " << delta0   << " (Jz relaxation const)" << endl;      
      }

      if(advScheme == defValue || gammaVal == defValue ||
	 geometry == defValue || NsubVal == defValue || 
         ZminVal == defValue) {
         cout << "ERROR: advScheme or gamma " << endl;
         cout << "or Nsub or geometry or Zmin" << endl;
         cout << "not declared in input file" << endl;
         exit (EXIT_FAILURE);
      } 
      
      advScheme0 = advScheme.asString();
      if(procID==0) {
         cout << "advection diff/interp scheme is " << advScheme0 << endl;
      }
      
      geometry0 = geometry.asString();
      if(geometry0=="CAR" || geometry0=="CYL") {
         if(procID==0) {
            cout << "geometry is " << geometry0 << endl;
         }
	 if(geometry0=="CYL") {
            hy_cc = Xgrid.Xcc;
            hy_ce = Xgrid.Xce2;
	 }
      }
      else {
         cout << "valid geometry types are CAR and CYL" << endl;
         exit (EXIT_FAILURE);
      }
      
      gamma0 = gammaVal.asDouble();
      if(procID==0) cout << "adiabatic coefficent = " << gamma0 << endl;
      if(gamma0 < 1.0) {
         printf("ERROR: adiabatic coefficient can't be < 1\n");
         exit (EXIT_FAILURE);
      }

      Zmin = ZminVal.asDouble();
      if(procID==0) cout << "Zmin = " << Zmin << endl;
      //Zbar;      
      Zbar.assign(nXcc,Zmin);

      Nsub = NsubVal.asInt();
      if(procID==0) cout << "Nsub = " << Nsub << endl;
      if(gamma0 < 1.0) {
         printf("ERROR: Nsub must be int >= 1\n");
         exit (EXIT_FAILURE);
      }
      Pthresh = Nthresh*Tthresh;
      Ethresh = Pthresh/(gamma0-1.0);

      etaVis0 = etaVisVal.asDouble();
      if(procID==0) cout << "viscosity = " << etaVis0 << endl;

   }
   else {
      cout << "value for key \"Physics\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   

   ////////////////////////////////////////////////////////////////////////
   //
   //   get initial profiles for variables
   //
   
   // initialize plasma density
   //
   const Json::Value Nvar = Phys.get("N",defValue);
   if(Nvar.isObject()) { 
      Xgrid.setInitialProfile(N,Nvar);
      if(procID==0) Xgrid.setXminBoundary(N, 0.0, 1.0);   
      if(procID==numProcs-1) Xgrid.setXmaxBoundary(N, 0.0, 1.0);   
      Xgrid.communicate(N);
   } else {
      cout << "value for Physics variable \"N\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   // initialize plasma velocity
   //
   const Json::Value Vvar = Phys.get("Vx",defValue);
   if(Vvar.isObject()) { 
      Xgrid.setInitialProfile(Vx,Vvar);
      if(procID==0) Xgrid.setXminBoundary(Vx, 0.0, -1.0);   
      if(procID==numProcs-1) Xgrid.setXmaxBoundary(Vx, 0.0, -1.0);   
      Xgrid.communicate(Vx);
   } else {
      cout << "value for Physics variable \"Vx\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   // initialize plasma temperature
   //
   const Json::Value Tvar = Phys.get("T",defValue);
   if(Tvar.isObject()) { 
      Xgrid.setInitialProfile(Ti,Tvar);
      Ti = Ti/Tscale; // input temp in eV
      if(procID==0) Xgrid.setXminBoundary(Ti, 0.0, 1.0);   
      if(procID==numProcs-1) Xgrid.setXmaxBoundary(Ti, 0.0, 1.0);   
      Xgrid.communicate(Ti);
      Te = Ti;        // initialize Te=Ti
   } else {
      cout << "value for Physics variable \"T\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   
   // initialize magnetic field
   //
   const Json::Value Bvar = Phys.get("By",defValue);
   if(Bvar.isObject()) { 
      Xgrid.setInitialProfile(By,Bvar);
      //B0 = 0;
      //if(procID==0) Xgrid.setXminBoundary(By, 0.0, 0.0);   
      //if(procID==numProcs-1) Xgrid.setXmaxBoundary_J(By,B0,0);
      //Xgrid.communicate(By);
   } else {
      cout << "value for Physics variable \"By\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

}

void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid, const int verbose)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   double Cmax, nue_izn_max;
   Cmax = max(Cspeed);
   nue_izn_max = max(abs(nue_izn));
   //cout << "Cmax = " << Cmax << endl;
   //cout << "abs(Vx) = " << max(abs(Vx))<< endl;
   //cout << "Cs = " << max(Cs)<< endl;

   const double dX = Xgrid.dX;
   double dtCFL_sound = dX/Cmax;
   double dt_izn = 1.0/nue_izn_max;
   double dtCFL_light = dX*sqrt(delta0);
   double dtmax = min(dtCFL_sound,dtCFL_light);
   dtmax = min(dt_izn,dtmax);
   
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);


   //dtSim = 5.0e-5;
   if(procID==0 && verbose) {
      cout << "sigma_0*dt/delta = " << dtSim/delta0/eta0 << endl;
      cout << "dtCFL_sound = " << dtCFL_sound << endl;
      cout << "dt_izn      = " << dt_izn << endl;
      cout << "dtCFL_light = " << dtCFL_light << endl;
      cout << "dtSim = " << dtSim << endl;
   }
}

