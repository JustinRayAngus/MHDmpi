/***
 * 
 * physics module for 2D dpf rundown using resistive MHD
 * 
 * dN/dt  + d(Mx)/dx + d(Mz)/dz = 0
 * dMx/dt + d(Vx*Mx + P)/dx + d(Vz*Mx)/dz = -Jz*By
 * dMz/dt + d(Vx*Mz)/dx + d(Vz*Mz + P)/dz = Jx*By
 * dE/dt  + d((E+P)*Vx)/dx +d((E+P)*Vz)/dz= JdotE 
 * dBy/dt + d(-Ez)/dx + d(Ex)/dz = 0
 * delta0*dEz/dt = Jz0 - Jz
 * delta0*dEx/dt = Jx0 - Jx
 * Le0or0sq*dJz/dt = N*[Ez + Vx*By] - eta*N*Jz
 * Le0or0sq*dJx/dt = N*[Ex - Vz*By] - eta*N*Jx
 *
 * Vx = Mx/N;
 * Vz = Mz/N;
 * P = (E - 0.5*N*(Vx^2+Vz^2))*(gamma0-1)
 * Pe = P/2.0
 * Pi = P-Pe
 * Ti = Pi/N
 * Te = Pe/N
 * eta = eta0/Te^1.5
 *
 * Jx0 = -d(By)/dz
 * Jz0 = d(By)/dx
 *
 * delta0 = (V0/cvac)^2
 * Le0or0sq = (Le0/r0)^2, Le0 = cvac/wpe0
 * gamma0 = 2/degFreedom + 1
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
#include "matrix2D.h"
#include "domainGrid.h"
#include "timeDomain.h"
#include "HDF5dataFile.h"

#include "Physics.h"


using namespace std;

string advScheme0;  // advection differencing scheme
string geometry0;   // CAR or CYL
double gamma0;      // adiabatic coefficient
double eta0;        // resistivity coefficient
double taue0;       // ele collision time coefficient
double taui0;       // ion collision time coefficient
double mM;          // me/Mi
double etaVis0;     // numerical viscosity coefficient
double delta0;      // relaxation const (V0/cvac)^2
double Le0or0sq;    // normalized electron skin depth squared
double B0, B00;     // boundary value of magnetic field
int Nsub;           // time-solver subcycle steps
matrix2D<double> N, Mx, Mz, E, By, Ez, Jz, Ex, Jx;   // time-evolving variables
matrix2D<double> eta, Cs, Vx, Vz, P, Pe, Pi, Te, Ti, Jz0, Jx0, Qvisc; // derived variables
matrix2D<double> etace_x, etace_z, Jzcc, Ezcc, Jxcc, Excc, VxBy_x, VzBy_z;
matrix2D<double> Nce_x, Nce_z;
matrix2D<double> Nold, Mxold, Mzold, Eold, Byold, Ezold, Exold, Jzold, Jxold;
matrix2D<double> FluxLimR_x, FluxLimL_x, FluxLimR_z, FluxLimL_z;
matrix2D<double> FluxR_x, FluxL_x, FluxR_z, FluxL_z;    
matrix2D<double> FluxN_x, FluxMx_x, FluxMz_x, FluxE_x; 
matrix2D<double> FluxN_z, FluxMx_z, FluxMz_z, FluxE_z; 
matrix2D<double> FluxBy_x, FluxEz_x;
matrix2D<double> FluxBy_z, FluxEx_z;
matrix2D<double> JdotE, taue, deltaP;

double Nscale, Tscale, Xscale, Amass, Iscale, dyIscale, dtIscale;
double Pscale, Tiscale, Tescale, Bscale, Jscale, Ezscale, Vscale, tscale, Mi;
double epsilonRel, meRel;

matrix2D<double> hy_cc, hy_ce; // y-dimension lame coefficient

// Set lower threshold values for N, T, and S
double Nthresh=1.0e-4, Tthresh=2.5e-2, Tmax = 1.0e6, Ethresh, Pthresh, Emax;

void computeFluxes(const domainGrid&, const int);
void setXminBoundary(matrix2D<double>&, const double, const double);
void setXminBoundary(matrix2D<double>&, const vector<double>);
void setXminExtrap(matrix2D<double>&, const int);
void setXmaxExtrap(matrix2D<double>&, const int);
void setXmaxBoundary(matrix2D<double>&, const double, const double);
void setXmaxBoundary(matrix2D<double>&, const vector<double>);
void setXmaxBy(matrix2D<double>&, const double);
void setXminBoundaryEz(matrix2D<double>&);
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
   const int nXcc = Xgrid.nXcc;
   const int nXce = Xgrid.nXce;
   const int nZcc = Xgrid.nZcc;
   const int nZce = Xgrid.nZce;
   
   
   N.initialize(nXcc,nZcc,0.0);
   Mx.initialize(nXcc,nZcc,0.0);
   Mz.initialize(nXcc,nZcc,0.0);
   E.initialize(nXcc,nZcc,0.0);
   By.initialize(nXcc,nZcc,0.0);
   //
   Nold.initialize(nXcc,nZcc,0.0);
   Mxold.initialize(nXcc,nZcc,0.0);
   Mzold.initialize(nXcc,nZcc,0.0);
   Eold.initialize(nXcc,nZcc,0.0);
   Byold.initialize(nXcc,nZcc,0.0);
   //
   P.initialize(nXcc,nZcc,0.0);
   Pi.initialize(nXcc,nZcc,0.0);
   Pe.initialize(nXcc,nZcc,0.0);
   Ti.initialize(nXcc,nZcc,0.0);
   Te.initialize(nXcc,nZcc,0.0);
   eta.initialize(nXcc,nZcc,0.0);
   Cs.initialize(nXcc,nZcc,0.0);
   Vx.initialize(nXcc,nZcc,0.0);
   Vz.initialize(nXcc,nZcc,0.0);
   Jzcc.initialize(nXcc,nZcc,0.0);
   Jxcc.initialize(nXcc,nZcc,0.0);
   Ezcc.initialize(nXcc,nZcc,0.0);
   Excc.initialize(nXcc,nZcc,0.0);
   Qvisc.initialize(nXcc,nZcc,0.0);
   //
   JdotE.initialize(nXcc,nZcc,0.0);
   deltaP.initialize(nXcc,nZcc,0.0);
   taue.initialize(nXcc,nZcc,0.0);

   // Ez and Jz are defined on cell edges
   //
   Ez.initialize(nXce,nZcc,0.0);
   Ezold.initialize(nXce,nZcc,0.0);
   Ex.initialize(nXcc,nZce,0.0);
   Exold.initialize(nXcc,nZce,0.0);
   Jz.initialize(nXce,nZcc,0.0);
   Jz0.initialize(nXce,nZcc,0.0);
   Jzold.initialize(nXce,nZcc,0.0);
   Jx.initialize(nXcc,nZce,0.0);
   Jx0.initialize(nXcc,nZce,0.0);
   Jxold.initialize(nXcc,nZce,0.0);
   etace_x.initialize(nXce,nZcc,0.0);
   etace_z.initialize(nXcc,nZce,0.0);
   Nce_x.initialize(nXce,nZcc,0.0);
   Nce_z.initialize(nXcc,nZce,0.0);
   VxBy_x.initialize(nXce,nZcc,0.0);
   VzBy_z.initialize(nXcc,nZce,0.0);
   //
   FluxLimL_x.initialize(nXce,nZcc,0.0);
   FluxLimR_x.initialize(nXce,nZcc,0.0);
   FluxLimL_z.initialize(nXcc,nZce,0.0);
   FluxLimR_z.initialize(nXcc,nZce,0.0);
   FluxN_x.initialize(nXce,nZcc,0.0);
   FluxMx_x.initialize(nXce,nZcc,0.0);
   FluxMz_x.initialize(nXce,nZcc,0.0);
   FluxE_x.initialize(nXce,nZcc,0.0);
   FluxN_z.initialize(nXcc,nZce,0.0);
   FluxMx_z.initialize(nXcc,nZce,0.0);
   FluxMz_z.initialize(nXcc,nZce,0.0);
   FluxE_z.initialize(nXcc,nZce,0.0);
   FluxBy_x.initialize(nXce,nZcc,0.0);
   FluxEz_x.initialize(nXcc,nZcc,0.0);
   FluxBy_z.initialize(nXcc,nZce,0.0);
   FluxEx_z.initialize(nXcc,nZcc,0.0);
   FluxR_x.initialize(nXce,nZcc,0.0);
   FluxL_x.initialize(nXce,nZcc,0.0);
   FluxR_z.initialize(nXcc,nZce,0.0);
   FluxL_z.initialize(nXcc,nZce,0.0);
   //
   hy_cc.initialize(nXcc,nZcc,1.0);
   hy_ce.initialize(nXce,nZcc,1.0);
   
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      Json::Value advScheme = Phys.get("advScheme",defValue);
      Json::Value geometry  = Phys.get("geometry",defValue);
      Json::Value gammaVal  = Phys.get("gammaC",defValue);
      Json::Value NsubVal   = Phys.get("Nsub",defValue);
      //Json::Value etaVal    = Phys.get("eta0",defValue);
      Json::Value etaVisVal    = Phys.get("etaVis0",defValue);
      
      //   get characteristic scales from input file
      // 
      Json::Value NscaleVal   = Phys.get("DensScale_invmc",defValue);
      Json::Value TscaleVal   = Phys.get("TempScale_eV",defValue);
      Json::Value XscaleVal   = Phys.get("SpatScale_m",defValue);
      Json::Value IscaleVal   = Phys.get("CurrScale_Amps",defValue);
      Json::Value dyIscaleVal = Phys.get("dyCurrScale_m",defValue);
      Json::Value dtIscaleVal = Phys.get("dtCurrScale_ns",defValue);
      Json::Value AmassVal    = Phys.get("Amass",defValue);
      Json::Value NthreshVal  = Phys.get("Nthresh",defValue);
      Json::Value epsilonRelVal = Phys.get("epsilonRel",defValue);
      Json::Value meRelVal    = Phys.get("meRel",defValue);
      //
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
      if(TscaleVal == defValue) {
         cout << "input ERROR: did not set TempScale_eV correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         Tscale = TscaleVal.asDouble();
         if(procID==0) cout << "temperature scale [eV] = " << Tscale << endl;
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
      if(dyIscaleVal == defValue) {
         cout << "input ERROR: did not set dyCurrScale_m correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         dyIscale = dyIscaleVal.asDouble();
         if(procID==0) cout << "in-plane current thickness [m] = " << dyIscale << endl;
      }
      //
      if(dtIscaleVal == defValue) {
         cout << "input ERROR: did not set dtCurrScale_ns correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         dtIscale = dtIscaleVal.asDouble();
         if(procID==0) cout << "current rise time [ns] = " << dtIscale << endl;
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
      //Pscale  = Nscale*Tscale*qe;     // pressure scale [J/m^3]
      Tiscale = 2.0*Tscale;
      Tescale = Tiscale;
      //Pscale  = 2.0*Nscale*Tscale*qe;     // pressure scale [J/m^3]
      Pscale  = Nscale*Tiscale*qe;     // pressure scale [J/m^3]
      Bscale  = pow(mu0*Pscale,0.5);        // magnetic field scale [T]
      Jscale  = Bscale/Xscale/mu0;          // current density scale [A/m^2]
      Mi      = Amass*amu;                  // ion mass [kg]
      Vscale  = pow(Pscale/Mi/Nscale,0.5);  // velocity scale [m/s]
      Ezscale = Vscale*Bscale;              // electric field scale [V/m]
      tscale  = Xscale/Vscale;              // time scale [s]
      double etascale  = Xscale*Xscale*mu0/tscale; // resistivity scale [Ohm-m]
      mM = me/Mi;
      double wpescale = 5.64e4*pow(Nscale/1.0e6,0.5); // ele plasma freq [rad/s]
      double wpiscale = wpescale*pow(me/Mi,0.5);    // ion plasma freq [rad/s]
      double wcescale = qe*Bscale/me;   // ele cyclotron freq [rad/s]
      double wciscale = qe*Bscale/Mi;   // ion cyclotron freq [rad/s]
      double tauescale = 3.44e5/10.0*pow(Tescale,1.5)/(Nscale/1.0e6); // collision time [s]
      double tauiscale = 2.09e7/10.0*pow(Tiscale,1.5)/(Nscale/1.0e6)*sqrt(Mi/Mp); // collision time [s]
      //
      double Lescale  = cvac/wpescale;    // ele inertial scale [m]
      double Liscale  = cvac/wpiscale;         // ion inertial scale [m]


      if(procID==0) {
         cout << endl;
         cout << "derived scales:" << endl;
         cout << "velocity scale [m/s] = " << Vscale << endl; 
         cout << "electric field scale [V/m] = " << Ezscale << endl; 
         cout << "magnetic field scale [T] = " << Bscale << endl; 
         cout << "Ti and Te scale [eV] = " << Tiscale << endl; 
         //cout << "pressure scale [J/m^3] = " << Pscale << endl; 
         cout << "time scale [s] = " << tscale << endl; 
         cout << "resistivity scale [Ohm-m] = " << etascale << endl; 
         cout << "ele plasma freq [rad/s] = " << wpescale << endl; 
         cout << "ion plasma freq [rad/s] = " << wpiscale << endl; 
         cout << "ele cyclotron freq [rad/s] = " << wcescale << endl; 
         cout << "ion cyclotron freq [rad/s] = " << wciscale << endl; 
         cout << "ele collision time [s] = " << tauescale << endl; 
         cout << "ion collision time [s] = " << tauiscale << endl; 
         cout << "ele inertial length [m] = " << Lescale << endl; 
         cout << "ion inertial length [m] = " << Liscale << endl; 
      }


      //   calculate dimensionless parameters
      //
      Json::Value etaVal = Phys.get("eta0",defValue);
      if(etaVal == defValue) {
         eta0   = 1.03e-4/10.0/pow(Tescale,1.5)/(Xscale*Xscale*mu0/tscale); // norm res
      } else {
	 eta0 = etaVal.asDouble();
      } 
      // 
      taue0 = tauescale/tscale;
      taui0 = tauiscale/tscale;
      delta0 = pow(Vscale/cvac,2.0)*epsilonRel;
      Le0or0sq = pow(Lescale/Xscale,2.0)*meRel;
      B00      = mu0*Iscale/dyIscale/Bscale;  // 
      if(procID==0) {
	 cout << endl;
         cout << "dimensionless parameters:" << endl;
         cout << "normalized resistivity = " << eta0 << endl;
	 if(etaVal != defValue) cout << "WARNING: USING eta0 FROM INPUT FILE !!!" << endl;
         cout << "taue/tscale = " << taue0 << endl;
         cout << "taui/tscale = " << taui0 << endl;
         cout << "wce*taue = " << wcescale*tauescale << endl;
         cout << "wci*taui = " << wciscale*tauiscale << endl;
         cout << "(Le0/r0)^2 = " << Le0or0sq << " (Ez relaxation const)"<< endl;      
         cout << "(V0/c)^2 = " << delta0 << " (Jz relaxation const)" << endl;      
         cout << "B(R)/Bscale = " << B00 << endl << endl;      
      }


      if(advScheme == defValue || gammaVal == defValue ||
	 geometry == defValue || NsubVal == defValue) {
         cout << "ERROR: advScheme or gamma " << endl;
         cout << "or Nsub or geometry" << endl;
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
      
      geometry0 = geometry.asString();
      if(geometry0=="CAR" || geometry0=="CYL") {
         if(procID==0) {
            cout << "geometry is " << geometry0 << endl;
         }
	 if(geometry0=="CYL") {
            for (auto i=0; i<nXcc; i++) {
               for (auto j=0; j<nZcc; j++) {
                  hy_cc(i,j) = Xgrid.Xcc.at(i);
                  if(i<nXce) hy_ce(i,j) = Xgrid.Xce.at(i);
               }
            }
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

   const Json::Value Pvar = Phys.get("P",defValue);
   if(Pvar.isObject()) { 
      Xgrid.setInitialProfile(P,Pvar);
      if(procID==0) setXminBoundary(P, 0.0, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(P, 0.0, 1.0);   
      setZboundaryPeriodic(P);
      Xgrid.communicate(P);
      Pe = P/2.0;
      Pi = P-Pe;
      Te = Pe/N;
      Ti = Pi/N;
      Cs = pow(gamma0*P/N,0.5);
   } else {
      cout << "value for Physics variable \"P\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
  
   const Json::Value deltaPvar = Phys.get("deltaP",defValue);
   if(deltaPvar.isObject()) { 
      Xgrid.setInitialProfile(deltaP,deltaPvar);
      if(procID==0) setXminBoundary(deltaP, 0.0, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(deltaP, 0.0, 1.0);   
      setZboundaryPeriodic(deltaP);
      Xgrid.communicate(deltaP);
   } else {
      cout << "value for Physics variable \"deltaP\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   P = P*(1.0 + deltaP);
   N = N*(1.0 + deltaP);

   const Json::Value Byvar = Phys.get("By",defValue);
   if(Byvar.isObject()) { 
      Xgrid.setInitialProfile(By,Byvar);
      B0 = 0;
      if(procID==0) setXminBoundary(By, 0.0, 0.0);   
      if(procID==numProcs-1) setXmaxBoundary(By, B0, 0.0);
      setZboundaryPeriodic(By);
      Xgrid.communicate(By);
      Byold = By;
   } else {
      cout << "value for Physics variable \"B\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   Xgrid.DDX(Jz,hy_cc*By);
   //Jz = Jz/hy_ce; 
   Xgrid.DDX(Jzcc,hy_cc*By);
   Jzcc = Jzcc/hy_cc; 
   Xgrid.communicate(Jz);
   Xgrid.communicate(Jzcc);
   Jz0 = Jz;
   Jzold = Jz;
   eta = eta0/pow(Te,1.5);
   taue = taue0*pow(Te,1.5)/N;
   Xgrid.InterpToCellEdges(VxBy_x,Vx*By,By,"C2",0);
   Xgrid.InterpToCellEdges(VzBy_z,Vz*By,By,"C2",1);
   Xgrid.InterpToCellEdges(etace_x,eta,eta,"C2",0);
   Xgrid.InterpToCellEdges(etace_z,eta,eta,"C2",1);
   Ez = etace_x*Jz-VxBy_x; 
   Ezold = Ez;
   Ex = etace_z*Jx+VzBy_z; 
   Exold = Ex;
  
   E = 0.5*(Mx*Mx+Mz*Mz)/N + P/(gamma0-1.0);
   Eold = E;

   // need inital flux calculation before add to output
   // and before first call to Physics.advance   
   computeFluxes(Xgrid, 1); // 1 for first order calculation   
  

   // add stuff to output files
   //
   dataFile.add(N, "N", 1);      // density 
   dataFile.add(Mx, "Mx", 1);      // momentum density 
   dataFile.add(Mz, "Mz", 1);      // momentum density 
   dataFile.add(By, "By", 1);      // magnetic field
   dataFile.add(E, "E", 1);      // total plasmplasmanergy
   dataFile.add(P, "P", 1);      // total pressure
   dataFile.add(Pi, "Pi", 1);    // ion pressure
   dataFile.add(Pe, "Pe", 1);    // ele pressure
   dataFile.add(Ti, "Ti", 1);    // ion temperature
   dataFile.add(Te, "Te", 1);    // ele temperature
   dataFile.add(eta, "eta", 1);  // resistivity
   dataFile.add(taue, "taue", 1);  // collision time
   dataFile.add(Vx, "Vx", 1);      // velocity
   dataFile.add(Vz, "Vz", 1);      // velocity
   dataFile.add(Jz, "Jz", 1);     // current density
   dataFile.add(Jx, "Jx", 1);     // current density
   dataFile.add(Jzcc, "Jzcc", 1); // current density at cell-center
   dataFile.add(Jxcc, "Jxcc", 1); // current density at cell-center
   dataFile.add(Jz0, "Jz0", 1);   // curl of B
   dataFile.add(Jx0, "Jx0", 1);   // curl of B
   dataFile.add(Ez, "Ez", 1);    // z-electric field
   dataFile.add(Ex, "Ex", 1);    // x-electric field
   dataFile.add(Cs,"Cs",1);      // sound speed
   dataFile.add(gamma0,"gamma0",0); 
   dataFile.add(delta0,"delta0",0); 
   dataFile.add(Le0or0sq,"Le0or0sq",0); 
   //
   dataFile.add(FluxN_x, "FluxN_x", 1);  
   dataFile.add(FluxMx_x, "FluxMx_x", 1);  
   dataFile.add(FluxE_x, "FluxE_x", 1);  
   dataFile.add(FluxBy_x, "FluxBy_x", 1);  
   dataFile.add(FluxEz_x, "FluxEz_x", 1);  
   //
   dataFile.add(FluxLimL_x, "FluxLimL_x", 1);  
   dataFile.add(FluxLimR_x, "FluxLimR_x", 1);  
   dataFile.add(FluxR_x, "FluxR_x", 1);
   dataFile.add(FluxL_x, "FluxL_x", 1);
   //
   dataFile.add(Iscale,"Iscale",0);
   dataFile.add(Nscale,"Nscale",0);
   dataFile.add(Tscale,"Tscale",0); 
   dataFile.add(Tiscale,"Tiscale",0); 
   dataFile.add(Tescale,"Tescale",0); 
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

} // end Physics.initilize


void Physics::advance(const domainGrid& Xgrid, const double dt)
{
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg = Xgrid.nXg;
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

	 N(i,j) = Nold(i,j) - thisdt*(FluxN_x(i,j)-FluxN_x(i-1,j))/hy_cc(i,j)/Xgrid.dX
	                    - thisdt*(FluxN_z(i,j)-FluxN_z(i,j-1))/Xgrid.dZ;
         Mx(i,j)= Mxold(i,j)- thisdt*(FluxMx_x(i,j)-FluxMx_x(i-1,j))/hy_cc(i,j)/Xgrid.dX
                            - thisdt*(FluxMx_z(i,j)-FluxMx_z(i,j-1))/Xgrid.dZ
		            - thisdt*Jzcc(i,j)*By(i,j);
	 if(geometry0=="CYL") {
            Mx(i,j) +=  thisdt*(P(i,j)/hy_cc(i,j) - 4.0/3.0*etaVis0*Nold(i,j)*Vx(i,j)/hy_cc(i,j)/hy_cc(i,j));
         }
         Mz(i,j) = Mzold(i,j) - thisdt*(FluxMz_x(i,j)-FluxMz_x(i-1,j))/hy_cc(i,j)/Xgrid.dX
                              - thisdt*(FluxMz_z(i,j)-FluxMz_z(i,j-1))/Xgrid.dZ
		              + thisdt*Jxcc(i,j)*By(i,j);
         E(i,j) = Eold(i,j) - thisdt*(FluxE_x(i,j)-FluxE_x(i-1,j))/hy_cc(i,j)/Xgrid.dX
                            - thisdt*(FluxE_z(i,j)-FluxE_z(i,j-1))/Xgrid.dZ
      		            + thisdt*JdotE(i,j);
	 By(i,j) = Byold(i,j) - thisdt*(FluxBy_x(i,j)-FluxBy_x(i-1,j))/Xgrid.dX
	                      - thisdt*(FluxBy_z(i,j)-FluxBy_z(i,j-1))/Xgrid.dZ;
	 
	 if(N(i,j)<=Nthresh) N(i,j) = Nthresh;
	 if(N(i,j)!=N(i,j)) {
            cout << "bout to go bad: N(i,j) = " << N(i,j) << endl;
	    exit (EXIT_FAILURE);
	 }
         Ethresh = 2.0*N(i,j)*Tthresh/(gamma0-1.0) + 0.5*(Mx(i,j)*Mx(i,j)+Mz(i,j)*Mz(i,j))/N(i,j);
	 Emax    = 2.0*N(i,j)*Tmax/(gamma0-1.0) + 0.5*(Mx(i,j)*Mx(i,j)+Mz(i,j)*Mz(i,j))/N(i,j);
	 if(E(i,j)<=Ethresh) E(i,j) = Ethresh;
	 if(E(i,j)>=Emax) E(i,j) = Emax;

	 }
      }
     
      timeDomain* tmesh = timeDomain::tmesh;
      //cout << "thist = " << tmesh->tSim << endl;
      double thist = tmesh->tSim;
      B0 = thist*3000.0;
      //if(B0>25.0) B0 = 25.0;
      if(B0>B00) B0 = B00;
      
      if(procID==0) {
         setXminBoundary(N, 0.0, 1.0);   
         setXminBoundary(Mx, 0.0, -1.0);   
         setXminBoundary(Mz, 0.0, 1.0);   
         setXminBoundary(E, 0.0, 1.0);
         //setXminExtrap(N,0);
         //setXminExtrap(Mx,0);
         //setXminExtrap(E,0);
         setXminBoundary(By, 0.0, -1.0);   
      }
     
      if(procID==numProcs-1) {
         setXmaxBoundary(N, 0.0, 1.0);   
         setXmaxBoundary(Mx, 0.0, -1.0);   
         setXmaxBoundary(Mz, 0.0, 1.0);   
         setXmaxBoundary(E, 0.0, 1.0);   
         //cout << "hy_cc(nXcc-nXg,1) = " << hy_cc(nXcc-nXg,1) << endl;	 
         setXmaxBoundary(By, B0/hy_cc(nXcc-nXg,1), 0.0);   
         //
         //setXmaxExtrap(N,0);
         //setXmaxExtrap(Mx,0);
         //setXmaxExtrap(E,0);
      }

      setZboundaryPeriodic(N);
      setZboundaryPeriodic(Mx);
      setZboundaryPeriodic(Mz);
      setZboundaryPeriodic(E);
      setZboundaryPeriodic(By);
      Xgrid.communicate(N);
      Xgrid.communicate(Mx);
      Xgrid.communicate(Mz);
      Xgrid.communicate(E);
      Xgrid.communicate(By);

      // compute fluxes using partialy updated fields at n+1/2 for 2nd order
      // or n+1 for 1st order
      computeFluxes(Xgrid, Nsub);


      // Now update electric field and current density
      //
      double d0, e0;
  
      for (auto i=nXg-1; i<nXcc-nXg; i++) {
	 for (auto j=nZg-1; j<nZcc-nZg; j++) {
            d0 = delta0/thisdt;	   
            e0 = Le0or0sq/thisdt/Nce_x(i,j);	   
            //
            Ez(i,j) = -(FluxEz_x(i+1,j)-FluxEz_x(i,j))/hy_ce(i,j)/Xgrid.dX + d0*Ezold(i,j)
	              - (e0*Jzold(i,j) + VxBy_x(i,j))/(e0 + etace_x(i,j));
            Ez(i,j) /= d0 + 1.0/(e0 + etace_x(i,j)); 
            //
            Jz(i,j) = e0*Jzold(i,j) + Ez(i,j) + VxBy_x(i,j);
            Jz(i,j) /= e0 + etace_x(i,j);
            //
            e0 = Le0or0sq/thisdt/Nce_z(i,j);	   
            Ex(i,j) = -(FluxEx_z(i,j+1)-FluxEx_z(i,j))/Xgrid.dZ + d0*Exold(i,j)
	              - (e0*Jxold(i,j) - VzBy_z(i,j))/(e0 + etace_z(i,j));
            Ex(i,j) /= d0 + 1.0/(e0 + etace_z(i,j)); 
            //
            Jx(i,j) = e0*Jxold(i,j) + Ex(i,j) - VzBy_z(i,j);
            Jx(i,j) /= e0 + etace_z(i,j);
	 }
      }

      if(procID==0) {
         setXminBoundaryEz(Ez);   
         //setXminBoundary(Jz,0.0,0.0);   

         // redefine Jz at r=0 since Ez used previously was 0/0 = NAN
         //
         for (auto j=nZg-1; j<nZcc-nZg; j++) {
            d0 = delta0/thisdt;	   
            e0 = Le0or0sq/thisdt/Nce_x(nXg-1,j);	   
            //
            Jz(nXg-1,j) = e0*Jzold(nXg-1,j) + Ez(nXg-1,j) + VxBy_x(nXg-1,j);
            Jz(nXg-1,j) /= e0 + etace_x(nXg-1,j);
         }

      }
      
      setZboundaryPeriodic(Ez);
      setZboundaryPeriodic(Jz);
      setZboundaryPeriodic(Ex);
      setZboundaryPeriodic(Jx);
      Xgrid.communicate(Ez);
      Xgrid.communicate(Jz);
      Xgrid.communicate(Ex);
      Xgrid.communicate(Jx);
   
   } // finish subcycle steps

   // update old fields
   Nold = N;
   Mxold = Mx;
   Mzold = Mz;
   Eold = E;
   Byold = By;
   Ezold = Ez;
   Exold = Ex;
   Jzold = Jz;
   Jxold = Jx;
   
   // compute fluxes using fully updated fields at n+1
   computeFluxes(Xgrid, 1);

} // end Physics.advance


void computeFluxes(const domainGrid& Xgrid, const int order)
{

   const int nXcc = Xgrid.nXcc;
   const int nXce = Xgrid.nXce;
   const int nZcc = Xgrid.nZcc;
   const int nZce = Xgrid.nZce;

   //cout << "ARE WE HERE YET? " << endl;
 
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
  
   matrix2D<double> CspeedBx, Cspeedx, FluxNcc_x, FluxMxcc_x, FluxMzcc_x, FluxEcc_x;
   matrix2D<double> CspeedBz, Cspeedz, FluxNcc_z, FluxMxcc_z, FluxMzcc_z, FluxEcc_z;
   matrix2D<double> FluxBy0cc_x, FluxBy0cc_z, Exprime, Ezprime;
   matrix2D<double> FluxVisc_x, FluxVisc_z, dVxdx, dVzdz, dPedx; 
  
   Cspeedx.initialize(nXcc,nZcc,0.0);
   Cspeedz.initialize(nXcc,nZcc,0.0);
   CspeedBx.initialize(nXcc,nZcc,0.0);
   CspeedBz.initialize(nXcc,nZcc,0.0);
   FluxNcc_x.initialize(nXcc,nZcc,0.0);
   FluxMxcc_x.initialize(nXcc,nZcc,0.0);
   FluxMzcc_x.initialize(nXcc,nZcc,0.0);
   FluxEcc_x.initialize(nXcc,nZcc,0.0);
   FluxNcc_z.initialize(nXcc,nZcc,0.0);
   FluxMxcc_z.initialize(nXcc,nZcc,0.0);
   FluxMzcc_z.initialize(nXcc,nZcc,0.0);
   FluxEcc_z.initialize(nXcc,nZcc,0.0);
   FluxBy0cc_x.initialize(nXcc,nZcc,0.0);
   FluxBy0cc_z.initialize(nXcc,nZcc,0.0);
   dPedx.initialize(nXcc,nZcc,0.0);
   Exprime.initialize(nXcc,nZce,0.0);
   Ezprime.initialize(nXce,nZcc,0.0);
   FluxVisc_x.initialize(nXce,nZcc,0.0);
   FluxVisc_z.initialize(nXcc,nZce,0.0);
   dVxdx.initialize(nXce,nZcc,0.0);
   dVzdz.initialize(nXcc,nZce,0.0);

   //  define derived variables
   //
   Vx  = Mx/N;
   Vz  = Mz/N;
   if(min(N)<0.0) cout << " N IS LESS THAN ZERO " << endl;
   P  = (E - 0.5*Vx*Mx - 0.5*Vz*Mz)*(gamma0-1.0);
   Pe  = P/2.0;
   Pi = P - Pe; 
   Ti  = Pi/N;
   Te  = Pe/N;
   if(min(Te)<0.0) cout << " Te IS LESS THAN ZERO " << endl;
   if(min(Ti)<0.0) cout << " Ti IS LESS THAN ZERO " << endl;
   Cs = sqrt(gamma0*P/N + By*By/N);
   //eta = eta0/Te/sqrt(Te);
   eta = eta0/Te/sqrt(Te)*(1.0+1000.0*pow(0.01/N,4.0));
   //eta = eta0/Te/sqrt(Te) + eta0*3.0e3*pow(Nthresh/N,2.0);
   //eta = eta0/Te/sqrt(Te) + eta0*3.0e7*pow(Nthresh/N,3.0);
   Xgrid.InterpToCellEdges(etace_x,eta,eta,"C2",0);
   Xgrid.InterpToCellEdges(etace_z,eta,eta,"C2",1);
   Xgrid.communicate(etace_x);
   Xgrid.communicate(etace_z);
   Xgrid.InterpToCellEdges(Nce_x,N,N,"C2",0);
   Xgrid.InterpToCellEdges(Nce_z,N,N,"C2",1);
   Xgrid.communicate(Nce_x);
   Xgrid.communicate(Nce_z);
   Xgrid.communicate(eta);


   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   CspeedBx  = abs(Vx); // adv flux jacobian
   CspeedBz  = abs(Vz); // adv flux jacobian
   Cspeedx  = abs(Vx) + Cs; // adv flux jacobian
   Cspeedz  = abs(Vz) + Cs; // adv flux jacobian
   FluxNcc_x = hy_cc*Mx;
   FluxNcc_z = Mz;
   FluxMxcc_x = hy_cc*(Mx*Vx + P);
   FluxMxcc_z = Mx*Vz;
   FluxMzcc_x = hy_cc*Mz*Vx;
   FluxMzcc_z = Mz*Vz + P;
   FluxEcc_x = hy_cc*(E + P)*Vx;
   FluxEcc_z = (E + P)*Vz;
   FluxBy0cc_x = Vx*By;
   FluxBy0cc_z = Vz*By;
   FluxEz_x = -hy_cc*By;
   FluxEx_z = By;
   
   // compute viscous terms
   //
   matrix2D<double> etaVisc;
   etaVisc.initialize(nXcc,nZcc,0.0); 
   etaVisc = N*etaVis0;
   //Xgrid.DDX(FluxVisc,-4.0*etaVisc/3.0*V);
   Xgrid.DDX(dVxdx,Vx);
   Xgrid.DDZ(dVzdz,Vz);
   //Xgrid.communicate(FluxVisc);
   Xgrid.communicate(dVxdx);
   Xgrid.communicate(dVzdz);
   //Qvisc = 4.0/3.0*etaVis0*(dVxdx*dVxdx + dVzdz*dVzdz);
   FluxVisc_x = -4.0/3.0*etaVis0*Nce_x*dVxdx*hy_ce;
   FluxVisc_z = -4.0/3.0*etaVis0*Nce_z*dVzdz;
   //FluxMxcc_x = FluxMxcc_x + FluxVisc;
   //FluxMzcc_z = FluxMzcc_z + FluxVisc_z;

   // compute advective flux using
   // specified scheme from input file
   //
   //Nsub = 2;
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN_x,FluxL_x,FluxR_x,FluxLimL_x,FluxLimR_x,
                           FluxNcc_x,Cspeedx,hy_cc*N,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxMx_x,FluxL_x,FluxR_x,FluxLimL_x,FluxLimR_x,
                           FluxMxcc_x,Cspeedx,hy_cc*Mx,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxMz_x,FluxL_x,FluxR_x,FluxLimL_x,FluxLimR_x,
                           FluxMzcc_x,Cspeedx,hy_cc*Mz,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxE_x,FluxL_x,FluxR_x,FluxLimL_x,FluxLimR_x,
                           FluxEcc_x,Cspeedx,hy_cc*E,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(VxBy_x,FluxL_x,FluxR_x,FluxLimL_x,FluxLimR_x,
                           FluxBy0cc_x,CspeedBx,By,"minmod",0,Nsub);
      //
      Xgrid.computeFluxTVD(FluxN_z,FluxL_z,FluxR_z,FluxLimL_z,FluxLimR_z,
                           FluxNcc_z,Cspeedz,N,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxMx_z,FluxL_z,FluxR_z,FluxLimL_z,FluxLimR_z,
                           FluxMxcc_z,Cspeedz,Mx,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxMz_z,FluxL_z,FluxR_z,FluxLimL_z,FluxLimR_z,
                           FluxMzcc_z,Cspeedz,Mz,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxE_z,FluxL_z,FluxR_z,FluxLimL_z,FluxLimR_z,
                           FluxEcc_z,Cspeedz,E,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(VzBy_z,FluxL_z,FluxR_z,FluxLimL_z,FluxLimR_z,
                           FluxBy0cc_z,CspeedBz,By,"minmod",1,Nsub);
   }
   else {
      Xgrid.InterpToCellEdges(FluxN_x,FluxNcc_x,Vx,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxMx_x,FluxMxcc_x,Vx,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxE_x,FluxEcc_x,Vx,advScheme0,0);
      //Xgrid.InterpToCellEdges(FluxBy_x,FluxBycc,Vx,advScheme0,0);
      //Xgrid.InterpToCellEdges(FluxEz_x,FluxEzcc,Ez,advScheme0,0);
   } 
   //Xgrid.InterpToCellEdges(VxBy_x,FluxBy0cc_x,By,"C2",0);
   //Xgrid.InterpToCellEdges(VzBy_z,FluxBy0cc_z,By,"C2",1);
   FluxMx_x += FluxVisc_x;
   FluxMz_z += FluxVisc_z;

   vector<double> P0;
   P0.assign(nZcc,0.0);
   if(procID==0) {
      setXminBoundary(FluxN_x, 0.0, 0.0);
      //cout << "hy_ce(1,1)" << hy_ce(1,1) << endl;   
      for (auto j=0; j<nZcc; j++) { 
         P0.at(j) = hy_ce(1,j)*(P(2,j)+P(1,j))/2.0;
      }
      //setXminBoundary(FluxMx_x, hy_ce(1,1)*(P(2,1)+P(1,1))/2.0, 0.0);   
      setXminBoundary(FluxMx_x, P0);   
      setXminBoundary(FluxMz_x, 0.0,0.0);   
      setXminBoundary(FluxE_x, 0.0, 0.0);
      setXminBoundary(VxBy_x, 0.0, 0.0);
   }
   if(procID==numProcs-1) {
      setXmaxBoundary(FluxN_x, 0.0, 0.0);
      setXmaxBoundary(FluxE_x, 0.0, 0.0);
      const int thisnX = P.size0();
      //cout << "hy_ce(thisnX-3,1)" << hy_ce(thisnX-3,1) << endl;   
      for (auto j=0; j<nZcc; j++) { 
         P0.at(j) = hy_ce(thisnX-3,j)*(P(thisnX-3,j)+P(thisnX-2,j))/2.0;
      }
      setXmaxBoundary(FluxMx_x, P0);   
      setXmaxBoundary(FluxMz_x, 0.0,0.0);   
      setXmaxBoundary(VxBy_x, 0.0, 0.0);
   }   
   //FluxBy = FluxBy-Eprime;
   FluxBy_x = -Ez;
   FluxBy_z = Ex;

   Xgrid.communicate(VxBy_x);
   Ezprime = Ez+VxBy_x;
   Exprime = Ex-VzBy_z;
   //Jz = (Ez + VxBy_x)/etace;
   //Xgrid.communicate(Jz);
   Xgrid.InterpToCellCenter(Ezcc,Ez);
   Xgrid.communicate(Ezcc);
   Xgrid.InterpToCellCenter(Jzcc,Jz);
   Xgrid.communicate(Jzcc);
   Xgrid.DDX(Jz0,hy_cc*By);
   Jz0 = Jz0/hy_ce;
   if(procID==0 && geometry0=="CYL") {
      for (auto j=0; j<nZcc; j++) { 
         P0.at(j) = 2.0*By(2,j)/hy_cc(2,j);
      }
      setXminBoundary(Jz0,P0);
      //setXminBoundary(Jz0,2.0*By(2,1)/hy_cc(2,1),0.0);
   } 
   Xgrid.communicate(Jz0);
   Xgrid.DDZ(Jx0,By);
   Jx0 = -1.0*Jx0;
   Xgrid.InterpToCellCenter(Jxcc,Jx);
   Xgrid.InterpToCellCenter(Excc,Ex);

   JdotE = Jzcc*Ezcc + Jxcc*Excc;
   taue = taue0*pow(Te,1.5)/N;
   //Xgrid.DDX(dPedx,Pe);
   //Xgrid.communicate(dPedx);   
   //NUdotE = -Vx*(Jzcc*By + dPedx);

   
   Xgrid.communicate(FluxN_x);   
   Xgrid.communicate(FluxMx_x);   
   Xgrid.communicate(FluxMz_x);   
   Xgrid.communicate(FluxE_x);   
   Xgrid.communicate(FluxBy_x);   
   
   //Xgrid.communicate(FluxL_x);   
   //Xgrid.communicate(FluxR_x);   
   //Xgrid.communicate(FluxLimR_x);   
   //Xgrid.communicate(FluxLimL_x);   

} // end computeFluxes


void setXminBoundary(matrix2D<double>& var, const double C0, const double C1)
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

void setXminBoundary(matrix2D<double>& var, const vector<double> C0)
{

   const int thisnZ = var.size1();
   domainGrid* mesh = domainGrid::mesh;
   const int ishift = mesh->nXg;

   for (auto i=0; i<ishift; i++) {
      for (auto j=0; j<thisnZ; j++) {
	 var(ishift-i-1,j) = C0.at(j);
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
   int ishift = thisnX - mesh->nXg;

   for (auto j=0; j<thisnZ; j++) {
      if(order==0) {
         var(ishift,j) = var(ishift-1,j);
         var(ishift+1,j) = var(ishift-1,j);
      }
      else {
         var(ishift,j) = 3.0*(var(ishift-1,j) - var(ishift-2,j)) + var(ishift-3,j);
         ishift = ishift + 1;
         var(ishift,j) = 3.0*(var(ishift-1,j) - var(ishift-2,j)) + var(ishift-3,j);
      }
   }

}

void setXmaxBoundary(matrix2D<double>& var, const double C0, const double C1)
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
   //var.back() = C;
   //cout << "var.size() = " var.size() << endl; 
      
}

void setXmaxBoundary(matrix2D<double>& var, const vector<double> C0)
{
   const int thisnX = var.size0();
   const int thisnZ = var.size1();
   domainGrid* mesh = domainGrid::mesh;
   const int ishift = thisnX-mesh->nXg;
   
   for (auto i=ishift; i<thisnX; i++) {
      for (auto j=0; j<thisnZ; j++) {
         var(i,j) = C0.at(j);
      }
   }
   //var.back() = C;
   //cout << "var.size() = " var.size() << endl; 
      
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

   const int thisnZ = var.size1();
   domainGrid* mesh = domainGrid::mesh;
   const int ishift = mesh->nXg;
   
   for (auto j=0; j<thisnZ; j++) {
      var(ishift-1,j) = (3.0*var(ishift,j) - var(ishift+1,j))/2.0;
      //var(ishift-1,j) = 4.0/3.0*var(ishift,j) - 1.0/3.0*var(ishift+1,j);
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
   //cout << "abs(V) = " << max(abs(V))<< endl;
   //cout << "Cs = " << max(Cs)<< endl;

   const double dX = Xgrid.dX;
   const double dZ = Xgrid.dZ;

   double dtCFL_sound = 0.5*dX*dZ/(dX+dZ)/Cmax;
   double dtCFL_light = dX*sqrt(delta0);
   double dtmax = min(dtCFL_sound,dtCFL_light);
   //double dtmax = dtCFL_sound;
   
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);


   //dtSim = 1.0e-4;
   if(procID==0 && verbose) {
      cout << "sigma_0*dt/delta = " << dtSim/delta0/eta0 << endl;
      cout << "dtSim = " << dtSim << endl;
   }

}



