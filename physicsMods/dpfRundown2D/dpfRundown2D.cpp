/***
 * 
 * physics module for 1D two Temp MHD dpf rundown
 * (1D railgun)
 * Includes Hall field (Ex) in direction of flow
 * assuming force balance in that direction.
 * Because 1D, there is no dEx/dz and so Hall field
 * only enters in how the ions gain energy
 *
 *
 * dN/dt  + d(Mx)/dx = 0
 * dMx/dt + d(Mx*Vx + P)/dx = -Jz*By
 * dEi/dt  + d((Ei+Pi)*Vx)/dx = NUdotE + Qie 
 * dEe/dt  + d((Ee+Pe)*Vx)/dx = Jz*Ez - Qie - NUdotE
 * dBy/dt + d(-Ez)/dx = 0
 * delta0*dEz/dt + d(-By)/dx = -Jz
 * Le0or0sq*dJz/dt = [Ez + Vx*By] - eta*Jz
 *
 * Vx = Mx/N;
 * Jz = [Ez + Ux*By]/eta
 * Pi = (Ei - 0.5*N*Vx^2)*(gamma0-1)
 * Pe = Ee*(gamma0-1)
 * P = Pi + Pe
 * Ti = Pi/N
 * Te = Pe/N
 * eta = eta0/Te^1.5
 * NUdotE = -Ux*(Jz*By + dPe/dx) (Ex from force balance)
 * Qie = 3.0*me/Mi*N/taue*(Te-Ti)
 * taue = taue0*Te^1.5/N
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
matrix2D<double> N, M, Ei, Ee, B, Ez, Jz;   // time-evolving variables
matrix2D<double> eta, Cs, V, P, Pe, Pi, Te, Ti, Jz0, Qvisc; // derived variables
matrix2D<double> etace, Jzcc, Ezcc, VBce;
matrix2D<double> Nold, Mold, Eiold, Eeold, Bold, Ezold, Jzold;
matrix2D<double> FluxRatio, FluxLim;
matrix2D<double> FluxR, FluxL;  // flux at cell-edges   
matrix2D<double> FluxN, FluxM, FluxEi, FluxEe, FluxB, FluxEz;
matrix2D<double> Qie, NUdotE, JdotE, taue;

double Nscale, Tscale, Xscale, Amass, Iscale, dyIscale, dtIscale;
double Pscale, Tiscale, Tescale, Bscale, Jscale, Ezscale, Vscale, tscale, Mi;
double epsilonRel, meRel;

matrix2D<double> hy_cc, hy_ce; // y-dimension lame coefficient

// Set lower threshold values for N, T, and S
double Nthresh=1.0e-4, Tthresh=2.5e-2, Ethresh, Pthresh;

void computeFluxes(const domainGrid&, const int);
void setXminBoundary(matrix2D<double>&, const double, const double);
void setXminExtrap(matrix2D<double>&, const int);
void setXmaxBoundary(matrix2D<double>&, const double, const double);
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
   M.initialize(nXcc,nZcc,0.0);
   Ei.initialize(nXcc,nZcc,0.0);
   Ee.initialize(nXcc,nZcc,0.0);
   B.initialize(nXcc,nZcc,0.0);
   //
   Nold.initialize(nXcc,nZcc,0.0);
   Mold.initialize(nXcc,nZcc,0.0);
   Eiold.initialize(nXcc,nZcc,0.0);
   Eeold.initialize(nXcc,nZcc,0.0);
   Bold.initialize(nXcc,nZcc,0.0);
   //
   P.initialize(nXcc,nZcc,0.0);
   Pi.initialize(nXcc,nZcc,0.0);
   Pe.initialize(nXcc,nZcc,0.0);
   Ti.initialize(nXcc,nZcc,0.0);
   Te.initialize(nXcc,nZcc,0.0);
   eta.initialize(nXcc,nZcc,0.0);
   Cs.initialize(nXcc,nZcc,0.0);
   V.initialize(nXcc,nZcc,0.0);
   Jzcc.initialize(nXcc,nZcc,0.0);
   Ezcc.initialize(nXcc,nZcc,0.0);
   Qvisc.initialize(nXcc,nZcc,0.0);
   //
   Qie.initialize(nXcc,nZcc,0.0);
   NUdotE.initialize(nXcc,nZcc,0.0);
   JdotE.initialize(nXcc,nZcc,0.0);
   taue.initialize(nXcc,nZcc,0.0);

   // Ez and Jz are defined on cell edges
   //
   Ez.initialize(nXce,nZcc,0.0);
   Ezold.initialize(nXce,nZcc,0.0);
   Jz.initialize(nXce,nZcc,0.0);
   Jzold.initialize(nXce,nZcc,0.0);
   Jz0.initialize(nXce,nZcc,0.0);
   etace.initialize(nXce,nZcc,0.0);
   VBce.initialize(nXce,nZcc,0.0);
   //
   FluxRatio.initialize(nXce,nZcc,0.0);
   FluxLim.initialize(nXce,nZcc,0.0);
   FluxN.initialize(nXce,nZcc,0.0);
   FluxM.initialize(nXce,nZcc,0.0);
   FluxEe.initialize(nXce,nZcc,0.0);
   FluxEi.initialize(nXce,nZcc,0.0);
   FluxB.initialize(nXce,nZcc,0.0);
   FluxEz.initialize(nXcc,nZcc,0.0);
   FluxR.initialize(nXce,nZcc,0.0);
   FluxL.initialize(nXce,nZcc,0.0);
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
            //hy_cc = Xgrid.Xcc;
            //hy_ce = Xgrid.Xce;
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
 
   const Json::Value Vvar = Phys.get("V",defValue);
   if(Vvar.isObject()) { 
      Xgrid.setInitialProfile(V,Vvar);
      if(procID==0) setXminBoundary(V, 0.0, -1.0);   
      if(procID==numProcs-1) setXmaxBoundary(V, 0.0, -1.0);   
      setZboundaryPeriodic(V);
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

   const Json::Value Bvar = Phys.get("B",defValue);
   if(Bvar.isObject()) { 
      Xgrid.setInitialProfile(B,Bvar);
      B0 = 0;
      if(procID==0) setXminBoundary(B, 0.0, 0.0);   
      if(procID==numProcs-1) setXmaxBoundary(B, B0, 0.0);
      setZboundaryPeriodic(B);
      Xgrid.communicate(B);
      Bold = B;
   } else {
      cout << "value for Physics variable \"B\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   Xgrid.DDX(Jz,hy_cc*B);
   //Jz = Jz/hy_ce; 
   Xgrid.DDX(Jzcc,hy_cc*B);
   Jzcc = Jzcc/hy_cc; 
   Xgrid.communicate(Jz);
   Xgrid.communicate(Jzcc);
   Jz0 = Jz;
   Jzold = Jz;
   eta = eta0/pow(Te,1.5);
   taue = taue0*pow(Te,1.5)/N;
   Xgrid.InterpToCellEdges(VBce,V*B,B,"C2",0);
   Xgrid.InterpToCellEdges(etace,eta,eta,"C2",0);
   Ez = etace*Jz-VBce; 
   Ezold = Ez;


   Ei = 0.5*M*M/N + Pi/(gamma0-1.0);
   Ee = Pe/(gamma0-1.0);
   Eiold = Ei;
   Eeold = Ee;

   // need inital flux calculation before add to output
   // and before first call to Physics.advance   
   computeFluxes(Xgrid, 1); // 1 for first order calculation   
  

   // add stuff to output files
   //
   dataFile.add(N, "N", 1);      // density 
   dataFile.add(M, "M", 1);      // momentum density 
   dataFile.add(B, "B", 1);      // magnetic field
   dataFile.add(Ei, "Ei", 1);    // total ion energy
   dataFile.add(Ee, "Ee", 1);    // total ele energy
   dataFile.add(P, "P", 1);      // total pressure
   dataFile.add(Pi, "Pi", 1);    // ion pressure
   dataFile.add(Pe, "Pe", 1);    // ele pressure
   dataFile.add(Ti, "Ti", 1);    // ion temperature
   dataFile.add(Te, "Te", 1);    // ele temperature
   dataFile.add(eta, "eta", 1);  // resistivity
   dataFile.add(taue, "taue", 1);  // collision time
   dataFile.add(V, "V", 1);      // velocity
   dataFile.add(Jz, "J", 1);     // current density
   dataFile.add(Jzcc, "Jcc", 1); // current density at cell-center
   dataFile.add(Jz0, "J0", 1);   // curl of B
   dataFile.add(Ez, "Ez", 1);    // z-electric field
   dataFile.add(Cs,"Cs",1);      // sound speed
   dataFile.add(gamma0,"gamma0",0); 
   dataFile.add(delta0,"delta0",0); 
   dataFile.add(Le0or0sq,"Le0or0sq",0); 
   //
   dataFile.add(FluxN, "FluxN", 1);  
   dataFile.add(FluxM, "FluxM", 1);  
   dataFile.add(FluxEi, "FluxEi", 1);  
   dataFile.add(FluxEe, "FluxEe", 1);  
   dataFile.add(FluxB, "FluxB", 1);  
   dataFile.add(FluxEz, "FluxEz", 1);  
   //
   dataFile.add(FluxRatio, "FluxRatio", 1);  
   dataFile.add(FluxLim, "FluxLim", 1);  
   dataFile.add(FluxR, "FluxR", 1);
   dataFile.add(FluxL, "FluxL", 1);
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
   double thisdt, expFact, Jsp;
   for (auto n=1; n<Nsub+1; n++) {
      thisdt = dt*n/Nsub;

      for (auto i=nXg; i<nXcc-nXg; i++) {
	 for (auto j=nZg; j<nZcc-nZg; j++) {

	 N(i,j) = Nold(i,j) - thisdt*(FluxN(i,j)-FluxN(i-1,j))/hy_cc(i,j)/Xgrid.dX;
         M(i,j) = Mold(i,j) - thisdt*(FluxM(i,j)-FluxM(i-1,j))/hy_cc(i,j)/Xgrid.dX
		- thisdt*Jzcc(i,j)*B(i,j);
	 if(geometry0=="CYL") M(i,j) = M(i,j) + thisdt*P(i,j)/hy_cc(i,j);
         Ei(i,j) = Eiold(i,j) - thisdt*(FluxEi(i,j)-FluxEi(i-1,j))/hy_cc(i,j)/Xgrid.dX
      		  + thisdt*(NUdotE(i,j) + Qie(i,j));
         Ee(i,j) = Eeold(i,j) - thisdt*(FluxEe(i,j)-FluxEe(i-1,j))/hy_cc(i,j)/Xgrid.dX
      		  + thisdt*(JdotE(i,j) - Qie(i,j) - NUdotE(i,j));
	 B(i,j) = Bold(i,j) - thisdt*(FluxB(i,j)-FluxB(i-1,j))/Xgrid.dX;
	 
	 if(N(i,j)<=Nthresh) N(i,j) = Nthresh;
	 if(N(i,j)!=N(i,j)) cout << "bout to go bad: N(i,j) = " << N(i,j) << endl;
	 Ethresh = Pthresh/(gamma0-1.0) + 0.5*M(i,j)*M(i,j)/N(i,j);
	 if(Ei(i,j)<=Ethresh) Ei(i,j) = Ethresh;
	 Ethresh = Pthresh/(gamma0-1.0);
	 if(Ee(i,j)<=Ethresh) Ee(i,j) =Ethresh;

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
         setXminBoundary(M, 0.0, -1.0);   
         setXminBoundary(Ei, 0.0, 1.0);
         setXminBoundary(Ee, 0.0, 1.0);
         //setXminExtrap(N,0);
         //setXminExtrap(M,0);
         //setXminExtrap(E,0);
         setXminBoundary(B, 0.0, -1.0);   
      }
     
      if(procID==numProcs-1) {
         setXmaxBoundary(N, 0.0, 1.0);   
         setXmaxBoundary(M, 0.0, -1.0);   
         setXmaxBoundary(Ei, 0.0, 1.0);   
         setXmaxBoundary(Ee, 0.0, 1.0); 
         //cout << "hy_cc(nXcc-nXg,1) = " << hy_cc(nXcc-nXg,1) << endl;	 
         setXmaxBoundary(B, B0/hy_cc(nXcc-nXg,1), 0.0);   
      }

      setZboundaryPeriodic(N);
      setZboundaryPeriodic(M);
      setZboundaryPeriodic(Ei);
      setZboundaryPeriodic(Ee);
      setZboundaryPeriodic(B);
      Xgrid.communicate(N);
      Xgrid.communicate(M);
      Xgrid.communicate(Ei);
      Xgrid.communicate(Ee);
      Xgrid.communicate(B);

      // compute fluxes using partialy updated fields at n+1/2 for 2nd order
      // or n+1 for 1st order
      computeFluxes(Xgrid, Nsub);


      // Now update electric field and current density
      //
      double d0, e0;
  
      for (auto i=nXg-1; i<nXcc-nXg; i++) {
	 for (auto j=nZg-1; j<nZcc-nZg; j++) {
            d0 = delta0/thisdt;	   
            e0 = Le0or0sq/thisdt/N(i,j);	   
            //
            Ez(i,j) = -(FluxEz(i+1,j)-FluxEz(i,j))/hy_ce(i,j)/Xgrid.dX + d0*Ezold(i,j)
	             - (e0*Jzold(i,j) + VBce(i,j))/(e0 + etace(i,j));
            Ez(i,j) /= d0 + 1.0/(e0 + etace(i,j)); 
            //
            Jz(i,j) = e0*Jzold(i,j) + Ez(i,j) + VBce(i,j);
            Jz(i,j) /= e0 + etace(i,j);
	 }
      }

      //cout << "Xce(nXg-1) =" << Xgrid.Xce.at(nXg-1) << endl;
      //cout << "Xce(nMax-Xg-1) =" << Xgrid.Xce.at(nMax-nXg-1) << endl;
      if(procID==0) {
         setXminBoundaryEz(Ez);   
      }
      //if(procID==numProcs-1) {
      //   setXmaxBoundary(Ez, 0.0, 1.0);   
      //   // setXmaxExtrap(Ez);   
      //}
      
      setZboundaryPeriodic(Ez);
      setZboundaryPeriodic(Jz);
      Xgrid.communicate(Ez);
      Xgrid.communicate(Jz);
   
   } // finish subcycle steps

   // update old fields
   Nold = N;
   Mold = M;
   Eiold = Ei;
   Eeold = Ee;
   Bold = B;
   Ezold = Ez;
   Jzold = Jz;
   
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
  
   matrix2D<double> Cspeed, FluxNcc, FluxMcc, FluxEicc, FluxEecc;
   matrix2D<double> Cspeed2, FluxB0cc, Eprime;
   matrix2D<double> FluxVisc, dVdx, dPedx; 
  
   Cspeed.initialize(nXcc,nZcc,0.0);
   Cspeed2.initialize(nXcc,nZcc,0.0);
   FluxNcc.initialize(nXcc,nZcc,0.0);
   FluxMcc.initialize(nXcc,nZcc,0.0);
   FluxEicc.initialize(nXcc,nZcc,0.0);
   FluxEecc.initialize(nXcc,nZcc,0.0);
   FluxB0cc.initialize(nXcc,nZcc,0.0);
   dPedx.initialize(nXcc,nZcc,0.0);
   Eprime.initialize(nXce,nZcc,0.0);
   FluxVisc.initialize(nXcc,nZcc,0.0);
   dVdx.initialize(nXcc,nZcc,0.0);

   //  define derived variables
   //
   V  = M/N;
   if(min(N)<0.0) cout << " N IS LESS THAN ZERO " << endl;
   Pi  = (Ei - 0.5*V*M)*(gamma0-1.0);
   Pe  = Ee*(gamma0-1.0);
   P = Pe + Pi; 
   Ti  = Pi/N;
   Te  = Pe/N;
   if(min(Te)<0.0) cout << " Te IS LESS THAN ZERO " << endl;
   if(min(Ti)<0.0) cout << " Ti IS LESS THAN ZERO " << endl;
   Cs = sqrt(gamma0*P/N + B*B/N);
   //eta = eta0/Te/sqrt(Te);
   eta = eta0/Te/sqrt(Te)*(1.0+1000.0*pow(0.01/N,4.0));
   //eta = eta0/Te/sqrt(Te)*(1.0+1000.0*pow(0.001/N,4));
   Xgrid.InterpToCellEdges(etace,eta,eta,"C2",0);
   Xgrid.communicate(etace);
   Xgrid.communicate(eta);


   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeed  = abs(V) + Cs; // adv flux jacobian
   FluxNcc = hy_cc*M;
   FluxMcc = hy_cc*(M*V + P);
   FluxEicc = hy_cc*(0.5*V*M + Pi*gamma0/(gamma0-1.0) )*V;
   FluxEecc = hy_cc*Pe*gamma0/(gamma0-1.0)*V;
   FluxB0cc = V*B;
   FluxEz = -hy_cc*B;
   
   // compute viscous terms
   //
   matrix2D<double> etaVisc;
   etaVisc.initialize(nXcc,nZcc,0.0); 
   etaVisc = N*etaVis0;
   //Xgrid.DDX(FluxVisc,-4.0*etaVisc/3.0*V);
   Xgrid.DDX(dVdx,V);
   //Xgrid.communicate(FluxVisc);
   Xgrid.communicate(dVdx);
   Qvisc = 4.0/3.0*etaVisc*dVdx*dVdx;
   FluxVisc = -4.0/3.0*etaVisc*dVdx;
   FluxMcc = FluxMcc + FluxVisc;

   // compute advective flux using
   // specified scheme from input file
   //
   //Nsub = 2;
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxNcc,Cspeed,hy_cc*N,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxM,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxMcc,Cspeed,hy_cc*M,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxEi,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxEicc,Cspeed,hy_cc*Ei,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxEe,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxEecc,Cspeed,hy_cc*Ee,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(VBce,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxB0cc,Cspeed,B,"minmod",0,Nsub);
   }
   else {
      Xgrid.InterpToCellEdges(FluxN,FluxNcc,V,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxM,FluxMcc,V,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxEi,FluxEicc,V,advScheme0,0);
      Xgrid.InterpToCellEdges(FluxEe,FluxEecc,V,advScheme0,0);
      //Xgrid.InterpToCellEdges(FluxB,FluxBcc,V,advScheme0,0);
      //Xgrid.InterpToCellEdges(FluxEz,FluxEzcc,Ez,advScheme0,0);
   } 
   //Xgrid.InterpToCellEdges(VBce,FluxB0cc,B,"C2",0);
   //FluxM = FluxM + FluxVisc;

   if(procID==0) {
      setXminBoundary(FluxN, 0.0, 0.0);
      //cout << "hy_ce.at(1)" << hy_ce.at(1) << endl;   
      setXminBoundary(FluxM, hy_ce(1,1)*(P(2,1)+P(1,1))/2.0, 0.0);   
      setXminBoundary(FluxEi, 0.0, 0.0);
      setXminBoundary(FluxEe, 0.0, 0.0);
      setXminBoundary(VBce, 0.0, 0.0);
   }
   if(procID==numProcs-1) {
      setXmaxBoundary(FluxN, 0.0, 0.0);
      setXmaxBoundary(FluxEi, 0.0, 0.0);
      setXmaxBoundary(FluxEe, 0.0, 0.0);
      const int thisnX = P.size0();
      //cout << "hy_ce.at(thisnX-3)" << hy_ce.at(thisnX-3) << endl;   
      double P0 = hy_ce(thisnX-3,1)*(P(thisnX-3,1)+P(thisnX-2,1))/2.0;
      setXmaxBoundary(FluxM, P0, 0.0);   
      setXmaxBoundary(VBce, 0.0, 0.0);
   }   
   //FluxB = FluxB-Eprime;
   FluxB = -Ez;

   Xgrid.communicate(VBce);
   Eprime = Ez+VBce;
   //Jz = (Ez + VBce)/etace;
   //Xgrid.communicate(Jz);
   Xgrid.InterpToCellCenter(Ezcc,Ez);
   Xgrid.communicate(Ezcc);
   Xgrid.InterpToCellCenter(Jzcc,Jz);
   Xgrid.communicate(Jzcc);
   Xgrid.DDX(Jz0,hy_cc*B);
   Jz0 = Jz0/hy_ce;
   if(procID==0 && geometry0=="CYL") {
      setXminBoundary(Jz0,2.0*B(2,1)/hy_cc(2,1),0.0);
   } 
   Xgrid.communicate(Jz0);

   JdotE = Jzcc*Ezcc;
   taue = taue0*pow(Te,1.5)/N;
   Qie   = 3.0*mM/taue*(Te-Ti);
   Xgrid.DDX(dPedx,Pe);
   Xgrid.communicate(dPedx);   
   NUdotE = -V*(Jzcc*B + dPedx);

   Xgrid.communicate(FluxN);   
   Xgrid.communicate(FluxM);   
   Xgrid.communicate(FluxEi);   
   Xgrid.communicate(FluxEe);   
   Xgrid.communicate(FluxB);   
   
   Xgrid.communicate(FluxL);   
   Xgrid.communicate(FluxR);   
   Xgrid.communicate(FluxRatio);   
   Xgrid.communicate(FluxLim);   

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

void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   matrix2D<double> Cchar(Cs);
   double Cmax;
   Cchar = abs(V)+Cs;
   Cmax = max(Cchar);
   //cout << "Cmax = " << Cmax << endl;
   //cout << "abs(V) = " << max(abs(V))<< endl;
   //cout << "Cs = " << max(Cs)<< endl;

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

