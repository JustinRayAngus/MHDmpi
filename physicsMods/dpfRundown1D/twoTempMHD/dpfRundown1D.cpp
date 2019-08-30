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
 * dN/dt  + d(Mx)/dx = 0
 * dMx/dt + d(Mx*Vx + P)/dx = -Jz*By
 * dEi/dt  + d((Ei+Pi)*Vx)/dx = NUdotE + Qie 
 * dEe/dt  + d((Ee+Pe)*Vx)/dx = Jz*Ez - Qie - NUdotE
 * dBy/dt + d(-Ez)/dx = 0
 * delta0*dEz/dt + d(-By)/dx = -Jz
 * Le0or0sq*dJz/dt = [Ez + Vx*By] - eta*Jz
 *
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
vector<double> N, M, Ei, Ee, B, Ez, Jz;   // time-evolving variables
vector<double> eta, Cs, V, P, Pe, Pi, Te, Ti, Jz0, Qvisc; // derived variables
vector<double> etace, Jzcc, Ezcc, VBce, Nce;
vector<double> Nold, Mold, Eiold, Eeold, Bold, Ezold, Jzold;
vector<double> FluxRatio, FluxLim;
vector<double> FluxR, FluxL;  // flux at cell-edges   
vector<double> FluxN, FluxM, FluxEi, FluxEe, FluxB, FluxEz;
vector<double> Qie, NUdotE, JdotE, taue, nue_spi, nue_vac;

double Nscale, Xscale, Amass, Iscale, dtIscale;
double Pscale, Tscale, Bscale, Jscale, Ezscale, Vscale, tscale, Mi;
double epsilonRel, meRel;

vector<double> hy_cc, hy_ce; // y-dimension lame coefficient

// Set lower threshold values for N, T, and S
double Nthresh=1.0e-4, Tthresh=2.5e-2, Ethresh, Pthresh;

// Set conditions for vaccum resistivity model
double NvacC=0.01, NvacP = 4;

void computeFluxes(const domainGrid&, const int);
void setXminBoundary(vector<double>&, const double, const double);
void setXminExtrap(vector<double>&);
void setXmaxBoundary(vector<double>&, const double, const double);
void setXminBoundaryEz(vector<double>&);
void setXmaxBy(vector<double>&, const double);

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
   //
   N.assign(nXcc,0.0);
   M.assign(nXcc,0.0);
   Ei.assign(nXcc,0.0);
   Ee.assign(nXcc,0.0);
   B.assign(nXcc,0.0);
   //
   Nold.assign(nXcc,0.0);
   Mold.assign(nXcc,0.0);
   Eiold.assign(nXcc,0.0);
   Eeold.assign(nXcc,0.0);
   Bold.assign(nXcc,0.0);
   //
   P.assign(nXcc,0.0);
   Pi.assign(nXcc,0.0);
   Pe.assign(nXcc,0.0);
   Ti.assign(nXcc,0.0);
   Te.assign(nXcc,0.0);
   eta.assign(nXcc,0.0);
   Cs.assign(nXcc,0.0);
   V.assign(nXcc,0.0);
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
   etace.assign(nXce,0.0);
   Nce.assign(nXce,0.0);
   VBce.assign(nXce,0.0);
   //
   FluxRatio.assign(nXce,0.0);
   FluxLim.assign(nXce,0.0);
   FluxN.assign(nXce,0.0);
   FluxM.assign(nXce,0.0);
   FluxEe.assign(nXce,0.0);
   FluxEi.assign(nXce,0.0);
   FluxB.assign(nXce,0.0);
   FluxEz.assign(nXcc,0.0);
   FluxR.assign(nXce,0.0);
   FluxL.assign(nXce,0.0);
   //
   hy_cc.assign(nXcc,1.0);
   hy_ce.assign(nXce,1.0);
   
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
      Json::Value etaVisVal   = Phys.get("etaVis0",defValue);
      
      //   get characteristic scales from input file
      // 
      Json::Value NscaleVal   = Phys.get("DensScale_invmc",defValue);
      Json::Value XscaleVal   = Phys.get("SpatScale_m",    defValue);
      Json::Value IscaleVal   = Phys.get("CurrScale_Amps", defValue);
      Json::Value dtIscaleVal = Phys.get("dtCurrScale",    defValue);
      Json::Value AmassVal    = Phys.get("Amass",          defValue);
      Json::Value NthreshVal  = Phys.get("Nthresh",        defValue);
      Json::Value epsilonRelVal = Phys.get("epsilonRel",defValue);
      Json::Value meRelVal    = Phys.get("meRel",defValue);
      Json::Value NvacCVal    = Phys.get("NvacC",        defValue);
      Json::Value NvacPVal    = Phys.get("NvacP",        defValue);
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
         if(procID==0) cout << "current rise time / tpinch = " << dtIscale << endl;
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
         eta0   = 1.03e-4*10.0/pow(Tscale,1.5)/(Xscale*Xscale*mu0/tscale); // norm res
      } else {
	 eta0 = etaVal.asDouble();
      } 
      // 
      taue0 = tauescale/tscale;
      taui0 = tauiscale/tscale;
      delta0 = pow(Vscale/cvac,2.0)*epsilonRel;
      Le0or0sq = pow(Lescale/Xscale,2.0)*meRel;
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
         cout << "(Le0/r0)^2 = " << Le0or0sq << " (Ez relaxation const)" << endl;      
         cout << "(V0/c)^2   = " << delta0   << " (Jz relaxation const)" << endl;      
      }


      if(advScheme == defValue || gammaVal == defValue ||
	 geometry == defValue || NsubVal == defValue) {
         cout << "ERROR: advScheme or gamma " << endl;
         cout << "or Nsub or geometry" << endl;
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
            hy_ce = Xgrid.Xce;
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
      P = P/Tscale; // input pressure uses eV temp... should clean this up
      if(procID==0) setXminBoundary(P, 0.0, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(P, 0.0, 1.0);   
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
      if(procID==numProcs-1) setXmaxBy(B,B0);
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
   Xgrid.InterpToCellEdges(VBce,V*B,B,"C2");
   Xgrid.InterpToCellEdges(etace,eta,eta,"C2");
   Xgrid.InterpToCellEdges(Nce,N,N,"C2");
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
   dataFile.add(nue_spi, "nue_spi", 1);  // spitzer coll freq
   dataFile.add(nue_vac, "nue_vac", 1);  // vac coll freq
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
   const int nXcc = N.size();
   const int nXg = Xgrid.nXg;
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
	 
	 N.at(i) = Nold.at(i) - thisdt*(FluxN.at(i)-FluxN.at(i-1))/hy_cc.at(i)/Xgrid.dX;
         M.at(i) = Mold.at(i) - thisdt*(FluxM.at(i)-FluxM.at(i-1))/hy_cc.at(i)/Xgrid.dX
		 - thisdt*Jzcc.at(i)*B.at(i);
	 if(geometry0=="CYL") M.at(i) = M.at(i) + thisdt*P.at(i)/hy_cc.at(i);
         Ei.at(i) = Eiold.at(i) - thisdt*(FluxEi.at(i)-FluxEi.at(i-1))/hy_cc.at(i)/Xgrid.dX
      		  + thisdt*(NUdotE.at(i) + Qie.at(i));
         Ee.at(i) = Eeold.at(i) - thisdt*(FluxEe.at(i)-FluxEe.at(i-1))/hy_cc.at(i)/Xgrid.dX
      		  + thisdt*(JdotE.at(i) - Qie.at(i) - NUdotE.at(i));
	 B.at(i) = Bold.at(i) - thisdt*(FluxB.at(i)-FluxB.at(i-1))/Xgrid.dX;
	 
	 if(N.at(i)<=Nthresh) N.at(i) = Nthresh;
	 if(N.at(i)!=N.at(i)) cout << "bout to go bad: N.at(i) = " << N.at(i) << endl;
	 Ethresh = Pthresh/(gamma0-1.0) + 0.5*M.at(i)*M.at(i)/N.at(i);
	 if(Ei.at(i)<=Ethresh) Ei.at(i) =Ethresh;
	 Ethresh = Pthresh/(gamma0-1.0);
	 if(Ee.at(i)<=Ethresh) Ee.at(i) =Ethresh;

      }

      timeDomain* tmesh = timeDomain::tmesh;
      //cout << "thist = " << tmesh->tSim << endl;
      double thist = tmesh->tSim;
      B0 = thist/dtIscale*B00;
      if(B0>B00) B0 = B00;
      //if(thist>=0.03) B0 = B00*pow(0.03/thist,2);

      if(procID==0) {
         //setXminBoundary(N, N.at(2), 0.0);   
         //setXminBoundary(M, -M.at(2), 0.0);   
         //setXminBoundary(E, E.at(2), 0.0);
         setXminBoundary(N, 0.0, 1.0);   
         setXminBoundary(M, 0.0, -1.0);   
         setXminBoundary(Ei, 0.0, 1.0);
         setXminBoundary(Ee, 0.0, 1.0);
         //setXminExtrap(N);
         //setXminExtrap(M);
         //setXminExtrap(E);
         setXminBoundary(B, 0.0, -1.0);   
	 if(N.at(0)<Nthresh || N.at(1)<Nthresh) {
            N.at(0) = Nthresh;
            N.at(1) = Nthresh;
            //cout << "Are we here?" << endl;
	 }
      }
     
      if(procID==numProcs-1) {
         setXmaxBoundary(N, 0.0, 1.0);   
         setXmaxBoundary(M, 0.0, -1.0);   
         setXmaxBoundary(Ei, 0.0, 1.0);   
         setXmaxBoundary(Ee, 0.0, 1.0); 
         //cout << "hy_cc.at(nXcc-nXg) = " << hy_cc.at(nXcc-nXg) << endl;	 
         setXmaxBy(B,B0);   
      }

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
         d0 = delta0/thisdt;	   
         e0 = Le0or0sq/thisdt/Nce.at(i);	   
         //
         Ez.at(i) = -(FluxEz.at(i+1)-FluxEz.at(i))/hy_ce.at(i)/Xgrid.dX + d0*Ezold.at(i)
	          - (e0*Jzold.at(i) + VBce.at(i))/(e0 + etace.at(i));
         Ez.at(i) /= d0 + 1.0/(e0 + etace.at(i)); 
         //
         Jz.at(i) = e0*Jzold.at(i) + Ez.at(i) + VBce.at(i);
         Jz.at(i) /= e0 + etace.at(i);
      }

      //cout << "Xce(nXg-1) =" << Xgrid.Xce.at(nXg-1) << endl;
      //cout << "Xce(nXcc-Xg-1) =" << Xgrid.Xce.at(nXcc-nXg-1) << endl;
      if(procID==0) {
         setXminBoundaryEz(Ez);   
      
         // need to redefine Jz on boundary after BC is set
         //
         d0 = delta0/thisdt;	   
         e0 = Le0or0sq/thisdt/Nce.at(nXg-1);	   
         Jz.at(nXg-1) = e0*Jzold.at(nXg-1) + Ez.at(nXg-1) + VBce.at(nXg-1);
         Jz.at(nXg-1) /= e0 + etace.at(nXg-1);
      }
      if(procID==numProcs-1) {
         Ez.at(nXcc-nXg) = Ez.at(nXcc-nXg-1);
      }   
 
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

   const int nCE = FluxN.size();
   const int nCC = N.size();
   const int nXg = Xgrid.nXg;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   vector<double> Cspeed, FluxNcc, FluxMcc, FluxEicc, FluxEecc;
   vector<double> FluxB0cc, Eprime, nue, taue_spi, Te_eV;
   vector<double> FluxVisc, dVdx, dPedx; 
   FluxNcc.assign(nCC,0.0);
   FluxMcc.assign(nCC,0.0);
   FluxEicc.assign(nCC,0.0);
   FluxEecc.assign(nCC,0.0);
   FluxB0cc.assign(nCC,0.0);
   nue.assign(nCC,0.0);
   taue_spi.assign(nCC,0.0);
   Te_eV.assign(nCC,0.0);
   dPedx.assign(nCC,0.0);
   Eprime.assign(nCE,0.0);
   FluxVisc.assign(nCC,0.0);
   dVdx.assign(nCC,0.0);

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
   Cs = sqrt(gamma0*P/N + B*B/(N+1.0e-4));  // 1.0e-4 here for time step

   
   //  define collision times and resistivity
   //
   Te_eV = Te*Tscale; 
   nue_spi = pow(Tscale,1.5)/taue0*N*( 1.0/pow(Te_eV,1.5) );
   nue_vac = pow(Tscale,1.5)/taue0*N*( 0.01*pow(NvacC/N,NvacP) );
   nue = nue_spi + nue_vac;
   double nueR_min = pow(Tscale,1.5)/taue0*(1.0/pow(0.1,1.5));
   const int thisnX = nue.size();
   if(procID==numProcs-1) { // set resistivity at insulator boundary
      if(nue.at(thisnX-nXg+1) < nueR_min) setXmaxBoundary(nue, nueR_min, 0.0);
   }
   /*  // trying different vac res method
   nue = nue_spi;
   for (auto i=0; i<thisnX; i++) {
      if(N.at(i)<=0.01) {
         nue.at(i) = nueR_min;
      }
   }
   */
   eta = eta0*taue0*nue/N;   
   taue = 1.0/nue; // dont use vac taue for thermalization...too small time step needed
   taue_spi = taue0*pow(Te,1.5)/N;

   Xgrid.InterpToCellEdges(etace,eta,eta,"C2");
   Xgrid.InterpToCellEdges(Nce,N,N,"C2");
   Xgrid.communicate(etace);
   Xgrid.communicate(eta);


   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeed   = abs(V) + Cs; // adv flux jacobian
   FluxNcc = hy_cc*M;
   FluxMcc = hy_cc*(M*V + P);
   FluxEicc = hy_cc*(0.5*V*M + Pi*gamma0/(gamma0-1.0) )*V;
   FluxEecc = hy_cc*Pe*gamma0/(gamma0-1.0)*V;
   FluxB0cc = V*B;
   FluxEz = -hy_cc*B;
   
   // compute viscous terms
   //
   vector<double> etaVisc;
   etaVisc.assign(nCC,0.0); 
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
                           FluxNcc,Cspeed,hy_cc*N,"minmod",Nsub);
      Xgrid.computeFluxTVD(FluxM,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxMcc,Cspeed,hy_cc*M,"minmod",Nsub);
      Xgrid.computeFluxTVD(FluxEi,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxEicc,Cspeed,hy_cc*Ei,"minmod",Nsub);
      Xgrid.computeFluxTVD(FluxEe,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxEecc,Cspeed,hy_cc*Ee,"minmod",Nsub);
      Xgrid.computeFluxTVD(VBce,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxB0cc,abs(V),B,"minmod",Nsub);
   }
   else {
      Xgrid.computeFluxTVDsimple(FluxN,FluxL,FluxR,FluxNcc,Cspeed,hy_cc*N,advScheme0);
      Xgrid.computeFluxTVDsimple(FluxM,FluxL,FluxR,FluxMcc,Cspeed,hy_cc*M,advScheme0);
      Xgrid.computeFluxTVDsimple(FluxEi,FluxL,FluxR,FluxEicc,Cspeed,hy_cc*Ei,advScheme0);
      Xgrid.computeFluxTVDsimple(FluxEe,FluxL,FluxR,FluxEecc,Cspeed,hy_cc*Ee,advScheme0);
      Xgrid.computeFluxTVDsimple(VBce,FluxL,FluxR,FluxB0cc,abs(V),B,advScheme0);
   } 
   //Xgrid.InterpToCellEdges(VBce,FluxB0cc,B,"C2");
   //FluxM = FluxM + FluxVisc;

   if(procID==0) {
      setXminBoundary(FluxN, 0.0, 0.0);
      //cout << "hy_ce.at(1)" << hy_ce.at(1) << endl;   
      setXminBoundary(FluxM, hy_ce.at(nXg-1)*(P.at(nXg)+P.at(nXg-1))/2.0, 0.0);   
      setXminBoundary(FluxEi, 0.0, 0.0);
      setXminBoundary(FluxEe, 0.0, 0.0);
      setXminBoundary(VBce, 0.0, 0.0);
   }
   if(procID==numProcs-1) {
      setXmaxBoundary(FluxN, 0.0, 0.0);
      setXmaxBoundary(FluxEi, 0.0, 0.0);
      setXmaxBoundary(FluxEe, 0.0, 0.0);
      const int thisnX = P.size();
      //cout << "hy_ce.at(thisnX-3)" << hy_ce.at(thisnX-3) << endl;   
      double P0 = hy_ce.at(thisnX-nXg-1)*(P.at(thisnX-nXg-1)+P.at(thisnX-nXg))/2.0;
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
      setXminBoundary(Jz0,2.0*B.at(nXg)/hy_cc.at(nXg),0.0);
   } 
   Xgrid.communicate(Jz0);

   JdotE = Jzcc*Ezcc;
   Qie   = 2.0/(gamma0-1.0)*mM/taue_spi*N*(Te-Ti); // use taue_spi here for time-step
   Xgrid.DDX(dPedx,Pe);
   Xgrid.communicate(dPedx);   
   NUdotE = -V*(Jzcc*B + dPedx);

   Xgrid.communicate(FluxN);   
   Xgrid.communicate(FluxM);   
   Xgrid.communicate(FluxEi);   
   Xgrid.communicate(FluxEe);   
   Xgrid.communicate(FluxB);   
   
   //Xgrid.communicate(FluxL);   
   //Xgrid.communicate(FluxR);   
   //Xgrid.communicate(FluxRatio);   
   //Xgrid.communicate(FluxLim);   

} // end computeFluxes


void setXminBoundary(vector<double>& var, const double C0, const double C1)
{
   
   domainGrid* mesh = domainGrid::mesh;
   const int ishift = mesh->nXg;

   for (auto i=0; i<ishift; i++) {
      //var.at(i) = C0;
      var.at(ishift-i-1) = C0 + C1*var.at(ishift+i);
   }
   /*
   if(C0==0) {
      var.at(1) = 2.0*var.at(2) - var.at(3);
      var.at(0) = var.at(1);
   }
   */
}

void setXminExtrap(vector<double>& var)
{
   
   //var.at(1) = 2.0*var.at(2) - var.at(3);
   //var.at(0) = 2.0*var.at(1) - var.at(2);
   var.at(1) = 3.0*(var.at(2) - var.at(3)) + var.at(4);
   var.at(0) = 3.0*(var.at(1) - var.at(2)) + var.at(3);
   //var.at(0) = var.at(1);
}

void setXmaxBoundary(vector<double>& var, const double C0, const double C1)
{
   const int thisnX = var.size();
   domainGrid* mesh = domainGrid::mesh;
   const int ishift = thisnX-mesh->nXg;
   
   for (auto i=ishift; i<thisnX; i++) {
      var.at(i) = C0 + C1*var.at(2*ishift-i-1);
   }
   //var.back() = C;
   //cout << "var.size() = " var.size() << endl; 
      
}

void setXminBoundaryEz(vector<double>& var)
{
   domainGrid* mesh = domainGrid::mesh;
   const int ishift = mesh->nXg;
   
   var.at(ishift-1) = (3.0*var.at(ishift) - var.at(ishift+1))/2.0;
      
}

void setXmaxBy(vector<double>& var, const double B0)
{
   const int thisnX = var.size();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift = thisnX-mesh->nXg;
   
   for (auto i=ishift; i<thisnX; i++) {
      var.at(i) = B0/hy_cc.at(i);
   }

}


void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid, const int verbose)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   vector<double> Cchar;
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
   if(procID==0 && verbose) {
      cout << "sigma_0*dt/delta = " << dtSim/delta0/eta0 << endl;
      cout << "dtSim = " << dtSim << endl;
   }
}

