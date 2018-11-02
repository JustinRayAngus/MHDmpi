/***
 * 
 * physics module for 1D resistive MHD dpf rundown
 * (1D railgun)
 *
 * dN/dt  + d(Mx)/dx = 0
 * dMx/dt + d(Mx*Vx + P)/dx = -Jz*By
 * dE/dt  + d((E+P)*Vx)/dx = Jz*Ez 
 * dBy/dt + d(-Ez)/dx = 0
 * delta0*dEz/dt + d(-By)/dx = -Jz
 * Le0or0sq*dJz/dt = [Ez + Vx*By] - eta*Jz
 *
 * Vx = Mx/N;
 * E  = 0.5*N*Vx^2 + P/(gamma0-1);
 * J = [E + U*B]/eta
 * P = (E - 0.5*N*Vx^2 - By^2/2)*(gamma0-1)
 * T = P/2.0/N
 * eta = eta0/T^1.5
 *
 * delta0 = (V0/cvac)^2
 * Le0or0sq = (Le0/r0)^2, Le0 = cvac/wpe0
 * gamma0 = 2/degFreedom + 1
 * eta0 = 1.03/10/T0^1.5/(r0^2*mu0/t0)
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

string advScheme0;  // advection differencing scheme
double gamma0;      // adiabatic coefficient
double eta0;        // resistivity coefficient
double etaVis0;     // numerical viscosity coefficient
double delta0;      // relaxation const (V0/cvac)^2
double Le0or0sq;    // normalized electron skin depth squared
double B0;          // boundary value of magnetic field
int Nsub;           // time-solver subcycle steps
vector<double> N, M, E, B, Ez, Jz;   // time-evolving variables
vector<double> eta, Cs, V, P, T, Jz0, Qvisc; // derived variables
vector<double> etace, Jzcc, Ezcc, VBce;
vector<double> Nold, Mold, Eold, Bold, Ezold, Jzold;
vector<double> FluxRatio, FluxLim;
vector<double> FluxR, FluxL;  // flux at cell-edges   
vector<double> FluxN, FluxM, FluxE, FluxB, FluxEz;

double Nscale, Tscale, Xscale, Amass, Iscale, dyIscale, dtIscale;
double epsilonRel, meRel;

// Set lower threshold values for N, T, and S
double Nthresh=1.0e-4, Tthresh=2.5e-2, Ethresh, Pthresh;

void computeFluxes(const domainGrid&, const int);
void setXminBoundary(vector<double>&, const double, const double);
void setXminExtrap(vector<double>&);
void setXmaxBoundary(vector<double>&, const double, const double);

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
   E.assign(nXcc,0.0);
   B.assign(nXcc,0.0);
   //
   Nold.assign(nXcc,0.0);
   Mold.assign(nXcc,0.0);
   Eold.assign(nXcc,0.0);
   Bold.assign(nXcc,0.0);
   //
   P.assign(nXcc,0.0);
   T.assign(nXcc,0.0);
   eta.assign(nXcc,0.0);
   Cs.assign(nXcc,0.0);
   V.assign(nXcc,0.0);
   Jzcc.assign(nXcc,0.0);
   Ezcc.assign(nXcc,0.0);
   Qvisc.assign(nXcc,0.0);

   // Ez and Jz are defined on cell edges
   //
   Ez.assign(nXce,0.0);
   Ezold.assign(nXce,0.0);
   Jz.assign(nXce,0.0);
   Jzold.assign(nXce,0.0);
   Jz0.assign(nXce,0.0);
   etace.assign(nXce,0.0);
   VBce.assign(nXce,0.0);
   //
   FluxRatio.assign(nXce,0.0);
   FluxLim.assign(nXce,0.0);
   FluxN.assign(nXce,0.0);
   FluxM.assign(nXce,0.0);
   FluxE.assign(nXce,0.0);
   FluxB.assign(nXce,0.0);
   FluxEz.assign(nXcc,0.0);
   FluxR.assign(nXce,0.0);
   FluxL.assign(nXce,0.0);
   //
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      Json::Value advScheme = Phys.get("advScheme",defValue);
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
      double me  = 9.1094e-31; // electron mass [kg]
      double mu0 = pi*4.0e-7;  // permeability of free space [H/m]
      double ep0 = 8.8542e-12; // permittivity of free space [F/m]
      double cvac= 2.9979e8;   // speed of light [m/s]

      //   calculate derived parameter scales
      //
      double Pscale  = 2.0*(Nscale)*Tscale*qe;     // pressure scale [J/m^3]
      double Bscale  = pow(mu0*Pscale,0.5);        // magnetic field scale [T]
      double Mi      = Amass*amu;                  // ion mass [kg]
      double Vscale  = pow(Pscale/Mi/Nscale,0.5);  // velocity scale [m/s]
      double Ezscale = Vscale*Bscale;              // electric field scale [V/m]
      double tscale  = Xscale/Vscale;              // time scale [s]
      double etascale  = Xscale*Xscale*mu0/tscale; // resistivity scale [Ohm-m]

      if(procID==0) {
         cout << endl;
         cout << "derived scales:" << endl;
         cout << "velocity scale [m/s] = " << Vscale << endl; 
         cout << "electric field scale [V/m] = " << Ezscale << endl; 
         cout << "magnetic field scale [T] = " << Bscale << endl; 
         //cout << "pressure scale [J/m^3] = " << Pscale << endl; 
         cout << "time scale [s] = " << tscale << endl; 
         cout << "resistivity scale [Ohm-m] = " << etascale << endl; 
      }


      //   calculate dimensionless parameters
      //
      double wpe0 = 5.64e4*pow(Nscale/1.0e6,0.5); // ele plasma freq [rad/s]
      double wpi0 = wpe0*pow(me/Mi,0.5);    // ion plasma freq [rad/s]
      double Le0  = cvac/wpe0;              // ele inertial scale [m]
      double Li0  = cvac/wpi0;              // ion inertial scale [m]
      //
      Json::Value etaVal = Phys.get("eta0",defValue);
      if(etaVal == defValue) {
         eta0   = 1.03e-4/10.0*pow(Tscale,1.5)/(Xscale*Xscale*mu0/tscale); // norm res
      } else {
	 eta0 = etaVal.asDouble();
      } 
      // 
      delta0 = pow(Vscale/cvac,2.0)*epsilonRel;
      Le0or0sq = pow(Le0/Xscale,2.0)*meRel; // (Le0/r0)^2
      if(procID==0) {
	 cout << endl;
         cout << "dimensionless parameters:" << endl;
         cout << "normalized resistivity = " << eta0 << endl;
	 if(etaVal != defValue) cout << "WARNING: USING eta0 FROM INPUT FILE !!!" << endl;
         cout << "(Le0/r0)^2 = " << Le0or0sq << " (Ez relaxation const)"<< endl;      
         cout << "(V0/c)^2 = " << delta0 << " (Jz relaxation const)" << endl << endl;      
      }


      if(advScheme == defValue || gammaVal == defValue ||
	 NsubVal == defValue) {
         cout << "ERROR: advScheme or gamma " << endl;
         cout << "or Nsub " << endl;
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
      Pthresh = 2.0*Nthresh*Tthresh;
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
      if(procID==0) setXminBoundary(P, 0.0, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(P, 0.0, 1.0);   
      Xgrid.communicate(P);
      T = P/N/2.0;
      Cs = pow(gamma0*P/N,0.5);
   } else {
      cout << "value for Physics variable \"P\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   const Json::Value Bvar = Phys.get("B",defValue);
   if(Bvar.isObject()) { 
      Xgrid.setInitialProfile(B,Bvar);
      if(procID==0) setXminBoundary(B, B0, 0.0);   
      if(procID==numProcs-1) setXmaxBoundary(B, 0.0, 0.0);
      Xgrid.communicate(B);
      Bold = B;
   } else {
      cout << "value for Physics variable \"B\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

   Xgrid.DDX(Jz,B); 
   Xgrid.DDX(Jzcc,B); 
   Xgrid.communicate(Jz);
   Xgrid.communicate(Jzcc);
   Jz0 = Jz;
   Jzold = Jz;
   if(min(T)<0.0) cout << " T IS LESS THAN ZERO " << endl;
   eta = eta0/pow(T,1.5);
   Xgrid.InterpToCellEdges(VBce,V*B,B,"C2");
   Xgrid.InterpToCellEdges(etace,eta,eta,"C2");
   Ez = etace*Jz-VBce; 
   Ezold = Ez;

   E = 0.5*M*M/N + P/(gamma0-1.0) + B*B/2.0;
   Eold = E;

   // need inital flux calculation before add to output
   // and before first call to Physics.advance   
   computeFluxes(Xgrid, 1); // 1 for first order calculation   
  

   // add stuff to output files
   //
   dataFile.add(N, "N", 1);      // density 
   dataFile.add(M, "M", 1);      // momentum density 
   dataFile.add(B, "B", 1);      // magnetic field
   dataFile.add(E, "E", 1);      // total energy
   dataFile.add(P, "P", 1);      // total pressure
   dataFile.add(T, "T", 1);      // temperature
   dataFile.add(eta, "eta", 1);  // resistivity
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
   dataFile.add(FluxE, "FluxE", 1);  
   dataFile.add(FluxB, "FluxB", 1);  
   dataFile.add(FluxEz, "FluxEz", 1);  
   //
   dataFile.add(FluxRatio, "FluxRatio", 1);  
   dataFile.add(FluxLim, "FluxLim", 1);  
   dataFile.add(FluxR, "FluxR", 1);
   dataFile.add(FluxL, "FluxL", 1);
   //
   dataFile.add(Nscale,"Nscale",0); 
   dataFile.add(Tscale,"Tscale",0); 
   dataFile.add(Xscale,"Xscale",0); 

} // end Physics.initilize


void Physics::advance(const domainGrid& Xgrid, const double dt)
{
   const int nMax = N.size();
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

      for (auto i=nXg; i<nMax-nXg; i++) {
	 
	 N.at(i) = Nold.at(i) - thisdt*(FluxN.at(i)-FluxN.at(i-1))/Xgrid.dX;
         M.at(i) = Mold.at(i) - thisdt*(FluxM.at(i)-FluxM.at(i-1))/Xgrid.dX
		 - thisdt*Jzcc.at(i)*B.at(i);
         E.at(i) = Eold.at(i) - thisdt*(FluxE.at(i)-FluxE.at(i-1))/Xgrid.dX
		 + thisdt*Jzcc.at(i)*Ezcc.at(i);
	 B.at(i) = Bold.at(i) - thisdt*(FluxB.at(i)-FluxB.at(i-1))/Xgrid.dX;
	 
	 if(N.at(i)<=Nthresh) N.at(i) = Nthresh;
	 if(N.at(i)!=N.at(i)) cout << "bout to go bad: N.at(i) = " << N.at(i) << endl;
	 Ethresh = Pthresh/(gamma0-1.0) + 0.5*M.at(i)*M.at(i)/N.at(i);
	 if(E.at(i)<=Ethresh) E.at(i) =Ethresh;

      }

      timeDomain* tmesh = timeDomain::tmesh;
      //cout << "thist = " << tmesh->tSim << endl;
      double thist = tmesh->tSim;
      B0 = thist*3000.0;
      if(B0>25.0) B0 = 25.0;
      
      if(procID==0) {
         //setXminBoundary(N, N.at(2), 0.0);   
         //setXminBoundary(M, -M.at(2), 0.0);   
         //setXminBoundary(E, E.at(2), 0.0);
         setXminBoundary(N, 0.0, 1.0);   
         setXminBoundary(M, 0.0, -1.0);   
         setXminBoundary(E, 0.0, 1.0);
         //setXminExtrap(N);
         //setXminExtrap(M);
         //setXminExtrap(E);
         setXminBoundary(B, 0.0, 1.0);   
	 if(N.at(0)<Nthresh || N.at(1)<Nthresh) {
            N.at(0) = Nthresh;
            N.at(1) = Nthresh;
            //cout << "Are we here?" << endl;
	 }
      }
     
      if(procID==numProcs-1) {
         setXmaxBoundary(N, 0.0, 1.0);   
         setXmaxBoundary(M, 0.0, -1.0);   
         setXmaxBoundary(E, 0.0, 1.0);   
         setXmaxBoundary(B, B0, 0.0);   
      }

      Xgrid.communicate(N);
      Xgrid.communicate(M);
      Xgrid.communicate(E);
      Xgrid.communicate(B);

      // compute fluxes using partialy updated fields at n+1/2 for 2nd order
      // or n+1 for 1st order
      computeFluxes(Xgrid, Nsub);


      // Now update electric field and current density
      //
      double d0, e0;
  
      for (auto i=nXg-1; i<nMax-nXg; i++) {
         d0 = delta0/thisdt;	   
         e0 = Le0or0sq/thisdt/N.at(i);	   
         //
         Ez.at(i) = (B.at(i+1)-B.at(i))/Xgrid.dX + d0*Ezold.at(i)
	          - (e0*Jzold.at(i) + VBce.at(i))/(e0 + etace.at(i));
         Ez.at(i) /= d0 + 1.0/(e0 + etace.at(i)); 
         //
         Jz.at(i) = e0*Jzold.at(i) + Ez.at(i) + VBce.at(i);
         Jz.at(i) /= e0 + etace.at(i);
      }

      /*
      //cout << "Xce(nXg-1) =" << Xgrid.Xce.at(nXg-1) << endl;
      //cout << "Xce(nMax-Xg-1) =" << Xgrid.Xce.at(nMax-nXg-1) << endl;
      if(procID==0) {
         setXminBoundary(Ez, Ez.at(1), 0.0);   
      }
      if(procID==numProcs-1) {
         setXmaxBoundary(Ez, 0.0, 1.0);   
         // setXmaxExtrap(Ez);   
      }
      */
      Xgrid.communicate(Ez);
      Xgrid.communicate(Jz);
   
   } // finish subcycle steps

   // update old fields
   Nold = N;
   Mold = M;
   Eold = E;
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
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   vector<double> Cspeed, FluxNcc, FluxMcc, FluxEcc;
   vector<double> Cspeed2, FluxB0cc, Eprime;
   vector<double> FluxVisc, dVdx; 
   FluxNcc.assign(nCC,0.0);
   FluxMcc.assign(nCC,0.0);
   FluxEcc.assign(nCC,0.0);
   FluxB0cc.assign(nCC,0.0);
   Eprime.assign(nCE,0.0);
   FluxVisc.assign(nCC,0.0);
   dVdx.assign(nCC,0.0);

   //  define derived variables
   //
   V  = M/N;
   if(min(N)<0.0) cout << " N IS LESS THAN ZERO " << endl;
   P  = (E - 0.5*V*M)*(gamma0-1.0);
   T  = P/2.0/N;
   if(min(T)<0.0) cout << " T IS LESS THAN ZERO " << endl;
   Cs = sqrt(gamma0*P/N + B*B/N);
   //eta = eta0/T/sqrt(T);
   eta = eta0/T/sqrt(T)*(1.0+1000.0*pow(0.01/N,4));
   Xgrid.InterpToCellEdges(etace,eta,eta,"C2");
   Xgrid.communicate(etace);
   Xgrid.communicate(eta);


   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeed  = abs(V) + Cs; // adv flux jacobian
   FluxNcc = M;
   FluxMcc = M*V + P;
   FluxEcc = (0.5*V*M + P*gamma0/(gamma0-1.0) )*V;
   FluxB0cc = V*B;
   FluxEz = -B;
   
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
                           FluxNcc,Cspeed,N,Nsub);
      Xgrid.computeFluxTVD(FluxM,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxMcc,Cspeed,M,Nsub);
      Xgrid.computeFluxTVD(FluxE,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxEcc,Cspeed,E,Nsub);
      Xgrid.computeFluxTVD(VBce,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxB0cc,Cspeed,B,Nsub);
   }
   else {
      Xgrid.InterpToCellEdges(FluxN,FluxNcc,V,advScheme0);
      Xgrid.InterpToCellEdges(FluxM,FluxMcc,V,advScheme0);
      Xgrid.InterpToCellEdges(FluxE,FluxEcc,V,advScheme0);
      //Xgrid.InterpToCellEdges(FluxB,FluxBcc,V,advScheme0);
      //Xgrid.InterpToCellEdges(FluxEz,FluxEzcc,Ez,advScheme0);
   } 
   //Xgrid.InterpToCellEdges(VBce,FluxB0cc,B,"C2");
   //FluxM = FluxM + FluxVisc;

   if(procID==0) {
      setXminBoundary(FluxN, 0.0, 0.0);   
      setXminBoundary(FluxM, (P.at(2)+P.at(1))/2.0, 0.0);   
      setXminBoundary(FluxE, 0.0, 0.0);
      setXminBoundary(VBce, 0.0, 0.0);
   }
   if(procID==numProcs-1) {
      setXmaxBoundary(FluxN, 0.0, 0.0);
      setXmaxBoundary(FluxE, 0.0, 0.0);
      const int thisnX = P.size();
      double P0 = (P.at(thisnX-3)+P.at(thisnX-2))/2.0;
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
   Xgrid.DDX(Jz0,B); 
   Xgrid.communicate(Jz0);

   Xgrid.communicate(FluxN);   
   Xgrid.communicate(FluxM);   
   Xgrid.communicate(FluxE);   
   
   Xgrid.communicate(FluxL);   
   Xgrid.communicate(FluxR);   
   Xgrid.communicate(FluxRatio);   
   Xgrid.communicate(FluxLim);   

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


void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid)
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
   if(procID==0) {
      cout << "sigma_0*dt/delta = " << dtSim/delta0/eta0 << endl;
      cout << "dtSim = " << dtSim << endl;
   }
}
