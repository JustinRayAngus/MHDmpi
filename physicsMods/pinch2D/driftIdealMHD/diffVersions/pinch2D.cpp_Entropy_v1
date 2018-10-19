/***
 * 
 * physics module for 2D zpinch with m=0 using drift-MHD 
 * with electron inertia and displacment current (relaxation)
 * and with electron pressure and with seperate electron and
 * ion heat equations with diamagnetic heat fluxes and thermalization
 *
 * dN/dt  + d(r*Mr)/r/dr + d(Mz)/dz = 0
 * dMr/dt + d(r*(Mr*Ur + P))/r/dr + d(Mr*Uz)/dz = -Jz*By + P/r
 * dMz/dt + d(r*Mz*Ur)/r/dr  + d(Mz*Uz + P)/dz = Jr*By
 * dSe/dt + d(r*Se*Uer)/r/dr + d(Se*Uez)/dz = (gamma0-1)*(-divqe + Qe)/N^(gamma0-1)
 * dSi/dt + d(r*Si*Ur)/r/dr  + d(Si*Uz)/dz  = (gamma0-1)*(-divqi + Qi)/N^(gamma0-1)  
 * dBy/dt + d(-Ez)/dr + d(Er)/dz = 0
 *
 * dEr/dt = 1/delta0*(Jr0 - Jr)
 * dEz/dt = 1/delta0*(Jz0 - Jz)
 *
 * dJr/dt = 1/epsilon0*[N*(Er - Er0) + lambda0*(Jz*By + dPe/dr)]
 * dJz/dt = 1/epsilon0*[N*(Ez - Ez0) - lambda0*(Jr*By - dPe/dz)]
 *
 * Ez0 = -Ur*By
 * Er0 = Uz*By
 * Jz0 = d(r*By)/r/dr
 * Jr0 = -d(By)/dz
 *
 * Pe = Se*N^(gamma0-1)
 * Pi = Si*N^(gamma0-1)
 * Te = Pe/N
 * Ti = Pi/N
 * P  = Pe + Pi
 *
 * qe = -gamma0/(gamma0-1)*lambda0*Pe*B x nabla Te / B^2
 * qi =  gamma0/(gamma0-1)*lambda0*Pi*B x nabla Ti / B^2
 *
 * Qi = 2/(gamma0-1)*nuT*N*(Te-Ti)
 * Qe = -Qi
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
double lambda0, epsilon0, delta0, taui0;
double nuTherm0, TempRatio0;
int Nsub, modelGyroVisc;
matrix2D<double> N, Mx, Mz, Se, Si, By, Ez, Ex, Jz, Jx;  // time-evolving variables
matrix2D<double> P0, deltaP;                    // initial perturbation
matrix2D<double> eta, Cs, Vx, Vz, P, T, S;     // derived variables
matrix2D<double> Ez0, Jz0, Ex0, Jx0;        // derived variables
matrix2D<double> eta_x, eta_z, Jzcc, Ezcc, Jxcc, VxBy_x, VzBy_z;
matrix2D<double> Nold, Mxold, Mzold, Seold, Siold, Byold;
matrix2D<double> Exold, Jxold, Ezold, Jzold;
matrix2D<double> FluxR_x, FluxL_x, FluxRatio_x, FluxLim_x;
matrix2D<double> FluxR_z, FluxL_z, FluxRatio_z, FluxLim_z;  
matrix2D<double> FluxN_x, FluxSe_x, FluxSi_x, FluxMx_x, FluxMz_x, FluxBy_x, FluxEz_x;
matrix2D<double> FluxN_z, FluxSe_z, FluxSi_z, FluxMx_z, FluxMz_z, FluxBy_z, FluxEx_z;
matrix2D<double> FluxEz_x2, FluxEx_z2, Jx02;
matrix2D<double> Fx, rcc, rce_z, rce_x;
//
matrix2D<double> Ve_drift_x, Ve_drift_z, Vi_drift_x, Vi_drift_z;
matrix2D<double> Pe, Pi, Te, Ti, divqe, divqi, Qe, Qi;
matrix2D<double> Vex, Vez, Vhallx, Vhallz;
matrix2D<double> qe_x, qe_z, qi_x, qi_z;
matrix2D<double> Fluxqe_x, Fluxqe_z, Fluxqi_x, Fluxqi_z;
//
matrix2D<double> eta0, eta1, eta3;
matrix2D<double> Qvis, divPvis_x, divPvis_z;


// Set lower threshold values for N, T, and S
double Nthresh=1.0e-4, Tthresh=2.5e-2, Sthresh;

void computeFluxes(const domainGrid&, const int);
void computeFluxes_E(const domainGrid&, const int);
void computeViscousStress(const domainGrid&);
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
   Se.initialize(nXcc,nZcc,0.0);
   Seold.initialize(nXcc,nZcc,0.0);
   Si.initialize(nXcc,nZcc,0.0);
   Siold.initialize(nXcc,nZcc,0.0);
   By.initialize(nXcc,nZcc,0.0);
   Byold.initialize(nXcc,nZcc,0.0);
   P.initialize(nXcc,nZcc,0.0);
   T.initialize(nXcc,nZcc,0.0);
   S.initialize(nXcc,nZcc,0.0);
   Cs.initialize(nXcc,nZcc,0.0);
   Vx.initialize(nXcc,nZcc,0.0);
   Vz.initialize(nXcc,nZcc,0.0);


   // extra stuff for drift physics
   //
   Ve_drift_x.initialize(nXcc,nZcc,0.0);
   Ve_drift_z.initialize(nXcc,nZcc,0.0);
   Vi_drift_x.initialize(nXcc,nZcc,0.0);
   Vi_drift_z.initialize(nXcc,nZcc,0.0);
   Pe.initialize(nXcc,nZcc,0.0);
   Te.initialize(nXcc,nZcc,0.0);
   Pi.initialize(nXcc,nZcc,0.0);
   Ti.initialize(nXcc,nZcc,0.0);
   Vex.initialize(nXcc,nZcc,0.0);
   Vez.initialize(nXcc,nZcc,0.0);
   Vhallx.initialize(nXcc,nZcc,0.0);
   Vhallz.initialize(nXcc,nZcc,0.0);
   divqe.initialize(nXcc,nZcc,0.0);   
   divqi.initialize(nXcc,nZcc,0.0);   
   Qe.initialize(nXcc,nZcc,0.0);   
   Qi.initialize(nXcc,nZcc,0.0);   
   qe_x.initialize(nXcc,nZcc,0.0);
   qe_z.initialize(nXcc,nZcc,0.0);
   qi_x.initialize(nXcc,nZcc,0.0);
   qi_z.initialize(nXcc,nZcc,0.0);
   
   eta0.initialize(nXcc,nZcc,0.0);
   eta1.initialize(nXcc,nZcc,0.0);
   eta3.initialize(nXcc,nZcc,0.0);
   Qvis.initialize(nXcc,nZcc,0.0);
   divPvis_x.initialize(nXcc,nZcc,0.0);
   divPvis_z.initialize(nXcc,nZcc,0.0);

   Jz.initialize(nXcc,nZcc,0.0);
   Ez.initialize(nXcc,nZcc,0.0);
   Jz0.initialize(nXcc,nZcc,0.0);
   Ez0.initialize(nXcc,nZcc,0.0);
   Jzold.initialize(nXcc,nZcc,0.0);
   Ezold.initialize(nXcc,nZcc,0.0);
   //
   Jx.initialize(nXcc,nZcc,0.0); 
   Ex.initialize(nXcc,nZcc,0.0);
   Jx0.initialize(nXcc,nZcc,0.0); 
   //Jx02.initialize(nXcc,nZcc,0.0); 
   Ex0.initialize(nXcc,nZcc,0.0);
   Jxold.initialize(nXcc,nZcc,0.0); 
   Exold.initialize(nXcc,nZcc,0.0);

   //
   FluxN_x.initialize(nXce, nZcc, 0.0);
   FluxMx_x.initialize(nXce,nZcc, 0.0);
   FluxMz_x.initialize(nXce,nZcc, 0.0);
   FluxSe_x.initialize(nXce, nZcc, 0.0);
   FluxSi_x.initialize(nXce, nZcc, 0.0);
   FluxBy_x.initialize(nXce,nZcc, 0.0);
   FluxEz_x.initialize(nXce,nZcc, 0.0);
   FluxEz_x2.initialize(nXce,nZcc, 0.0);
   Fluxqe_x.initialize(nXce,nZcc, 0.0);
   Fluxqi_x.initialize(nXce,nZcc, 0.0);
   
   FluxN_z.initialize(nXcc, nZce, 0.0);
   FluxMx_z.initialize(nXcc,nZce, 0.0);
   FluxMz_z.initialize(nXcc,nZce, 0.0);
   FluxSe_z.initialize(nXcc, nZce, 0.0);
   FluxSi_z.initialize(nXcc, nZce, 0.0);
   FluxBy_z.initialize(nXcc,nZce, 0.0);
   FluxEx_z.initialize(nXcc,nZce, 0.0);
   FluxEx_z2.initialize(nXcc,nZce, 0.0);
   Fluxqe_z.initialize(nXcc,nZce, 0.0);
   Fluxqi_z.initialize(nXcc,nZce, 0.0);
 
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
      Json::Value advScheme  = Phys.get("advScheme",defValue);
      Json::Value gammaVal   = Phys.get("gammaC",defValue);
      Json::Value lambdaVal  = Phys.get("lambda",defValue);
      Json::Value epsilonVal = Phys.get("epsilon",defValue);
      Json::Value deltaVal   = Phys.get("delta",defValue);
      Json::Value tauiVal    = Phys.get("taui",defValue);
      Json::Value modelGyroViscVal    = Phys.get("modelGyroVisc",defValue);
      Json::Value nuThermVal = Phys.get("nuTherm",defValue);
      Json::Value TempRatioVal  = Phys.get("TempRatio",defValue);
      Json::Value NsubVal    = Phys.get("Nsub",defValue);
      if(advScheme == defValue || gammaVal == defValue ||
	 NsubVal == defValue ) {
         cout << "ERROR: advScheme or gamma " << endl;
         cout << "or Nsub is " << endl;
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

      nuTherm0 = nuThermVal.asDouble();
      if(procID==0) cout << "thermalization coefficent = " << nuTherm0 << endl;
      if(nuTherm0 < 0.0) {
         printf("ERROR: thermalization coefficient can't be < 0 \n");
         exit (EXIT_FAILURE);
      }
      
      TempRatio0 = TempRatioVal.asDouble();
      if(procID==0) cout << "Ti0/Te0 = " << TempRatio0 << endl;
     
 
      lambda0 = lambdaVal.asDouble();
      if(procID==0) cout << "Li/L0 = " << lambda0 << endl;
      if(lambda0 < 0.0) {
         printf("ERROR: ion skin depth must be positive \n");
         exit (EXIT_FAILURE);
      }
      
      epsilon0 = epsilonVal.asDouble();
      if(procID==0) cout << "(Le/L0)^2 = " << epsilon0 << endl;
      if(epsilon0 > 1.0 || epsilon0 < 0.0) {
         printf("ERROR: electron skin depth should not be larger than L0 \n");
         exit (EXIT_FAILURE);
      }
      
      delta0 = deltaVal.asDouble();
      if(procID==0) cout << "(V0/Clight)^2 = " << delta0 << endl;
      if(delta0 > 1.0) {
         printf("ERROR: V0 cannot be larger than speed of light \n");
         exit (EXIT_FAILURE);
      }
      
      taui0 = tauiVal.asDouble();
      if(procID==0) cout << "taui0/t0 = " << taui0 << endl;
      
      Nsub = NsubVal.asInt();
      if(procID==0) cout << "Nsub = " << Nsub << endl;
      if(gamma0 < 1.0) {
         printf("ERROR: Nsub must be int >= 1\n");
         exit (EXIT_FAILURE);
      }
      Sthresh = Nthresh*Tthresh/pow(Nthresh,gamma0-1);
      
      modelGyroVisc = modelGyroViscVal.asInt();
      if(procID==0) {
         if(modelGyroVisc==0) { 
            cout << "Not modeling gyro viscosity"<< endl;
         }
         else if(modelGyroVisc==1) { 
            cout << "modeling gyro viscosity"<< endl;
         }
         else{   
            printf("ERROR: modelGyroVisc should be 0 or 1 \n");
         }
      }

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
      cout << "value for Physics variable \"P\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   
   const Json::Value Vzvar = Phys.get("Vz",defValue);
   if(Vzvar.isObject()) { 
      Xgrid.setInitialProfile(Vz,Vzvar);
      if(procID==0) setXminExtrap(Vz, 0);   
      if(procID==numProcs-1) setXmaxExtrap(Vz,0); 
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
   rce_z.initialize(nXcc,nZce,0.0);
   rce_x.initialize(nXce,nZcc,0.0);
   for (auto i=0; i<nXcc; i++) {
      for (auto j=0; j<nZcc; j++) {
	 rcc(i,j) = Xgrid.Xcc.at(i);
	 if(i<nXce) rce_x(i,j) = Xgrid.Xce.at(i);
	 if(j<nZce) rce_z(i,j) = Xgrid.Xcc.at(i);
      }
   }
   
   // calculate current density 
   //
   Xgrid.DDX(Jz,rcc*By); 
   Jz /= rcc;
   setXmaxBoundary(Jz,0.0,0.0);
   Xgrid.communicate(Jz);
   Jz0 = Jz;
   Jzold = Jz;

   // add perturbation to pressure profile
   //
   P0 = P;
   P *= 1.0+deltaP;
   N = P/2.0/T;
   Nold = N;
   Pe = P/(1.0+TempRatio0);
   Pi = P-Pe;
   Se = Pe/pow(N,gamma0-1.0);
   Si = Pi/pow(N,gamma0-1.0);
   Seold  = Se;
   Siold  = Si;
   S = Se + Si;
   Cs = pow(gamma0*P/N + By*By/N,0.5);
   
   // set initial electric field from Hall term and pressure
   //
   Te = Pe/N;
   Ti = Pi/N;
   T = (Te + Ti)/2.0;
   matrix2D<double> dPedx,dPedz;
   dPedx.initialize(nXcc, nZcc, 0.0);
   dPedz.initialize(nXcc, nZcc, 0.0);
   Xgrid.DDX(dPedx,Pe);
   Xgrid.DDZ(dPedz,Pe);

   //Ez0 = -Vx*By;
   Ex0 =  Vz*By;
   Ex = Ex0 - lambda0/N*(By*Jz + dPedx);
   //Ez = Ez0 + lambda0/N*(By*Jx - dPedz);
   setXminExtrap(Ez,0);
   setXmaxBoundary(Ez, 0.0, -1.0);   
   setZboundaryPeriodic(Ex);  
   setZboundaryPeriodic(Ez);  
   Xgrid.communicate(Ez);
   Xgrid.communicate(Ex);
  
   // set initial electron velocity and ion momentum
   //
   Vhallz = 0.0 - lambda0/N*Jz;
   //Vhallx = 0.0 - lambda0/N*Jx;
   Vez = Vz + Vhallz;
   //Vex = Vx + Vhallx;
   Mz = N*Vz;
   Mzold = Mz;


   computeFluxes(Xgrid, 2);  


   // add stuff to output files
   //
   dataFile.add(rcc, "rcc", 0);  // force density x-direction 
   dataFile.add(Fx, "Fx", 1);  // force density x-direction 
   dataFile.add(N,  "N",  1);  // density 
   dataFile.add(deltaP,  "deltaP",  0);  // density perturbation 
   dataFile.add(P0, "P0",  0);  // initial unperturbed pressure profile 
   dataFile.add(Mx, "Mx", 1);  // momentum density 
   dataFile.add(Mz, "Mz", 1);  // momentum density 
   dataFile.add(S,  "S",  1);  // entropy density
   dataFile.add(Se, "Se", 1);  // entropy density
   dataFile.add(Si, "Si", 1);  // entropy density
   dataFile.add(By, "By", 1);  // magnetic field
   dataFile.add(P,  "P",  1);  // pressure
   dataFile.add(Pe, "Pe", 1);  // pressure
   dataFile.add(Pi, "Pi", 1);  // pressure
   dataFile.add(T,  "T",  1);  // temperature
   dataFile.add(Te, "Te", 1);  // temperature
   dataFile.add(Ti, "Ti", 1);  // temperature
   dataFile.add(Vx, "Vx", 1);  // x-velocity
   dataFile.add(Vz, "Vz", 1);  // z-velocity
   dataFile.add(Jz, "Jz", 1);  // current density
   dataFile.add(Ez, "Ez", 1);  // z-electric field
   dataFile.add(Jx, "Jx", 1);  // current density
   dataFile.add(Ex, "Ex", 1);  // x-electric field
   dataFile.add(Jz0, "Jz0", 1);  // current density
   dataFile.add(Ez0, "Ez0", 1);  // z-electric field
   dataFile.add(Jx0, "Jx0", 1);  // current density
   //dataFile.add(Jx02, "Jx02", 1);  // current density
   dataFile.add(Ex0, "Ex0", 1);  // x-electric field
   dataFile.add(Cs, "Cs", 1);  // sound speed
   dataFile.add(gamma0,"gamma0",0); 
   dataFile.add(lambda0,"lambda0",0);   // Li/L0 
   dataFile.add(epsilon0,"epsilon0",0); // (Le/L0)^2
   dataFile.add(delta0,"delta0",0);     // (V0/C)^2
   dataFile.add(nuTherm0,"nuTherm0",0); 
   //
   dataFile.add(FluxN_x,  "FluxN_x",  1);  
   dataFile.add(FluxMx_x, "FluxMx_x", 1);  
   dataFile.add(FluxMz_x, "FluxMz_x", 1);  
   dataFile.add(FluxSe_x, "FluxSe_x", 1);  
   dataFile.add(FluxBy_x, "FluxBy_x", 1);  
   dataFile.add(FluxEz_x, "FluxEz_x", 1);  
   dataFile.add(FluxEx_z, "FluxEx_z", 1);  
   dataFile.add(FluxEz_x2, "FluxEz_x2", 1);  
   dataFile.add(FluxEx_z2, "FluxEx_z2", 1);  
   //
   dataFile.add(FluxRatio_x, "FluxRatio_x", 1);  
   dataFile.add(FluxLim_x, "FluxLim_x", 1);  
   dataFile.add(FluxR_x, "FluxR_x", 1);
   dataFile.add(FluxL_x, "FluxL_x", 1);
   //
   dataFile.add(divqe,"divqe",1);
   dataFile.add(divqi,"divqi",1);
   dataFile.add(Qe,"Qe",1);
   dataFile.add(Qi,"Qi",1);
   //
   dataFile.add(Qvis,"Qvis",1);  

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
   
   matrix2D<double> dPedx, dPedz, Nhalf;
   dPedx.initialize(nXcc,nZcc,0.0);
   dPedz.initialize(nXcc,nZcc,0.0);
   Nhalf.initialize(nXcc,nZcc,0.0);

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   double a, b, c, d, e, f, detA, Az, Ax;

   //domainGrid* mesh = domainGrid::mesh;

   // Explicit forward advance from n to n+1 using subcycling in time 
   //  
   Nhalf = N;
   double thisdt;
   for (auto n=1; n<Nsub+1; n++) {
      thisdt = dt*n/Nsub;

      // Update N, M, S, and B terms from n to n+1
      for (auto i=nXg; i<nXcc-nXg; i++) {
         for (auto j=nZg; j<nZcc-nZg; j++) {

         Fx(i,j) = - (FluxMx_x(i,j)-FluxMx_x(i-1,j))/Xcc.at(i)/Xgrid.dX
		   + P(i,j)/Xcc.at(i) - Jz(i,j)*By(i,j);
	 if(procID==numProcs-1 && i==nXcc-nXg-1) Fx(i,j) = Fx(i-1,j)/3.0;
	 //if(i==nXcc-nXg-1) Fx(i,j) = 2.0*Fx(i-1,j)-Fx(i-2,j); // doesn't work well

	 N(i,j)  = Nold(i,j)  - thisdt*(FluxN_x(i,j)  - FluxN_x(i-1,j))/Xcc.at(i)/Xgrid.dX
	                      - thisdt*(FluxN_z(i,j)  - FluxN_z(i,j-1))/Xgrid.dZ;
         Mx(i,j) = Mxold(i,j) + thisdt*Fx(i,j)
                              - thisdt*divPvis_x(i,j)
                              - thisdt*(FluxMx_z(i,j) - FluxMx_z(i,j-1))/Xgrid.dZ;
         Mz(i,j) = Mzold(i,j) - thisdt*(FluxMz_x(i,j) - FluxMz_x(i-1,j))/Xcc.at(i)/Xgrid.dX
                              - thisdt*(FluxMz_z(i,j) - FluxMz_z(i,j-1))/Xgrid.dZ
                              - thisdt*divPvis_z(i,j)
                              + thisdt*Jx(i,j)*By(i,j);
         Se(i,j) = Seold(i,j) - thisdt*(FluxSe_x(i,j) - FluxSe_x(i-1,j))/Xcc.at(i)/Xgrid.dX
                              - thisdt*(FluxSe_z(i,j) - FluxSe_z(i,j-1))/Xgrid.dZ
                              + thisdt*(-divqe(i,j) + Qe(i,j))
                                      *(gamma0-1.0)/pow(Nhalf(i,j),gamma0-1);
         Si(i,j) = Siold(i,j) - thisdt*(FluxSi_x(i,j) - FluxSi_x(i-1,j))/Xcc.at(i)/Xgrid.dX
                              - thisdt*(FluxSi_z(i,j) - FluxSi_z(i,j-1))/Xgrid.dZ
                              + thisdt*(-divqi(i,j) + Qi(i,j) + Qvis(i,j))
                                      *(gamma0-1.0)/pow(Nhalf(i,j),gamma0-1.0);
	 By(i,j) = Byold(i,j) - thisdt*(FluxBy_x(i,j) - FluxBy_x(i-1,j))/Xgrid.dX
	                      - thisdt*(FluxBy_z(i,j) - FluxBy_z(i,j-1))/Xgrid.dZ;
	 
	 //if(N(i,j)<=Nthresh) N(i,j) = Nthresh;
	 //if(M(i,j)<0.0) M(i,j) = 0.0;
	 //if(S(i,j)<=Sthresh) S(i,j) = Sthresh;

         }
      }

      timeDomain* tmesh = timeDomain::tmesh;
      //cout << "thist = " << tmesh->tSim << endl;
      //double thist = tmesh->tSim;
      
      if(procID==0) {
         setXminExtrap(N, 0);   
         setXminBoundary(Mx,0.0,-1.0);   
         setXminExtrap(Mz,0);   
         //setXminBoundary(Mz,0.0,1.0);   
         setXminExtrap(Se,0);
         setXminExtrap(Si,0);
         setXminBoundary(By,0.0,-1.0);   
      }
     
      if(procID==numProcs-1) {
         setXmaxBoundary(Mx, 0.0, -1.0);   
         setXmaxExtrap(N, 0);   
         setXmaxExtrap(Mz,0);   
         setXmaxExtrap(Se,0);   
         setXmaxExtrap(Si,0);   
         setXmaxBy(By,-1);
      }
      setZboundaryPeriodic(N);
      setZboundaryPeriodic(Mx);
      setZboundaryPeriodic(Mz);
      setZboundaryPeriodic(Se);
      setZboundaryPeriodic(Si);
      setZboundaryPeriodic(By);
      //
      Xgrid.communicate(Fx);
      Xgrid.communicate(N);
      Xgrid.communicate(Mx);
      Xgrid.communicate(Mz);
      Xgrid.communicate(Se);
      Xgrid.communicate(Si);
      Xgrid.communicate(By);


      // Now update E and J, which are done implicitly
      //
      computeFluxes_E(Xgrid, 2); // compute electric field fluxes 
      /*
      Xgrid.DDX(Jz0,rcc*By); 
      Jz0 /= rcc;
      setXmaxBoundary(Jz0,0.0,0.0);
      Xgrid.communicate(Jz0);
   
      Xgrid.DDZ(Jx0,By);
      Xgrid.communicate(Jx0);
      Jx0 *= -1.0; 
      */

      Ez0 = -Vx*By;
      Ex0 = Vz*By;
      Xgrid.DDX(dPedx,Pe); // electron pressure gradient 
      Xgrid.DDZ(dPedz,Pe); // electron pressure gradient 
      for (auto i=nXg; i<nXcc-nXg; i++) {
         for (auto j=nZg; j<nZcc-nZg; j++) {
	 
         Jz0(i,j) = -(FluxEz_x(i,j) - FluxEz_x(i-1,j))/Xgrid.dX/Xcc.at(i);
         Jx0(i,j) = -(FluxEx_z(i,j) - FluxEx_z(i,j-1))/Xgrid.dZ;
         //Jx02(i,j) = -(FluxEx_z2(i,j) - FluxEx_z2(i,j-1))/Xgrid.dZ;

         // formulate matrix elements for J update
         // a*Jx^n+1 + b*Jz^n+1 = e
         // c*Jx^n+1 + d*Jz^n+1 = f
         //
         a = 1.0+thisdt*thisdt*N(i,j)/epsilon0/delta0;
         d = a;
         b = -thisdt/epsilon0*lambda0*By(i,j);
         c = -b;
         detA = a*d-b*c;
         e = Jxold(i,j) + thisdt*N(i,j)/epsilon0*(Exold(i,j)-Ex0(i,j))
           + thisdt/epsilon0*lambda0*dPedx(i,j)
           + thisdt*thisdt*N(i,j)/epsilon0/delta0*Jx0(i,j);
         f = Jzold(i,j) + thisdt*N(i,j)/epsilon0*(Ezold(i,j)-Ez0(i,j))
           + thisdt/epsilon0*lambda0*dPedz(i,j)
           + thisdt*thisdt*N(i,j)/epsilon0/delta0*Jz0(i,j);
         //
         Jx(i,j) = (e*d - b*f)/detA;
         Jz(i,j) = (f*a - e*c)/detA;
      
         // use new J values to update E
         //
         Ex(i,j) = Exold(i,j) + thisdt/delta0*(Jx0(i,j) - Jx(i,j));
         Ez(i,j) = Ezold(i,j) + thisdt/delta0*(Jz0(i,j) - Jz(i,j));
         

         // neglect electron inertia => NO GOOD, ~1/lambda0/By
         //
         //Jx(i,j) =  N(i,j)/lambda0*(Ez(i,j)-Ez0(i,j))/By(i,j) + dPedx(i,j)/By(i,j);
         //Jz(i,j) = -N(i,j)/lambda0*(Ex(i,j)-Ex0(i,j))/By(i,j) - dPedz(i,j)/By(i,j);
 
         /*
         Jx(i,j) = Jx0(i,j);
         Jz(i,j) = Jz0(i,j);

         c = N(i,j)*thisdt/(delta0*By(i,j)*lambda0);
         Ax = Exold(i,j) + thisdt/delta0*(Jx(i,j) - dPedz(i,j)/By(i,j) 
                                        + Ez0(i,j)/N(i,j)/lambda0/By(i,j));
         Az = Ezold(i,j) + thisdt/delta0*(Jz(i,j) + dPedx(i,j)/By(i,j) 
                                        - Ex0(i,j)/N(i,j)/lambda0/By(i,j));
   
         // update E
         //
         Ex(i,j) = (Ax - c*Az)/(1.+c*c);
         Ez(i,j) = (Az + c*Ax)/(1.+c*c);
         */

         }
      }
      if(procID==0) {
         setXminExtrap(Ex,0);
         setXminExtrap(Ez,0);
         setXminExtrap(Jx,0);
         setXminExtrap(Jz,0);
      }
      if(procID==numProcs-1) {
         setXmaxExtrap(Ex, 0);   
         setXmaxBoundary(Ez, 0.0, -1.0);   
         setXmaxExtrap(Jx,0);
         setXmaxExtrap(Jz,0);
      }
      setZboundaryPeriodic(Ex);
      setZboundaryPeriodic(Ez);
      setZboundaryPeriodic(Jx);
      setZboundaryPeriodic(Jz);
      //
      Xgrid.communicate(Ex);
      Xgrid.communicate(Ez);
      Xgrid.communicate(Jx);
      Xgrid.communicate(Jz);
      
      Xgrid.communicate(Jx0);
      Xgrid.communicate(Jz0);
      //Xgrid.communicate(Jx02);


      // compute fluxes using partialy updated fields at n+1/2 for 2nd order
      // or n+1 for 1st order
      computeFluxes(Xgrid, Nsub); // compute second order fluxes at n+1/2

   } // finish subcycle steps for N, M, S, and B

   // update old fields
   Nold = N;
   Mxold = Mx;
   Mzold = Mz;
   Seold = Se;
   Siold = Si;
   Byold = By;
   Exold = Ex;
   Ezold = Ez;
   Jxold = Jx;
   Jzold = Jz;
   
   // compute fluxes using fully updated fields at n+1
   //computeFluxes(Xgrid, 2);

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
   
   matrix2D<double> FluxNcc_x, FluxMxcc_x, FluxMzcc_x, FluxSecc_x, FluxSicc_x, FluxBycc_x;
   matrix2D<double> FluxNcc_z, FluxMxcc_z, FluxMzcc_z, FluxSecc_z, FluxSicc_z, FluxBycc_z;
   matrix2D<double> FluxSecc_hall_x, FluxSecc_hall_z, FluxSecc_x0, FluxSecc_z0;
   matrix2D<double> Cspeedx, Cspeedz, Cspeede, Cvac; 
 
   Cspeedx.initialize(nXcc,nZcc,0.0);
   Cspeedz.initialize(nXcc,nZcc,0.0);
   Cspeede.initialize(nXcc,nZcc,0.0);
   double Cvac0 = 1.0/sqrt(delta0);
   Cvac.initialize(nXcc,nZcc,Cvac0);
   FluxNcc_x.initialize(nXcc,nZcc,0.0);
   FluxNcc_z.initialize(nXcc,nZcc,0.0);
   FluxMxcc_x.initialize(nXcc,nZcc,0.0);
   FluxMxcc_z.initialize(nXcc,nZcc,0.0);
   FluxMzcc_x.initialize(nXcc,nZcc,0.0);
   FluxMzcc_z.initialize(nXcc,nZcc,0.0);
   FluxSecc_x.initialize(nXcc,nZcc,0.0);
   FluxSecc_x0.initialize(nXcc,nZcc,0.0);
   FluxSecc_hall_x.initialize(nXcc,nZcc,0.0);
   FluxSicc_x.initialize(nXcc,nZcc,0.0);
   FluxSecc_z.initialize(nXcc,nZcc,0.0);
   FluxSecc_z0.initialize(nXcc,nZcc,0.0);
   FluxSecc_hall_z.initialize(nXcc,nZcc,0.0);
   FluxSicc_z.initialize(nXcc,nZcc,0.0);
   FluxBycc_x.initialize(nXcc,nZcc,0.0);
   FluxBycc_z.initialize(nXcc,nZcc,0.0);
   

   matrix2D<double> Ez0, Ex0, Ezprime, Exprime, FluxBycc_x0, FluxBycc_z0;
   matrix2D<double> FluxBycc_xp, FluxBycc_zp, FluxBy_xp, FluxBy_zp;
   matrix2D<double> FluxSe_hall_x, FluxSe_hall_z;

   Ez0.initialize(nXcc,nZcc,0.0);
   Ex0.initialize(nXcc,nZcc,0.0);
   Ezprime.initialize(nXcc,nZcc,0.0);
   Exprime.initialize(nXcc,nZcc,0.0);
   FluxBycc_x0.initialize(nXcc,nZcc,0.0);
   FluxBycc_z0.initialize(nXcc,nZcc,0.0);
   FluxBycc_xp.initialize(nXcc,nZcc,0.0);
   FluxBycc_zp.initialize(nXcc,nZcc,0.0);
   FluxBy_xp.initialize(nXce,nZcc,0.0);
   FluxBy_zp.initialize(nXcc,nZce,0.0);
   FluxSe_hall_x.initialize(nXce,nZcc,0.0);
   FluxSe_hall_z.initialize(nXcc,nZce,0.0);
   

   //  define derived variables
   //
   Vx = Mx/N;
   Vz = Mz/N;
   Pe  = Se*pow(N,gamma0-1.0);
   Pi  = Si*pow(N,gamma0-1.0);
   Te  = Pe/N;
   Ti  = Pi/N;
   S = Se + Si;
   P = Pe + Pi;
   T = (Te + Ti)/2.0;
   Cs = sqrt(gamma0*P/N + By*By/N);

   
   // set drift velocities
   //
   matrix2D<double> dByrdx, dBydz;
   dByrdx.initialize(nXcc,nZcc,0.0);
   dBydz.initialize(nXcc,nZcc,0.0);
   Xgrid.DDZ(dBydz,By);
   Xgrid.DDX(dByrdx,rcc/By);
   Xgrid.communicate(dBydz);
   Xgrid.communicate(dByrdx);

   Ve_drift_x = -lambda0*Te*dBydz/By/By;
   Ve_drift_z = -lambda0*Te*dByrdx/rcc;
   Vi_drift_x =  lambda0*Ti*dBydz/By/By;
   Vi_drift_z =  lambda0*Ti*dByrdx/rcc;
   
   
   //   compute entropy density source terms
   //
   matrix2D<double> dTedx, dTedz, dTidx, dTidz;
   matrix2D<double> dPedx, dPedz, dPidx, dPidz;
   matrix2D<double> Sqe_x, Sqe_z, Sqi_x, Sqi_z;
   matrix2D<double> Fluxqe_x_cc, Fluxqe_z_cc, Fluxqi_x_cc, Fluxqi_z_cc;
   matrix2D<double> omegatau, omegatau2, omegatau4, omegatauFactor;
   dTedx.initialize(nXcc,nZcc,0.0);
   dTedz.initialize(nXcc,nZcc,0.0);
   dTidx.initialize(nXcc,nZcc,0.0);
   dTidz.initialize(nXcc,nZcc,0.0);
   dPedx.initialize(nXcc,nZcc,0.0);
   dPedz.initialize(nXcc,nZcc,0.0);
   dPidx.initialize(nXcc,nZcc,0.0);
   dPidz.initialize(nXcc,nZcc,0.0);
   Sqe_x.initialize(nXcc,nZcc,0.0);
   Sqe_z.initialize(nXcc,nZcc,0.0);
   Sqi_x.initialize(nXcc,nZcc,0.0);
   Sqi_z.initialize(nXcc,nZcc,0.0);
   Fluxqe_x_cc.initialize(nXcc,nZcc,0.0);
   Fluxqe_z_cc.initialize(nXcc,nZcc,0.0);
   Fluxqi_x_cc.initialize(nXcc,nZcc,0.0);
   Fluxqi_z_cc.initialize(nXcc,nZcc,0.0);
   omegatau.initialize(nXcc,nZcc,0.0);
   omegatau2.initialize(nXcc,nZcc,0.0);
   omegatau4.initialize(nXcc,nZcc,0.0);
   omegatauFactor.initialize(nXcc,nZcc,0.0);
   Xgrid.DDZ(dTedz,Te);
   Xgrid.DDX(dTedx,Te);
   Xgrid.DDZ(dTidz,Ti);
   Xgrid.DDX(dTidx,Ti);
   Xgrid.DDZ(dPedz,Pe);
   Xgrid.DDX(dPedx,Pe);
   Xgrid.DDZ(dPidz,Pi);
   Xgrid.DDX(dPidx,Pi);
   Xgrid.communicate(dTedx);
   Xgrid.communicate(dTidx);
   Xgrid.communicate(dPedx);
   Xgrid.communicate(dPidx);


   //omegatau = 0.0213*1.84e3*By*pow(Te,1.5); // omega_ce*tau_ei
   //omegatau2 = omegatau*omegatau;
   //omegatau4 = omegatau2*omegatau2;
   //omegatauFactor  = omegatau4 + 2.0*21.67/5.0*omegatau2;
   //omegatauFactor /= omegatau4 + 14.79*omegatau2 + 3.7703;   
   
   omegatau = 1.0e2*By; // omega_ce*tau_ei
   omegatau2 = omegatau*omegatau;
   omegatau4 = omegatau2*omegatau2;
   //omegatauFactor  = omegatau2;
   //omegatauFactor /= omegatau2 + 3.7703;   
   
   omegatauFactor  = omegatau2*(omegatau2 + 4.65*(gamma0-1.0)/gamma0);
   omegatauFactor /= omegatau4 + 2.70*omegatau2 + 0.677;


   qe_x = -gamma0/(gamma0-1.0)*Pe*lambda0/By*dTedz*omegatauFactor;  
   qe_z =  gamma0/(gamma0-1.0)*Pe*lambda0/By*dTedx*omegatauFactor;  
   qi_x =  gamma0/(gamma0-1.0)*Pi*lambda0/By*dTidz*omegatauFactor;  
   qi_z = -gamma0/(gamma0-1.0)*Pi*lambda0/By*dTidx*omegatauFactor;  

   Fluxqe_x_cc = qe_x*rcc;
   Fluxqe_z_cc = qe_z;
   Fluxqi_x_cc = qi_x*rcc;
   Fluxqi_z_cc = qi_z;
   
   Xgrid.DDX(Sqe_x,qe_x*rcc);
   Xgrid.DDZ(Sqe_z,qe_z);
   Xgrid.DDX(Sqi_x,qi_x*rcc);
   Xgrid.DDZ(Sqi_z,qi_z);

   divqe = Sqe_x/rcc + Sqe_z;
   divqi = Sqi_x/rcc + Sqi_z;
   
   Xgrid.communicate(divqe);
   Xgrid.communicate(divqi);

   //Qi = 3.0*nuTherm0*N*(Te-Ti);
   Qi = 2.0/(gamma0-1.0)*nuTherm0*N*(Te-Ti);
   Qe = -Qi;

   //   compute ion gyro-viscosicity fluxes
   //
   matrix2D<double> dVxdx, dVzdx, dVxdz, dVzdz;
   matrix2D<double> Wrz, Wzz, Wrr, divV;
   dVxdx.initialize(nXcc,nZcc,0.0);
   dVzdx.initialize(nXcc,nZcc,0.0);
   dVxdz.initialize(nXcc,nZcc,0.0);
   dVzdz.initialize(nXcc,nZcc,0.0);
   Wrz.initialize(nXcc,nZcc,0.0);
   Wrr.initialize(nXcc,nZcc,0.0);
   Wzz.initialize(nXcc,nZcc,0.0);
   divV.initialize(nXcc,nZcc,0.0);

   if(modelGyroVisc) {
      Xgrid.DDX(dVxdx,Vx);
      Xgrid.DDX(dVzdx,Vz);
      Xgrid.DDZ(dVxdz,Vx);
      Xgrid.DDZ(dVzdz,Vz);
      setXminBoundary(dVxdx, 0.0, 0.0);   
      setXmaxBoundary(dVxdx, 0.0, 0.0);   
      setXminBoundary(dVzdx, 0.0, 0.0);   
      setXmaxBoundary(dVzdx, 0.0, 0.0);   
      Xgrid.communicate(dVxdx);
      Xgrid.communicate(dVzdx);
      Wrz = dVxdz + dVzdx;
      Wrr = 2.0*dVxdx - 2.0/3.0*divV;
      Wzz = 2.0*dVzdz - 2.0/3.0*divV;
      eta3 = 0.0*lambda0*Pi/By/2.0*omegatauFactor;   
   }

   Vhallz = 0.0 - lambda0/N*Jz;
   Vhallx = 0.0 - lambda0/N*Jx;
   Vez = Vz + Vhallz;
   Vex = Vx + Vhallx;
  
   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeedx = sqrt(Vx*Vx) + Cs;    
   Cspeedz = sqrt(Vz*Vz) + Cs;    
   Cspeede = sqrt(Vhallx*Vhallx + Vhallz*Vhallz);

   /*
   for (auto i=0; i<nXcc; i++) {
      for (auto j=0; j<nZcc; j++) {
         Cspeed(i,j) = max(Cspeed0i(i,j),Cspeed0e(i,j));
      }
   }
   */

   //   compute cell-center fluxes
   //
   FluxNcc_x = rcc*Mx;
   FluxNcc_z = Mz;
   FluxMxcc_x = rcc*(Mx*Vx + P - eta3*Wrz);
   FluxMxcc_z = Mx*Vz - eta3/2.0*(Wzz-Wrr);
   FluxMzcc_x = rcc*(Mz*Vx - eta3/2.0*(Wzz-Wrr));
   FluxMzcc_z = Mz*Vz + P + eta3*Wrz;
   FluxSicc_x  = rcc*Vx*Si;
   FluxSecc_z  = Vz*Se;
   FluxSecc_z0 = Vz*Se;
   FluxSicc_z = Vz*Si;
   FluxSecc_x0 = rcc*Vx*Se;
   FluxSecc_z0 = Vz*Se;
   FluxSecc_hall_x = rcc*Vhallx*Se;
   FluxSecc_hall_z = Vhallz*Se;
   
   FluxBycc_x = -Ez;
   FluxBycc_z = Ex;
  
   Ez0 = -Vx*By;
   Ex0 = Vz*By;
   
   Ezprime = Ez-Ez0;
   Exprime = Ex-Ex0;
   
   FluxBycc_x0 = -Ez0;
   FluxBycc_z0 = Ex0;
   FluxBycc_xp = -Ezprime;
   FluxBycc_zp = Exprime;


   // compute advective flux using
   // specified scheme from input file
   //
   //
   Nsub = 2; // hard code this for force balance at initiation
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN_x, FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxNcc_x, Cspeedx,rcc*N, 0,Nsub);
      Xgrid.computeFluxTVD(FluxMx_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMxcc_x,Cspeedx,rcc*Mx,0,Nsub);
      Xgrid.computeFluxTVD(FluxMz_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMzcc_x,Cspeedx,rcc*Mz,0,Nsub);
      Xgrid.computeFluxTVD(FluxSe_x, FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxSecc_x0,Cspeedx,rcc*Se,0,Nsub);
      Xgrid.computeFluxTVD(FluxSi_x, FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxSicc_x,Cspeedx,rcc*Si,0,Nsub);
      Xgrid.computeFluxTVD(FluxBy_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxBycc_x0,Cspeedx,By,   0,Nsub);
      Xgrid.computeFluxTVD(Fluxqe_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           Fluxqe_x_cc,Cspeedx,rcc*Se,  0,Nsub);
      Xgrid.computeFluxTVD(Fluxqi_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           Fluxqi_x_cc,Cspeedx,rcc*Si,  0,Nsub);
      //
      Xgrid.computeFluxTVD(FluxN_z, FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxNcc_z, Cspeedz,N, 1,Nsub);
      Xgrid.computeFluxTVD(FluxMx_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxMxcc_z,Cspeedz,Mx,1,Nsub);
      Xgrid.computeFluxTVD(FluxMz_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxMzcc_z,Cspeedz,Mz,1,Nsub);
      Xgrid.computeFluxTVD(FluxSe_z, FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxSecc_z0,Cspeedz,Se,1,Nsub);
      Xgrid.computeFluxTVD(FluxSi_z, FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxSicc_z,Cspeedz,Si,1,Nsub);
      Xgrid.computeFluxTVD(FluxBy_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxBycc_z0,Cspeedz,By,1,Nsub);
      Xgrid.computeFluxTVD(Fluxqe_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           Fluxqe_z_cc,Cspeedz,Se,  1,Nsub);
      Xgrid.computeFluxTVD(Fluxqi_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           Fluxqi_z_cc,Cspeedz,Si,  1,Nsub);
   }
   else {
      Xgrid.InterpToCellEdges(FluxN_x,  FluxNcc_x,  Vx, advScheme0, 0);
      Xgrid.InterpToCellEdges(FluxMx_x, FluxMxcc_x, Vx, advScheme0, 0);
      Xgrid.InterpToCellEdges(FluxMz_x, FluxMzcc_x, Vx, advScheme0, 0);
      Xgrid.InterpToCellEdges(FluxSe_x, FluxSecc_x, Vx, advScheme0, 0);
      Xgrid.InterpToCellEdges(FluxSi_x, FluxSicc_x, Vx, advScheme0, 0);
      Xgrid.InterpToCellEdges(FluxBy_x, FluxBycc_x, Vx, advScheme0, 0);
      //
      Xgrid.InterpToCellEdges(FluxN_z,  FluxNcc_z,  Vz, advScheme0, 1);
      Xgrid.InterpToCellEdges(FluxMx_z, FluxMxcc_z, Vz, advScheme0, 1);
      Xgrid.InterpToCellEdges(FluxMz_z, FluxMzcc_z, Vz, advScheme0, 1);
      Xgrid.InterpToCellEdges(FluxSe_z, FluxSecc_z, Vz, advScheme0, 1);
      Xgrid.InterpToCellEdges(FluxSi_z, FluxSicc_z, Vz, advScheme0, 1);
      Xgrid.InterpToCellEdges(FluxBy_z, FluxBycc_z, Vz, advScheme0, 1);
   } 
   Xgrid.InterpToCellEdges(FluxBy_xp, FluxBycc_xp, Vx, "C2", 0);
   Xgrid.InterpToCellEdges(FluxBy_zp, FluxBycc_zp, Vz, "C2", 1);
   FluxBy_x = FluxBy_x + FluxBy_xp;   
   FluxBy_z = FluxBy_z + FluxBy_zp;   
   

   // compute Hall velocity flux for electron heat eqn
   //
   //Xgrid.computeFluxTVD(FluxSe_hall_x, FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
   //                     FluxSecc_hall_x,Cspeed0e,rcc*Se,0,Nsub);
   //Xgrid.computeFluxTVD(FluxSe_hall_z, FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
   //                     FluxSecc_hall_z,Cspeed0e,Se,1,Nsub);
   Xgrid.InterpToCellEdges(FluxSe_hall_x, FluxSecc_hall_x, Vhallx, "C2", 0);
   Xgrid.InterpToCellEdges(FluxSe_hall_z, FluxSecc_hall_z, Vhallz, "C2", 1);
   FluxSe_x += FluxSe_hall_x;
   FluxSe_z += FluxSe_hall_z;


   if(procID==0)  {
      setXminBoundary(FluxN_x,  0.0, 0.0);   
      setXminBoundary(FluxMx_x, 0.0, 0.0);   
      setXminBoundary(FluxMz_x, 0.0, 0.0);   
      setXminBoundary(FluxSe_x, 0.0, 0.0);
      setXminBoundary(FluxSi_x, 0.0, 0.0);
      setXminBoundary(Fluxqe_x, 0.0, 0.0);
      setXminBoundary(Fluxqi_x, 0.0, 0.0);
      //setXminBoundary(FluxBy_x, 0.0, 0.0);
   }   
   if(procID==numProcs-1)  {
      setXmaxBoundary(FluxN_x,  0.0, 0.0);   
      //setXmaxBoundary(FluxMx_x, (P.at(2)+P.at(1))/2.0, 0.0);   
      setXmaxBoundary(FluxMz_x, 0.0, 0.0);   
      setXmaxBoundary(FluxSe_x, 0.0, 0.0);
      setXmaxBoundary(FluxSi_x, 0.0, 0.0);
      setXmaxBoundary(Fluxqe_x, 0.0, 0.0);
      setXmaxBoundary(Fluxqi_x, 0.0, 0.0);
      setXmaxBoundary(FluxBy_x, 0.0, 0.0);
   }   

   Xgrid.communicate(FluxL_x);   
   Xgrid.communicate(FluxR_x);   
   Xgrid.communicate(FluxRatio_x);   
   Xgrid.communicate(FluxLim_x);   

   // update stress tensor stuff
   //
   computeViscousStress(Xgrid);      
   
} // end computeFluxes

void computeViscousStress(const domainGrid& Xgrid)
{
   const int nXcc = Xgrid.nXcc;
   const int nXce = Xgrid.nXce;
   const int nXg = Xgrid.nXg;
   const int nZcc = Xgrid.nZcc;
   const int nZce = Xgrid.nZce;
   const int nZg = Xgrid.nZg;
   
   vector<double> Xcc, Xce;
   Xcc = Xgrid.Xcc;
   Xce = Xgrid.Xce;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   matrix2D<double> dVxdx, dVzdz, dVxdz, dVzdx, divV;
   matrix2D<double> dVxdx_cex, dVzdz_cex, dVxdz_cex, dVzdx_cex, divV_cex;
   matrix2D<double> eta0_cex, eta1_cex, eta3_cex;
   matrix2D<double> dVxdx_cez, dVzdz_cez, dVxdz_cez, dVzdx_cez, divV_cez;
   matrix2D<double> eta0_cez, eta1_cez, eta3_cez;
   matrix2D<double> rdivV_cex, Vx_cex, Vz_cez, Vx_cez;
   matrix2D<double> Wrr, Wzz, Wrz, Wthth, rWthth_cex, Wthth_cez;
   matrix2D<double> W0thth, W0rr, W0zz;
   matrix2D<double> W1rr, W1rz, W1zz;
   matrix2D<double> FluxP0x_x, FluxP0z_z, QP0, divP0x, divP0z;
   matrix2D<double> FluxP1x_x, FluxP1x_z, FluxP1z_x, FluxP1z_z, QP1, divP1x, divP1z;   
 

   dVxdx.initialize(nXcc,nZcc,0.0);
   dVzdz.initialize(nXcc,nZcc,0.0);
   dVxdz.initialize(nXcc,nZcc,0.0);
   dVzdx.initialize(nXcc,nZcc,0.0);
   divV.initialize(nXcc,nZcc,0.0);
   Wthth.initialize(nXcc,nZcc,0.0);
   Wrr.initialize(nXcc,nZcc,0.0);
   Wzz.initialize(nXcc,nZcc,0.0);
   Wrz.initialize(nXcc,nZcc,0.0);
   W0thth.initialize(nXcc,nZcc,0.0);
   W0rr.initialize(nXcc,nZcc,0.0);
   W0zz.initialize(nXcc,nZcc,0.0);
   QP0.initialize(nXcc,nZcc,0.0);
   divP0x.initialize(nXcc,nZcc,0.0);
   divP0z.initialize(nXcc,nZcc,0.0);
   W1rr.initialize(nXcc,nZcc,0.0);
   W1rz.initialize(nXcc,nZcc,0.0);
   W1zz.initialize(nXcc,nZcc,0.0);
   QP1.initialize(nXcc,nZcc,0.0);
   divP1x.initialize(nXcc,nZcc,0.0);
   divP1z.initialize(nXcc,nZcc,0.0);//
   //
   // 
   dVxdx_cex.initialize(nXce,nZcc,0.0);
   dVzdz_cex.initialize(nXce,nZcc,0.0);
   dVxdz_cex.initialize(nXce,nZcc,0.0);
   dVzdx_cex.initialize(nXce,nZcc,0.0);
   divV_cex.initialize(nXce,nZcc,0.0);
   eta0_cex.initialize(nXce,nZcc,0.0);
   eta1_cex.initialize(nXce,nZcc,0.0);
   eta3_cex.initialize(nXce,nZcc,0.0);
   rWthth_cex.initialize(nXce,nZcc,0.0);
   rdivV_cex.initialize(nXce,nZcc,0.0);
   Vx_cex.initialize(nXce,nZcc,0.0);
   //
   dVxdx_cez.initialize(nXcc,nZce,0.0);
   dVzdz_cez.initialize(nXcc,nZce,0.0);
   dVxdz_cez.initialize(nXcc,nZce,0.0);
   dVzdx_cez.initialize(nXcc,nZce,0.0);
   divV_cez.initialize(nXcc,nZce,0.0);
   eta0_cez.initialize(nXcc,nZce,0.0);
   eta1_cez.initialize(nXcc,nZce,0.0);
   eta3_cez.initialize(nXcc,nZce,0.0);
   Wthth_cez.initialize(nXcc,nZce,0.0);
   Vx_cez.initialize(nXcc,nZce,0.0);
   //
   FluxP0x_x.initialize(nXce,nZcc,0.0);
   FluxP0z_z.initialize(nXcc,nZce,0.0);
   FluxP1x_x.initialize(nXce,nZcc,0.0);
   FluxP1x_z.initialize(nXcc,nZce,0.0);
   FluxP1z_x.initialize(nXce,nZcc,0.0);
   FluxP1z_z.initialize(nXcc,nZce,0.0);

   
   //   set BCs for Vx and Vz
   //
   setXminBoundary(Vx, 0.0, -1.0);   
   setXmaxBoundary(Vx, 0.0, -1.0);   
   setXminBoundary(Vz, 0.0, 1.0);   
   setXmaxBoundary(Vz, 0.0, 1.0);   
   Xgrid.communicate(Vx);   
   Xgrid.communicate(Vz);   
   //double taui0 = 1.0e-3;

   //   compute tensor values at cell-center
   //
   eta0 = 0.96*Pi*taui0;
   eta1 = 0.96*Pi*taui0;
   Xgrid.DDX(dVxdx,Vx);
   Xgrid.DDZ(dVzdz,Vz);
   Xgrid.DDX(dVzdx,Vz);
   Xgrid.DDZ(dVxdz,Vx);
   Xgrid.DDX(divV, rcc*Vx);
   divV = divV/rcc + dVzdz;
   //
   Wrr = 2.0*dVxdx - 2.0/3.0*divV;
   Wthth = 2.0*Vx/rcc - 2.0/3.0*divV;
   Wzz = 2.0*dVzdz - 2.0/3.0*divV;
   Wrz = dVxdz + dVzdx;
   //
   W0rr   = -0.5*Wthth;
   W0thth = Wthth;
   W0zz   = -0.5*Wthth;
   W1rr   = 0.5*(Wrr - Wzz);
   W1zz   = -0.5*(Wrr - Wzz);
   W1rz   = Wrz;


   //   compute values at cell-edge in r
   //
   Xgrid.DDX(dVxdx_cex,Vx);
   Xgrid.DDX(dVzdx_cex,Vz);
   Xgrid.InterpToCellEdges(eta0_cex, eta0, Vx, "C2", 0);
   Xgrid.InterpToCellEdges(eta1_cex, eta1, Vx, "C2", 0);
   Xgrid.InterpToCellEdges(dVzdz_cex, dVzdz, Vx, "C2", 0);
   Xgrid.InterpToCellEdges(dVxdz_cex, dVxdz, Vx, "C2", 0);
   Xgrid.InterpToCellEdges(Vx_cex, Vx, Vx, "C2", 0);
   Xgrid.DDX(rdivV_cex,rcc*Vx);
   rdivV_cex = rdivV_cex + rce_x*dVzdz_cex;
   rWthth_cex = 2.0*Vx_cex - 2.0/3.0*rdivV_cex;
  
   //   compute values at cell-edge in z
   //
   Xgrid.DDZ(dVzdz_cez,Vz);
   Xgrid.DDZ(dVxdz_cez,Vx);
   Xgrid.InterpToCellEdges(eta0_cez, eta0, Vz, "C2", 1);
   Xgrid.InterpToCellEdges(eta1_cez, eta1, Vz, "C2", 1);
   Xgrid.InterpToCellEdges(dVxdx_cez, dVxdx, Vz, "C2", 1);
   Xgrid.InterpToCellEdges(dVzdx_cez, dVzdx, Vz, "C2", 1);
   Xgrid.InterpToCellEdges(Vx_cez, Vx, Vz, "C2", 1);
   Xgrid.DDX(divV_cez, rce_z*Vx_cez);
   divV_cez = divV_cez/rce_z + dVzdz_cez;
   Wthth_cez = 2.0*Vx_cez/rce_z - 2.0/3.0*divV_cez;


   //   compute divergence of stress tensor (divP)
   //   and heat source (Q = -P:nablaV)
   //
   FluxP0x_x = eta0_cex/2.0*rWthth_cex;
   FluxP0z_z = eta0_cez/2.0*Wthth_cez;
   FluxP1x_x = -eta1_cex*rce_x*(dVxdx_cex - dVzdz_cex);
   FluxP1x_z = -eta1_cez*(dVxdz_cez + dVzdx_cez);
   FluxP1z_x = -eta1_cex*rce_x*(dVxdz_cex + dVzdx_cex);
   FluxP1z_z = eta1_cez*(dVxdx_cez - dVzdz_cez);
   for (auto i=nXg; i<nXcc-nXg; i++) {
      for (auto j=nZg; j<nZcc-nZg; j++) {
         divP0x(i,j)  = (FluxP0x_x(i,j) - FluxP0x_x(i-1,j))/Xcc.at(i)/Xgrid.dX
                      + eta0(i,j)*Wthth(i,j)/Xcc.at(i);
         divP0z(i,j)  = (FluxP0z_z(i,j) - FluxP0z_z(i,j-1))/Xgrid.dZ;
         //
         divP1x(i,j)  = (FluxP1x_x(i,j) - FluxP1x_x(i-1,j))/Xcc.at(i)/Xgrid.dX
                      + (FluxP1x_z(i,j) - FluxP1x_z(i,j-1))/Xgrid.dZ;
         divP1z(i,j)  = (FluxP1z_x(i,j) - FluxP1z_x(i-1,j))/Xcc.at(i)/Xgrid.dX
                      + (FluxP1z_z(i,j) - FluxP1z_z(i,j-1))/Xgrid.dZ;
      }
   }
   QP0 = eta0/2.0*(W0rr*W0rr + W0thth*W0thth + W0zz*W0zz);
   QP1 = eta1/2.0*(W1rr*W1rr + 2.0*W1rz*W1rz + W1zz*W1zz);


   //  update total viscous terms for force and heat eqs
   //
   divPvis_x = divP0x + divP1x;
   divPvis_z = divP0z + divP1z;
   Qvis = QP0 + QP1;

} // end computeViscousStress


void computeFluxes_E(const domainGrid& Xgrid, const int order)
{

   const int nXcc = Xgrid.nXcc;
   //const int nXce = Xgrid.nXce;
   const int nZcc = Xgrid.nZcc;
   //const int nZce = Xgrid.nZce;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   matrix2D<double> FluxEzcc_x, FluxExcc_z, Cspeed0, Cspeed, Cspeed2, Cvac; 
 
   Cspeed0.initialize(nXcc,nZcc,0.0);
   Cspeed.initialize(nXcc,nZcc,0.0);
   Cspeed2.initialize(nXcc,nZcc,0.0);
   double Cvac0 = 1.0/sqrt(delta0);
   Cvac.initialize(nXcc,nZcc,Cvac0);
   FluxEzcc_x.initialize(nXcc,nZcc,0.0);
   FluxExcc_z.initialize(nXcc,nZcc,0.0);
   
   //  define derived variables
   //
   Vx = Mx/N;
   Vz = Mz/N;
   Pe  = Se*pow(N,gamma0-1.0);
   Pi  = Si*pow(N,gamma0-1.0);
   Te  = Pe/N;
   Ti  = Pi/N;
   P = Pe + Pi;
   Cs = sqrt(gamma0*P/N + By*By/N);

   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeed0  = sqrt(Vx*Vx+Vz*Vz); // adv flux jacobian for magnetic field
   Cspeed   = Cspeed0 + Cs;      // adv flux jacobian for all other variables
   

   FluxEzcc_x = -rcc*By;
   FluxExcc_z = By;
 

   //Xgrid.computeFluxTVD(FluxEz_x2,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
   //                     FluxEzcc_x,Cvac, rcc*Ez, 0,2);
   //Xgrid.computeFluxTVD(FluxEx_z2,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
   //                     FluxExcc_z,Cvac, Ex,     1,2);
      
   Xgrid.InterpToCellEdges(FluxEz_x, FluxEzcc_x, Vx, "C2", 0);
   Xgrid.InterpToCellEdges(FluxEx_z, FluxExcc_z, Vz, "C2", 1);
   
   if(procID==0)  {
     // setXnminBoundary(FluxEz_x, 0.0, 0.0);
     // setXminBoundary(FluxEz_x2, 0.0, 0.0);
   }   
  
   Xgrid.communicate(FluxEz_x);
   Xgrid.communicate(FluxEx_z);
   Xgrid.communicate(FluxEz_x2);
   Xgrid.communicate(FluxEx_z2);

   //FluxEx_z = FluxEx_z2;
   //FluxEz_x = FluxEz_x2;

}
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
   
   matrix2D<double> Cchar_x(Cs), Cchar_z(Cs);
   matrix2D<double> Cchar_Vex(Cs), Cchar_Vez(Cs);
   double Cmax_x, Cmax_z, Cmax_Vex, Cmax_Vez, Cmax, Clight;
   Cchar_x += abs(Vx);
   Cchar_z += abs(Vz);
   Cchar_Vex += abs(Vx-lambda0/N*Jx);
   Cchar_Vez += abs(Vz-lambda0/N*Jz);
   Cmax_x = max(Cchar_x);
   Cmax_z = max(Cchar_z);
   Cmax_Vex = max(Cchar_Vex);
   Cmax_Vez = max(Cchar_Vez);
   Cmax = max(Cmax_x,Cmax_z);
   Cmax = max(Cmax,Cmax_Vex);
   Cmax = max(Cmax,Cmax_Vez);
   Clight = 1.0/sqrt(delta0);
   //cout << "Cmax = " << Cmax << endl;

   const double dX = Xgrid.dX;
   const double dZ = Xgrid.dZ;
   
   double dtCFL_sound = 0.5*dX*dZ/(dX+dZ)/Cmax;
   double dtCFL_light = 0.5*dX*dZ/(dX+dZ)/Clight;
   double dtmax = min(dtCFL_sound,dtCFL_light); // should always be light
   //double dtmax = dtCFL_sound;
   
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);


   //dtSim = 1.0e-4;
   if(procID==0) {
      cout << "dtSim = " << dtSim << endl;
      if(dtCFL_sound==dtmax)
      cout << "warning: dtSim set by CFL_sound"<< endl;
   }
}

