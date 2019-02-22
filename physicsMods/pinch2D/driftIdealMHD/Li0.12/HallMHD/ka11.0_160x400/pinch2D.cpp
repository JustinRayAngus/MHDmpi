/***
 * 
 * Energy conservative algorithm
 *
 * physics module for 2D zpinch with m=0 using drift-MHD 
 * with electron inertia and displacment current (relaxation)
 * and with electron pressure and with seperate electron and
 * ion heat equations with diamagnetic heat fluxes and thermalization
 *
 * dN/dt  + d(r*Mr)/r/dr + d(Mz)/dz = 0
 * dMr/dt + d(r*(Mr*Ur + P))/r/dr + d(Mr*Uz)/dz = -Jz*By + P/r
 * dMz/dt + d(r*Mz*Ur)/r/dr  + d(Mz*Uz + P)/dz = Jr*By
 * dEi/dt + d(r*Ur*(Ei+Pi))/r/dr  + d(Uz*(Ei+Pi))/dz  + divqi  =  N*UdotE/lambda0 + Qi  
 * dEe/dt + d(r*Uer*(Ee+Pe))/r/dr + d(Uez*(Ee+Pe))/dz + divqe  = -N*UdotE/lambda0 + JdotE + Qe  
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
 * Pe = Ee*(gamma0-1)
 * Pi = (Ei - 0.5*(Mr^2 + Mz^2)/N)*(gamma0-1)
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
double lambda0, epsilon0, delta0, taui0, tauiVis0;
double nuTherm0, TempRatio0;
int Nsub, modelDriftTerms, modelGyroVisc;
matrix2D<double> N, Mx, Mz, Ee, Ei, By, Ez, Ex, Jz, Jx;  // time-evolving variables
matrix2D<double> P0, deltaP;                    // initial perturbation
matrix2D<double> eta, Cs, Vx, Vz, P, T, S;     // derived variables
matrix2D<double> Ez0, Jz0, Ex0, Jx0;        // derived variables
matrix2D<double> eta_x, eta_z, Jzcc, Ezcc, Jxcc, VxBy_x, VzBy_z;
matrix2D<double> Nold, Mxold, Mzold, Eeold, Eiold, Byold;
matrix2D<double> Exold, Jxold, Ezold, Jzold;
matrix2D<double> FluxR_x, FluxL_x, FluxRatio_x, FluxLim_x;
matrix2D<double> FluxR_z, FluxL_z, FluxRatio_z, FluxLim_z;  
matrix2D<double> FluxN_x, FluxEe_x, FluxEi_x, FluxMx_x, FluxMz_x, FluxBy_x, FluxEz_x;
matrix2D<double> FluxN_z, FluxEe_z, FluxEi_z, FluxMx_z, FluxMz_z, FluxBy_z, FluxEx_z;
matrix2D<double> Flux_qe_x, Flux_qi_x, Flux_qe_z, Flux_qi_z;
matrix2D<double> FluxEz_x2, FluxEx_z2, Jx02;
matrix2D<double> Fx, rcc, rce_z, rce_x;
//
matrix2D<double> Ve_drift_x, Ve_drift_z, Vi_drift_x, Vi_drift_z;
matrix2D<double> Pe, Pi, Te, Ti, Se, Si, divqe, divqi, Qe, Qi;
matrix2D<double> Vex, Vez, Vhallx, Vhallz;
matrix2D<double> NUdotE, JdotE, Eisource, Eesource;
matrix2D<double> qix, qiz, qix0, qiz0, qixold, qizold;
matrix2D<double> qex, qez, qex0, qez0, qexold, qezold;
matrix2D<double> kappae_cex, kappai_cex, kappae_cez, kappai_cez;
matrix2D<double> kappae_coperp, kappai_coperp;
matrix2D<double> qix_coperp, qiz_coperp, qex_coperp, qez_coperp;
//
matrix2D<double> eta0, eta1, eta3;
matrix2D<double> Qvis, divPvis_x, divPvis_z;


// Set lower threshold values for N, T, and S
double Nthresh=1.0e-4, Tthresh=2.5e-2;
double mM=9.1094e-31/1.6726e-27;
double tauiEff_min = 1.0e-4, taueEff_min = 1.0e-4*sqrt(mM/2.0);

void computeFluxes(const domainGrid&, const int);
void computeFluxes_E(const domainGrid&, const int);
void updateCollisionalHeatFluxes(const domainGrid&);
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
   //
   const int nXg = Xgrid.nXg;
   const int nZg = Xgrid.nZg;


   N.initialize(nXcc,nZcc,0.0);
   Nold.initialize(nXcc,nZcc,0.0);
   deltaP.initialize(nXcc,nZcc,0.0);
   P0.initialize(nXcc,nZcc,0.0);
   Mx.initialize(nXcc,nZcc,0.0);
   Mxold.initialize(nXcc,nZcc,0.0);
   Mz.initialize(nXcc,nZcc,0.0);
   Mzold.initialize(nXcc,nZcc,0.0);
   Se.initialize(nXcc,nZcc,0.0);
   Ee.initialize(nXcc,nZcc,0.0);
   Eeold.initialize(nXcc,nZcc,0.0);
   Si.initialize(nXcc,nZcc,0.0);
   Ei.initialize(nXcc,nZcc,0.0);
   Eiold.initialize(nXcc,nZcc,0.0);
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
   //   
   qix.initialize(nXce,nZcc,0.0);
   qiz.initialize(nXcc,nZce,0.0);
   qix0.initialize(nXce,nZcc,0.0);
   qiz0.initialize(nXcc,nZce,0.0);
   qixold.initialize(nXce,nZcc,0.0);
   qizold.initialize(nXcc,nZce,0.0);
   //
   qex.initialize(nXce,nZcc,0.0);
   qez.initialize(nXcc,nZce,0.0);
   qex0.initialize(nXce,nZcc,0.0);
   qez0.initialize(nXcc,nZce,0.0);
   qexold.initialize(nXce,nZcc,0.0);
   qezold.initialize(nXcc,nZce,0.0);
   //
   kappae_cex.initialize(nXce,nZcc,0.0);
   kappai_cex.initialize(nXce,nZcc,0.0);
   kappae_cez.initialize(nXcc,nZce,0.0);
   kappai_cez.initialize(nXcc,nZce,0.0);
   kappae_coperp.initialize(nXcc,nZcc,0.0);
   kappai_coperp.initialize(nXcc,nZcc,0.0);
   //
   qix_coperp.initialize(nXce,nZcc,0.0);
   qiz_coperp.initialize(nXcc,nZce,0.0);
   qex_coperp.initialize(nXce,nZcc,0.0);
   qez_coperp.initialize(nXcc,nZce,0.0);
   
   
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
   JdotE.initialize(nXcc,nZcc,0.0);
   NUdotE.initialize(nXcc,nZcc,0.0);
   Eisource.initialize(nXcc,nZcc,0.0);
   Eesource.initialize(nXcc,nZcc,0.0);

   //
   FluxN_x.initialize(nXce, nZcc, 0.0);
   FluxMx_x.initialize(nXce,nZcc, 0.0);
   FluxMz_x.initialize(nXce,nZcc, 0.0);
   FluxEe_x.initialize(nXce, nZcc, 0.0);
   FluxEi_x.initialize(nXce, nZcc, 0.0);
   FluxBy_x.initialize(nXce,nZcc, 0.0);
   FluxEz_x.initialize(nXce,nZcc, 0.0);
   FluxEz_x2.initialize(nXce,nZcc, 0.0);
   
   FluxN_z.initialize(nXcc, nZce, 0.0);
   FluxMx_z.initialize(nXcc,nZce, 0.0);
   FluxMz_z.initialize(nXcc,nZce, 0.0);
   FluxEe_z.initialize(nXcc, nZce, 0.0);
   FluxEi_z.initialize(nXcc, nZce, 0.0);
   FluxBy_z.initialize(nXcc,nZce, 0.0);
   FluxEx_z.initialize(nXcc,nZce, 0.0);
   FluxEx_z2.initialize(nXcc,nZce, 0.0);

   Flux_qe_x.initialize(nXce,nZcc,0.0);
   Flux_qi_x.initialize(nXce,nZcc,0.0);
   Flux_qe_z.initialize(nXcc,nZce,0.0);
   Flux_qi_z.initialize(nXcc,nZce,0.0);
 
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
      Json::Value tauiVisVal = Phys.get("tauiVis",defValue);
      Json::Value modelGyroViscVal   = Phys.get("modelGyroVisc",defValue);
      Json::Value modelDriftTermsVal = Phys.get("modelDriftTerms",defValue);
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
      if(procID==0) {
         cout << "advection diff/interp scheme is " << advScheme0 << endl;
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
      
      tauiVis0 = tauiVisVal.asDouble();
      if(procID==0) cout << "tauiVis0/t0 = " << tauiVis0 << endl;
      
      Nsub = NsubVal.asInt();
      if(procID==0) cout << "Nsub = " << Nsub << endl;
      if(gamma0 < 1.0) {
         printf("ERROR: Nsub must be int >= 1\n");
         exit (EXIT_FAILURE);
      }
      
      modelDriftTerms = modelDriftTermsVal.asInt();
      if(procID==0) {
         if(modelDriftTerms==0) { 
            cout << "Not modeling drift terms => extended MHD"<< endl;
         }
         else if(modelDriftTerms==1) { 
            cout << "modeling drift terms"<< endl;
         }
         else{   
            printf("ERROR: modelDriftTerms should be 0 or 1 \n");
         }
      }

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
   S = Se + Si;
   Ee = Pe/(gamma0-1.0);
   Ei = 0.5*(Mx*Mx + Mz*Mz)/N + Pi/(gamma0-1.0);
   Eeold  = Ee;
   Eiold  = Ei;
 
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
   Ex = Ex0 - lambda0/N*(Jz*By + dPedx);
   //Ez = Ez0 + lambda0/N*(Jx*By - dPedz);
   //setXminExtrap(Ez,0);
   //setXmaxBoundary(Ez, 0.0, -1.0);   
   setZboundaryPeriodic(Ex);  
   setZboundaryPeriodic(Ez);  
   Xgrid.communicate(Ez);
   Xgrid.communicate(Ex);
  
   // set initial electron velocity and ion momentum
   //
   Vhallz = 0.0 - lambda0/N*Jz;
   Vhallx = 0.0 - lambda0/N*Jx;
   Vez = Vz + Vhallz;
   Vex = Vx + Vhallx;
   Mz = N*Vz;
   Mzold = Mz;

   JdotE = Jx*Ex + Jz*Ez;
   Qi = 2.0/(gamma0-1.0)*nuTherm0*N*(Te-Ti);
   if(lambda0==0) {
      Eisource = N*Vz*(Jx*By - dPedz) - N*Vx*(Jz*By + dPedx) + Qi;
   } else {
      NUdotE = N*(Vx*Ex + Vz*Ez);
      Eisource = NUdotE/lambda0 + Qi;
   }
   Eesource = JdotE - Eisource;


   // compute fluxes before first time step
   //
   computeFluxes(Xgrid, Nsub);  
   qix = qix0;
   qiz = qiz0;
   qixold = qix;
   qizold = qiz;
   //
   qex = qex0;
   qez = qez0;
   qexold = qex;
   qezold = qez;


   // add stuff to output files
   //
   dataFile.add(rcc, "rcc", 0);  // force density x-direction 
   dataFile.add(Fx, "Fx", 1);  // force density x-direction 
   dataFile.add(N,  "N",  1);  // density 
   dataFile.add(deltaP,  "deltaP",  0);  // density perturbation 
   dataFile.add(P0, "P0",  0);  // initial unperturbed pressure profile 
   dataFile.add(Mx, "Mx", 1);  // momentum density 
   dataFile.add(Mz, "Mz", 1);  // momentum density 
   dataFile.add(Ee, "Ee", 1);  // electron energy density
   dataFile.add(Ei, "Ei", 1);  // ion energy density
   dataFile.add(By, "By", 1);  // magnetic field
   dataFile.add(P,  "P",  1);  // pressure
   dataFile.add(Pe, "Pe", 1);  // pressure
   dataFile.add(Pi, "Pi", 1);  // pressure
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
   //dataFile.add(FluxSe_x, "FluxSe_x", 1);  
   dataFile.add(FluxBy_x, "FluxBy_x", 1);  
   dataFile.add(FluxEz_x, "FluxEz_x", 1);  
   dataFile.add(FluxEx_z, "FluxEx_z", 1);  
   //dataFile.add(FluxEz_x2, "FluxEz_x2", 1);  
   //dataFile.add(FluxEx_z2, "FluxEx_z2", 1);  
   //
   dataFile.add(FluxRatio_x, "FluxRatio_x", 1);  
   dataFile.add(FluxLim_x, "FluxLim_x", 1);  
   dataFile.add(FluxR_x, "FluxR_x", 1);
   dataFile.add(FluxL_x, "FluxL_x", 1);
   //
   //dataFile.add(divqe,"divqe",1);
   //dataFile.add(divqi,"divqi",1);
   dataFile.add(qix,"qix",1);
   dataFile.add(qiz,"qiz",1);
   dataFile.add(qix0,"qix0",1);
   dataFile.add(qiz0,"qiz0",1);
   dataFile.add(qex,"qex",1);
   dataFile.add(qez,"qez",1);
   dataFile.add(qex0,"qex0",1);
   dataFile.add(qez0,"qez0",1);
   dataFile.add(qix_coperp,"qix_coperp",1); // multiplied by r !!!
   dataFile.add(qiz_coperp,"qiz_coperp",1);
   dataFile.add(qex_coperp,"qex_coperp",1); // multiplied by r !!!
   dataFile.add(qez_coperp,"qez_coperp",1);
   dataFile.add(kappae_cex,"kappae_cex",1);
   dataFile.add(kappai_cex,"kappai_cex",1);
   dataFile.add(kappae_coperp,"kappae_coperp",1);
   dataFile.add(kappai_coperp,"kappai_coperp",1);
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
   double thisC_x, thisC_z;

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
         Ee(i,j) = Eeold(i,j) - thisdt*(FluxEe_x(i,j) - FluxEe_x(i-1,j))/Xcc.at(i)/Xgrid.dX
                              - thisdt*(FluxEe_z(i,j) - FluxEe_z(i,j-1))/Xgrid.dZ
                              + thisdt*(Eesource(i,j) - 0.0*divqe(i,j));
         Ei(i,j) = Eiold(i,j) - thisdt*(FluxEi_x(i,j) - FluxEi_x(i-1,j))/Xcc.at(i)/Xgrid.dX
                              - thisdt*(FluxEi_z(i,j) - FluxEi_z(i,j-1))/Xgrid.dZ
                              + thisdt*(Eisource(i,j) - 0.0*divqi(i,j) + Qvis(i,j));
	 By(i,j) = Byold(i,j) - thisdt*(FluxBy_x(i,j) - FluxBy_x(i-1,j))/Xgrid.dX
	                      - thisdt*(FluxBy_z(i,j) - FluxBy_z(i,j-1))/Xgrid.dZ;
	 
	 //if(N(i,j)<=Nthresh) N(i,j) = Nthresh;
	 //if(M(i,j)<0.0) M(i,j) = 0.0;

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
         setXminExtrap(Ee,0);
         setXminExtrap(Ei,0);
         setXminBoundary(By,0.0,-1.0);   
      }
     
      if(procID==numProcs-1) {
         setXmaxBoundary(Mx, 0.0, -1.0);   
         setXmaxExtrap(N, 0);   
         setXmaxExtrap(Mz,0);   
         setXmaxExtrap(Ee,0);   
         setXmaxExtrap(Ei,0);   
         setXmaxBy(By,-1);
      }
      setZboundaryPeriodic(N);
      setZboundaryPeriodic(Mx);
      setZboundaryPeriodic(Mz);
      setZboundaryPeriodic(Ee);
      setZboundaryPeriodic(Ei);
      setZboundaryPeriodic(By);
      //
      Xgrid.communicate(Fx);
      Xgrid.communicate(N);
      Xgrid.communicate(Mx);
      Xgrid.communicate(Mz);
      Xgrid.communicate(Ee);
      Xgrid.communicate(Ei);
      Xgrid.communicate(By);


      // Now update E and J, which are done implicitly
      //
      computeFluxes_E(Xgrid, Nsub); // compute electric field fluxes 
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
      for (auto i=nXg-1; i<nXcc-nXg; i++) {
         for (auto j=nZg-1; j<nZcc-nZg; j++) {
	 
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

         //   update heat fluxes using first order relaxation scheme
         //
         thisC_x = thisdt/delta0/kappai_cex(i,j);
         thisC_z = thisdt/delta0/kappai_cez(i,j);
         qix(i,j)  = qixold(i,j) + thisC_x*qix0(i,j);
         qix(i,j) /= (1.0 + thisC_x);
         qiz(i,j)  = qizold(i,j) + thisC_z*qiz0(i,j);
         qiz(i,j) /= (1.0 + thisC_z);
         //
         thisC_x = thisdt/delta0/kappae_cex(i,j);
         thisC_z = thisdt/delta0/kappae_cez(i,j);
         qex(i,j)  = qexold(i,j) + thisC_x*qex0(i,j);
         qex(i,j) /= (1.0 + thisC_x);
         qez(i,j)  = qezold(i,j) + thisC_z*qez0(i,j);
         qez(i,j) /= (1.0 + thisC_z);
         

         //   update heat fluxes using 2nd order du-fort frenkel scheme
         //
         /*
         if(n==Nsub) {
            thisC_x = thisdt/delta0/kappai_cex(i,j);
            thisC_z = thisdt/delta0/kappai_cez(i,j);
            qix(i,j)  = qixold(i,j)*(1.0-thisC_x/2.0) + thisC_x*qix0(i,j);
            qix(i,j) /= (1.0 + thisC_x/2.0);
            qiz(i,j)  = qizold(i,j)*(1.0-thisC_z/2.0) + thisC_z*qiz0(i,j);
            qiz(i,j) /= (1.0 + thisC_z/2.0);
            //
            thisC_x = thisdt/delta0/kappae_cex(i,j);
            thisC_z = thisdt/delta0/kappae_cez(i,j);
            qex(i,j)  = qexold(i,j)*(1.0-thisC_x/2.0) + thisC_x*qex0(i,j);
            qex(i,j) /= (1.0 + thisC_x/2.0);
            qez(i,j)  = qezold(i,j)*(1.0-thisC_z/2.0) + thisC_z*qez0(i,j);
            qez(i,j) /= (1.0 + thisC_z/2.0);
         }
         */

         }
      }
      if(procID==0) {
         setXminExtrap(Ex,0);
         setXminExtrap(Ez,0);
         setXminExtrap(Jx,0);
         setXminExtrap(Jz,0);
         setXminBoundary(qix,0.0,0.0);   
         setXminBoundary(qex,0.0,0.0);   
         setXminExtrap(qiz,0);
         setXminExtrap(qez,0);
      }
      if(procID==numProcs-1) {
         setXmaxExtrap(Ex, 0);   
         setXmaxBoundary(Ez, 0.0, -1.0);   
         setXmaxExtrap(Jx,0);
         setXmaxExtrap(Jz,0);
         setXmaxBoundary(qix,0.0,0.0);
         setXmaxExtrap(qiz,0);
         setXmaxBoundary(qex,0.0,0.0);
         setXmaxExtrap(qez,0);
      }
      setZboundaryPeriodic(Ex);
      setZboundaryPeriodic(Ez);
      setZboundaryPeriodic(Jx);
      setZboundaryPeriodic(Jz);
      setZboundaryPeriodic(qix);
      setZboundaryPeriodic(qiz);
      setZboundaryPeriodic(qex);
      setZboundaryPeriodic(qez);
      //
      Xgrid.communicate(Ex);
      Xgrid.communicate(Ez);
      Xgrid.communicate(Jx);
      Xgrid.communicate(Jz);
      Xgrid.communicate(qix);
      Xgrid.communicate(qiz);
      Xgrid.communicate(qex);
      Xgrid.communicate(qez);
      
      Xgrid.communicate(Jx0);
      Xgrid.communicate(Jz0);
      Xgrid.communicate(qix0);
      Xgrid.communicate(qiz0);
      Xgrid.communicate(qex0);
      Xgrid.communicate(qez0);


      // compute fluxes using partialy updated fields at n+1/2 for 2nd order
      // or n+1 for 1st order
      computeFluxes(Xgrid, Nsub); // compute second order fluxes at n+1/2

   } // finish subcycle steps for N, M, S, and B

   // update old fields
   Nold = N;
   Mxold = Mx;
   Mzold = Mz;
   Eeold = Ee;
   Eiold = Ei;
   Byold = By;
   Exold = Ex;
   Ezold = Ez;
   Jxold = Jx;
   Jzold = Jz;
   qixold = qix;
   qizold = qiz;
   qexold = qex;
   qezold = qez;
   
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
   
   matrix2D<double> FluxNcc_x, FluxMxcc_x, FluxMzcc_x, FluxEecc_x, FluxEicc_x, FluxBycc_x;
   matrix2D<double> FluxNcc_z, FluxMxcc_z, FluxMzcc_z, FluxEecc_z, FluxEicc_z, FluxBycc_z;
   matrix2D<double> FluxEecc_hall_x, FluxEecc_hall_z, FluxEecc_x0, FluxEecc_z0;
   matrix2D<double> CspeedBx, CspeedBz, Cspeedx, Cspeedz, Cspeedex, Cspeedez, Cvac; 
   matrix2D<double> FluxMx_gviscc_x, FluxMx_gviscc_z, FluxMz_gviscc_x, FluxMz_gviscc_z;
   matrix2D<double> FluxMx_gvis_x, FluxMx_gvis_z, FluxMz_gvis_x, FluxMz_gvis_z;
   matrix2D<double> nuTherm;
 
   CspeedBx.initialize(nXcc,nZcc,0.0);
   CspeedBz.initialize(nXcc,nZcc,0.0);
   Cspeedx.initialize(nXcc,nZcc,0.0);
   Cspeedz.initialize(nXcc,nZcc,0.0);
   Cspeedex.initialize(nXcc,nZcc,0.0);
   Cspeedez.initialize(nXcc,nZcc,0.0);
   double Cvac0 = 1.0/sqrt(delta0);
   Cvac.initialize(nXcc,nZcc,Cvac0);
   FluxNcc_x.initialize(nXcc,nZcc,0.0);
   FluxNcc_z.initialize(nXcc,nZcc,0.0);
   FluxMxcc_x.initialize(nXcc,nZcc,0.0);
   FluxMxcc_z.initialize(nXcc,nZcc,0.0);
   FluxMzcc_x.initialize(nXcc,nZcc,0.0);
   FluxMzcc_z.initialize(nXcc,nZcc,0.0);
   FluxEecc_x.initialize(nXcc,nZcc,0.0);
   FluxEecc_x0.initialize(nXcc,nZcc,0.0);
   FluxEecc_hall_x.initialize(nXcc,nZcc,0.0);
   FluxEicc_x.initialize(nXcc,nZcc,0.0);
   FluxEecc_z.initialize(nXcc,nZcc,0.0);
   FluxEecc_z0.initialize(nXcc,nZcc,0.0);
   FluxEecc_hall_z.initialize(nXcc,nZcc,0.0);
   FluxEicc_z.initialize(nXcc,nZcc,0.0);
   FluxBycc_x.initialize(nXcc,nZcc,0.0);
   FluxBycc_z.initialize(nXcc,nZcc,0.0);
   
   FluxMx_gviscc_x.initialize(nXcc,nZcc,0.0);
   FluxMx_gviscc_z.initialize(nXcc,nZcc,0.0);
   FluxMz_gviscc_x.initialize(nXcc,nZcc,0.0);
   FluxMz_gviscc_z.initialize(nXcc,nZcc,0.0);
   FluxMx_gvis_x.initialize(nXce,nZcc,0.0);
   FluxMx_gvis_z.initialize(nXcc,nZce,0.0);
   FluxMz_gvis_x.initialize(nXce,nZcc,0.0);
   FluxMz_gvis_z.initialize(nXcc,nZce,0.0);
   nuTherm.initialize(nXcc,nZcc,0.0);


   matrix2D<double> Ez0, Ex0, Ezprime, Exprime, FluxBycc_x0, FluxBycc_z0;
   matrix2D<double> FluxBycc_xp, FluxBycc_zp, FluxBy_xp, FluxBy_zp;
   matrix2D<double> FluxEe_hall_x, FluxEe_hall_z;

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
   FluxEe_hall_x.initialize(nXce,nZcc,0.0);
   FluxEe_hall_z.initialize(nXcc,nZce,0.0);
   

   //  define derived variables
   //
   Vx = Mx/N;
   Vz = Mz/N;
   Pe = Ee*(gamma0-1.0);
   Pi = (Ei - 0.5*(Mx*Mx+Mz*Mz)/N)*(gamma0-1.0);
   Se  = Pe/pow(N,gamma0-1.0);
   Si  = Pi/pow(N,gamma0-1.0);
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
   //matrix2D<double> dTedx, dTedz, dTidx, dTidz;
   matrix2D<double> dPedx, dPedz, dPidx, dPidz;
   //dTedx.initialize(nXcc,nZcc,0.0);
   //dTedz.initialize(nXcc,nZcc,0.0);
   //dTidx.initialize(nXcc,nZcc,0.0);
   //dTidz.initialize(nXcc,nZcc,0.0);
   dPedx.initialize(nXcc,nZcc,0.0);
   dPedz.initialize(nXcc,nZcc,0.0);
   dPidx.initialize(nXcc,nZcc,0.0);
   dPidz.initialize(nXcc,nZcc,0.0);
   //Xgrid.DDZ(dTedz,Te);
   //Xgrid.DDX(dTedx,Te);
   //Xgrid.DDZ(dTidz,Ti);
   //Xgrid.DDX(dTidx,Ti);
   Xgrid.DDZ(dPedz,Pe);
   Xgrid.DDX(dPedx,Pe);
   Xgrid.DDZ(dPidz,Pi);
   Xgrid.DDX(dPidx,Pi);
   //Xgrid.communicate(dTedx);
   //Xgrid.communicate(dTidx);
   //Xgrid.communicate(dPedx);
   //Xgrid.communicate(dPidx);

   
   Vhallz = 0.0 - lambda0/N*Jz;
   Vhallx = 0.0 - lambda0/N*Jx;
   Vez = Vz + Vhallz;
   Vex = Vx + Vhallx;
   
   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   //CspeedBx = sqrt(Vx*Vx);    
   //CspeedBz = sqrt(Vz*Vz);    
   Cspeedx = sqrt(Vx*Vx + 0.0*Vz*Vz) + Cs;    
   Cspeedz = sqrt(Vz*Vz + 0.0*Vx*Vx) + Cs;    
   Cspeedex = sqrt(Vex*Vex) + Cs;
   Cspeedez = sqrt(Vez*Vez) + Cs;
   CspeedBx = Cspeedx;  
   CspeedBz = Cspeedz;  
 

   //   compute cell-center fluxes
   //
   FluxNcc_x = rcc*Mx;
   FluxNcc_z = Mz;
   FluxMxcc_x = rcc*(Mx*Vx + P);
   FluxMxcc_z = Mx*Vz;
   FluxMzcc_x = rcc*Mz*Vx;
   FluxMzcc_z = Mz*Vz + P;
   FluxEicc_x = rcc*Vx*(Ei+Pi);
   FluxEicc_z = Vz*(Ei+Pi);
   //
   FluxEecc_x0 = rcc*Vx*(Ee+Pe);
   FluxEecc_z0 = Vz*(Ee+Pe);
   FluxEecc_hall_x = rcc*Vhallx*(Ee+Pe);
   FluxEecc_hall_z = Vhallz*(Ee+Pe);
   //FluxEecc_x0 += FluxEecc_hall_x;   
   //FluxEecc_z0 += FluxEecc_hall_z;   

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
                           FluxNcc_x, Cspeedx,rcc*N, "vanleer",0,Nsub);
      Xgrid.computeFluxTVD(FluxMx_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMxcc_x,Cspeedx,rcc*Mx,"vanleer",0,Nsub);
      Xgrid.computeFluxTVD(FluxMz_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMzcc_x,Cspeedx,rcc*Mz,"vanleer",0,Nsub);
      Xgrid.computeFluxTVD(FluxEe_x, FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxEecc_x0,Cspeedx,rcc*Ee,"vanleer",0,Nsub);
      Xgrid.computeFluxTVD(FluxEi_x, FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxEicc_x,Cspeedx,rcc*Ei,"vanleer",0,Nsub);
      Xgrid.computeFluxTVD(FluxBy_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxBycc_x0,Cspeedx,By,   "vanleer",0,Nsub);
      //
      Xgrid.computeFluxTVD(FluxN_z, FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxNcc_z, Cspeedz,N, "vanleer",1,Nsub);
      Xgrid.computeFluxTVD(FluxMx_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxMxcc_z,Cspeedz,Mx,"vanleer",1,Nsub);
      Xgrid.computeFluxTVD(FluxMz_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxMzcc_z,Cspeedz,Mz,"vanleer",1,Nsub);
      Xgrid.computeFluxTVD(FluxEe_z, FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxEecc_z0,Cspeedz,Ee,"vanleer",1,Nsub);
      Xgrid.computeFluxTVD(FluxEi_z, FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxEicc_z,Cspeedz,Ei,"vanleer",1,Nsub);
      Xgrid.computeFluxTVD(FluxBy_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxBycc_z0,Cspeedz,By,"vanleer",1,Nsub);
   }
   else {
      Xgrid.computeFluxTVDsimple(FluxN_x, FluxL_x,FluxR_x,FluxNcc_x, 
                                 Cspeedx,rcc*N, advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxMx_x,FluxL_x,FluxR_x,FluxMxcc_x,
                                 Cspeedx,rcc*Mx,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxMz_x,FluxL_x,FluxR_x,FluxMzcc_x,
                                 Cspeedx,rcc*Mz,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxEi_x,FluxL_x,FluxR_x,FluxEicc_x,
                                 Cspeedx,rcc*Ei,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxEe_x,FluxL_x,FluxR_x,FluxEecc_x0,
                                 Cspeedx,rcc*Ee,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxBy_x,FluxL_x,FluxR_x,FluxBycc_x0,
                                 CspeedBx,By,   advScheme0,0);
      //
      Xgrid.computeFluxTVDsimple(FluxN_z, FluxL_z,FluxR_z,FluxNcc_z, 
                                 Cspeedz,N,  advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxMx_z,FluxL_z,FluxR_z,FluxMxcc_z,
                                 Cspeedz,Mx, advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxMz_z,FluxL_z,FluxR_z,FluxMzcc_z,
                                 Cspeedz,Mz, advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxEi_z,FluxL_z,FluxR_z,FluxEicc_z,
                                 Cspeedz,Ei, advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxEe_z,FluxL_z,FluxR_z,FluxEecc_z0,
                                 Cspeedz,Ee, advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxBy_z,FluxL_z,FluxR_z,FluxBycc_z0,
                                 CspeedBz,By,advScheme0,1);
   } 
   Xgrid.InterpToCellEdges(FluxBy_xp, FluxBycc_xp, Vx, "C2", 0);
   Xgrid.InterpToCellEdges(FluxBy_zp, FluxBycc_zp, Vz, "C2", 1);
   FluxBy_x = FluxBy_x + FluxBy_xp;   
   FluxBy_z = FluxBy_z + FluxBy_zp;   
   

   if(modelDriftTerms) {
      // compute Hall velocity flux for electron heat eqn
      //
      //Xgrid.computeFluxTVD(FluxEe_hall_x, FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
      //                     FluxEecc_hall_x,Cspeed0e,rcc*Ee,0,Nsub);
      //Xgrid.computeFluxTVD(FluxEe_hall_z, FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
      //                     FluxEecc_hall_z,Cspeed0e,Ee,1,Nsub);
      Xgrid.InterpToCellEdges(FluxEe_hall_x, FluxEecc_hall_x, Vhallx, "C2", 0);
      Xgrid.InterpToCellEdges(FluxEe_hall_z, FluxEecc_hall_z, Vhallz, "C2", 1);
      FluxEe_x += FluxEe_hall_x;
      FluxEe_z += FluxEe_hall_z;

      /*
      if(modelGyroVisc) {
         FluxMx_gviscc_x = 0.0 - rcc*eta3*Wrz;
         FluxMx_gviscc_z = 0.0 - eta3/2.0*(Wzz-Wrr);
         FluxMz_gviscc_x = 0.0 - rcc*eta3/2.0*(Wzz-Wrr);
         FluxMz_gviscc_z = eta3*Wrz;
         Xgrid.InterpToCellEdges(FluxMx_gvis_x, FluxMx_gviscc_x, Vx, "C2", 0);
         Xgrid.InterpToCellEdges(FluxMx_gvis_z, FluxMx_gviscc_z, Vz, "C2", 1);
         Xgrid.InterpToCellEdges(FluxMz_gvis_x, FluxMz_gviscc_x, Vx, "C2", 0);
         Xgrid.InterpToCellEdges(FluxMz_gvis_z, FluxMz_gviscc_z, Vz, "C2", 1);
         FluxMx_x += FluxMx_gvis_x;
         FluxMx_z += FluxMx_gvis_z;
         FluxMz_x += FluxMz_gvis_x;
         FluxMz_z += FluxMz_gvis_z;
      }
      */

   }

   // update energy fluxes to include collisional heat flux
   //
   //updateCollisionalHeatFluxes(Xgrid);
   //FluxEe_x += Flux_qe_x;
   //FluxEi_x += Flux_qi_x;
   //FluxEe_z += Flux_qe_z;
   //FluxEi_z += Flux_qi_z;
 
   if(procID==0)  {
      setXminBoundary(FluxN_x,  0.0, 0.0);   
      setXminBoundary(FluxMx_x, 0.0, 0.0);   
      setXminBoundary(FluxMz_x, 0.0, 0.0);   
      setXminBoundary(FluxEe_x, 0.0, 0.0);
      setXminBoundary(FluxEi_x, 0.0, 0.0);
      //setXminBoundary(FluxBy_x, 0.0, 0.0);
   }   
   if(procID==numProcs-1)  {
      setXmaxBoundary(FluxN_x,  0.0, 0.0);   
      //setXmaxBoundary(FluxMx_x, (P.at(2)+P.at(1))/2.0, 0.0);   
      setXmaxBoundary(FluxMz_x, 0.0, 0.0);   
      setXmaxBoundary(FluxEe_x, 0.0, 0.0);
      setXmaxBoundary(FluxEi_x, 0.0, 0.0);
      setXmaxBoundary(FluxBy_x, 0.0, 0.0);
   }   
   
   /*
   updateCollisionalHeatFluxes(Xgrid);
   FluxEe_x += Flux_qe_x;
   FluxEi_x += Flux_qi_x;
   FluxEe_z += Flux_qe_z;
   FluxEi_z += Flux_qi_z;
   */

   //Xgrid.communicate(FluxL_x);   
   //Xgrid.communicate(FluxR_x);   
   //Xgrid.communicate(FluxRatio_x);   
   //Xgrid.communicate(FluxLim_x);   
   

   // compute energy source terms
   //
   JdotE = Jx*Ex + Jz*Ez;
   nuTherm = mM*N/(taui0*pow(2.0*Te,1.5)*sqrt(mM/2.0)) + nuTherm0;
   Qi = 2.0/(gamma0-1.0)*nuTherm*N*(Te-Ti);
   Qe = -Qi;
   if(lambda0==0) {
      Eisource = N*Vz*(Jx*By - dPedz) - N*Vx*(Jz*By + dPedx) + Qi;
   } else {
      NUdotE = N*(Vx*Exprime + Vz*Ezprime);
      Eisource = NUdotE/lambda0 + Qi;
   }
   Eesource = JdotE - Eisource;
   //Eesource += lambda0*gamma0/(gamma0-1.0)*(Jx*dTedx + Jz*dTedz);

   
   // update stress tensor stuff
   //
   computeViscousStress(Xgrid);      
  
 
} // end computeFluxes

void updateCollisionalHeatFluxes(const domainGrid& Xgrid)
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

   matrix2D<double> kappae, kappai;
   matrix2D<double> dTedx, dTidx, dTedz, dTidz;
   matrix2D<double> dTedx_cex, dTidx_cex, dTedz_cez, dTidz_cez;
   matrix2D<double> taui, taue, Otaui, Otaue, delta, Otau2, OtauFactor, tauEff;

   kappae.initialize(nXcc,nZcc,0.0);
   kappai.initialize(nXcc,nZcc,0.0);
   
   dTedx.initialize(nXcc,nZcc,0.0);
   dTidx.initialize(nXcc,nZcc,0.0);
   dTedz.initialize(nXcc,nZcc,0.0);
   dTidz.initialize(nXcc,nZcc,0.0);

   dTedx_cex.initialize(nXce,nZcc,0.0);
   dTidx_cex.initialize(nXce,nZcc,0.0);
   dTedz_cez.initialize(nXcc,nZce,0.0);
   dTidz_cez.initialize(nXcc,nZce,0.0);
   
   taui.initialize(nXcc,nZcc,0.0);
   taue.initialize(nXcc,nZcc,0.0);
   Otaui.initialize(nXcc,nZcc,0.0);
   Otaue.initialize(nXcc,nZcc,0.0);
   delta.initialize(nXcc,nZcc,0.0);
   Otau2.initialize(nXcc,nZcc,0.0);
   OtauFactor.initialize(nXcc,nZcc,0.0);
   tauEff.initialize(nXcc,nZcc,0.0);

   
   //   define collision times and omega*taui for e and i
   //
   taui = taui0*pow(2.0*Ti,1.5)/N;
   taue = taui0*pow(2.0*Te,1.5)/N*sqrt(mM/2.0);
   Otaui = By/lambda0*taui;
   Otaue = By/lambda0*taue/mM;

   //   calculate perp ion heat fluxes
   //
   Xgrid.DDX(dTidx_cex,Ti);
   Xgrid.DDZ(dTidz_cez,Ti);
   Xgrid.communicate(dTidx_cex);
   
   Otau2 = Otaui*Otaui;
   delta = Otau2*Otau2 + 2.70*Otau2 + 0.677;
   OtauFactor = (2.0*Otau2 + 2.645)/delta;
   tauEff = taui*OtauFactor;

   kappai = Pi*(tauEff+tauiEff_min);
   Xgrid.InterpToCellEdges(kappai_cex, kappai, Ti, "C2", 0);
   Xgrid.InterpToCellEdges(kappai_cez, kappai, Ti, "C2", 1);
   Xgrid.communicate(kappai_cex);
   
   qix0 = -kappai_cex*dTidx_cex;
   qiz0 = -kappai_cez*dTidz_cez;

   if(modelDriftTerms) {
      
      //   calculate coperp ion heat fluxes
      //
      Xgrid.DDX(dTidx,Ti);
      Xgrid.DDZ(dTidz,Ti);
      Xgrid.communicate(dTidx);

      OtauFactor = Otaui*(gamma0/(gamma0-1.0)*Otau2 + 4.65)/delta;
      tauEff = taui*OtauFactor;
      //tauEff = min(tauEff,1.0);
  
      kappai_coperp = Pi*tauEff;
      Xgrid.InterpToCellEdges(qix_coperp,  kappai_coperp*dTidz*rcc, Ti, "C2", 0);
      Xgrid.InterpToCellEdges(qiz_coperp, -kappai_coperp*dTidx, Ti, "C2", 1);
      //setXminBoundary(qix_coperp, 0.0, 0.0);
      Xgrid.communicate(qix_coperp);
      
   }
   
   //   calculate perp electron heat fluxes
   //
   Xgrid.DDX(dTedx_cex,Te);
   Xgrid.DDZ(dTedz_cez,Te);
   
   Otau2 = Otaue*Otaue;
   delta = Otau2*Otau2 + 14.79*Otau2 + 3.7703;
   OtauFactor = (4.664*Otau2 + 11.92)/delta;
   tauEff = taue*OtauFactor;

   kappae = Pe*(tauEff/mM+taueEff_min);
   Xgrid.InterpToCellEdges(kappae_cex, kappae, Te, "C2", 0);
   Xgrid.InterpToCellEdges(kappae_cez, kappae, Te, "C2", 1);
   Xgrid.communicate(kappae_cex);
   
   qex0 = -kappae_cex*dTedx_cex;
   qez0 = -kappae_cez*dTedz_cez;

   if(modelDriftTerms) {
   
      //   calculate coperp electron heat fluxes
      //
      Xgrid.DDX(dTedx,Te);
      Xgrid.DDZ(dTedz,Te);
      Xgrid.communicate(dTedx);

      OtauFactor = Otaue*(gamma0/(gamma0-1.0)*Otau2 + 21.67)/delta;
      tauEff = taue*OtauFactor;
      //tauEff = min(tauEff,1.0*mM);

      kappae_coperp = Pe/mM*tauEff;
      Xgrid.InterpToCellEdges(qex_coperp, -kappae_coperp*dTedz*rcc, Te, "C2", 0);
      Xgrid.InterpToCellEdges(qez_coperp,  kappae_coperp*dTedx, Te, "C2", 1);
      //setXminBoundary(qex_coperp, 0.0, 0.0);
      Xgrid.communicate(qex_coperp);

   }
   
   //  update total collisional heat flux terms for heat eqs
   //  NOTE that coperp-x terms are already multiplied by r
   //  For some reason (I think it has to do with r=0 boundary)
   //  interpolating to cell edge and then multiplying by r doesnt work!
   //
   
   Flux_qi_x = qix*rce_x + qix_coperp;
   Flux_qi_z = qiz + qiz_coperp;
   Flux_qe_x = qex*rce_x + qex_coperp;
   Flux_qe_z = qez + qez_coperp;

   //divqe = 

}


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
   //double tauiVis0 = 1.0e-3;

   //   compute tensor values at cell-center
   //
   eta0 = 0.96*Pi*tauiVis0;
   eta1 = 0.96*Pi*tauiVis0;
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
   Qvis -= Vx*divPvis_x + Vz*divPvis_z;

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
   const int nZg = mesh->nZg;
 
   for (auto i=0; i<thisnX; i++) {
      for (auto j=0; j<nZg; j++) {
         var(i,j) = var(i,thisnZ-2*nZg+j);
         var(i,thisnZ-nZg+j) = var(i,nZg+j);
      }
   }

}


void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid, const int verbose)
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
   if(procID==0 && verbose) {
      cout << "dtSim = " << dtSim << endl;
      if(dtCFL_sound==dtmax)
      cout << "warning: dtSim set by CFL_sound"<< endl;
   }
}

