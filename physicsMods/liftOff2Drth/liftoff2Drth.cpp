/***
 * 
 * physics module for dpf liftoff 2Drth (electroTherm instability)
 * with no ion motion
 *
 * dBx/dt = cvac*dEz/dy/x 
 * dBy/dt = -dEz/dx
 * epsilonRel/cvac*dEz/dt = d(x*By)/dx/x - 4*pi/cvac*Jz
 * dEe/dt = J^2/sig - 2*me/Mi*N0/(gamma0-1)/taue*(Te-Ti0) - div(qe)
 * dqex/dt = delta_qe*(qex0 - qex)
 * dqey/dt = delta_qe*(qey0 - qey)
 *
 * Jz   = sig*Ez
 * Ee   = N0*Te/(gamma0-1)
 * qex  = -kappae*dTe/dx
 * qey  = -kappae*dTe/dy/x
 *
 * taue   = taue0*Te^1.5
 * sig    = sig0*Te^1.5
 * kappae = kappae0*Te^2.5
 *
 * Relaxation Const for Heat Flux Equations:
 * delta_qe = cvac^2/epsilonRel*N0/kappa0 
 *
 * gamma0  = 2/degFreedom + 1
 * taue0 [s]  = 3.44e4*Clog/10*Te0^1.5/N0 (Clog=10)
 * Ve0 [cm/s] = 4.19e-7*sqrt(Te0)
 * wpe0 [1/s] = 5.64e4*sqrt(N0)
 * sig0 [1/s] = wpe0^2*taue0/4/pi
 * kappae0 [1/s/cm] = 3.2*N0*Ve0^2*taue0
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

string geometry0;   // CAR or CYL
double gamma0;      // adiabatic coefficient
double kappae0, sig0;        // plasma conductivity [1/ns]
double taue0;       // ele collision time coefficient [ns]
double meMi;
double B0, B00;     // boundary value of magnetic field
matrix2D<double> Ex, Ey, Ez, Exold, Eyold, Ezold;       // electric field
matrix2D<double> Bx, By, Bz, Bxold, Byold, Bzold;       // magnetic field
matrix2D<double> Ee, Eeold, Te, SourceEe, divqe;        // electron energy density
matrix2D<double> qex, qey, qex0, qey0, qexold, qeyold;  // Heat flux terms
//matrix2D<double> eta, etace_x, etace_y, etace_xy;
matrix2D<double> taue, sig, sigce_x, sigce_y, sigce_xy; // conductivity terms
matrix2D<double> kappae, kappae_x, kappae_y;            // heat conductivity
matrix2D<double> Jx, Jy, Jz;                            // current density
matrix2D<double> etaJsq, Ex_cc, Ey_cc, Ez_cc, Stherm;   // electron heat source/sink
matrix2D<double> deltaTe;                               // initial perturbation profile

double N0, Ti0, Iscale_kA, dtIscale_ns;
double epsilonRel;

matrix2D<double> hy_cc, hy_cey, hy_ce; // y-dimension lame coefficient

// Set lower threshold values for density and temperature
//
double Tthresh=1.0e-1, Tmax = 1.0e4;
double cvac= 2.9979e1;   // speed of light [cm/ns]
double pi  = 3.14159265; // pi

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
   const int nYcc = Xgrid.nZcc;
   const int nYce = Xgrid.nZce;
   
   
   Bx.initialize(nXce,nYcc,0.0);
   By.initialize(nXcc,nYce,0.0);
   Bz.initialize(nXcc,nYcc,0.0);
   Bxold.initialize(nXce,nYcc,0.0);
   Byold.initialize(nXcc,nYce,0.0);
   Bzold.initialize(nXcc,nYcc,0.0);
   //
   Ex.initialize(nXcc,nYce,0.0);
   Ey.initialize(nXce,nYcc,0.0);
   Ez.initialize(nXce,nYce,0.0);
   Exold.initialize(nXcc,nYce,0.0);
   Eyold.initialize(nXce,nYcc,0.0);
   Ezold.initialize(nXce,nYce,0.0);
   
   Ee.initialize(nXcc,nYcc,0.0);
   Eeold.initialize(nXcc,nYcc,0.0);
   SourceEe.initialize(nXcc,nYcc,0.0);
   etaJsq.initialize(nXcc,nYcc,0.0);
   Stherm.initialize(nXcc,nYcc,0.0);
   
   Te.initialize(nXcc,nYcc,0.0);
   sig.initialize(nXcc,nYcc,0.0);
   kappae.initialize(nXcc,nYcc,0.0);
   deltaTe.initialize(nXcc,nYcc,0.0);
   taue.initialize(nXcc,nYcc,0.0);
   
   qex.initialize(nXce,nYcc,0.0);
   qey.initialize(nXcc,nYce,0.0);
   qex0.initialize(nXce,nYcc,0.0);
   qey0.initialize(nXcc,nYce,0.0);
   qexold.initialize(nXce,nYcc,0.0);
   qeyold.initialize(nXcc,nYce,0.0);

   Jx.initialize(nXcc,nYce,0.0);
   Jy.initialize(nXce,nYcc,0.0);
   Jz.initialize(nXce,nYce,0.0);
   
   Ex_cc.initialize(nXcc,nYcc,0.0);
   Ey_cc.initialize(nXcc,nYcc,0.0);
   Ez_cc.initialize(nXcc,nYcc,0.0);
   
   sigce_x.initialize(nXce,nYcc,0.0);
   sigce_y.initialize(nXcc,nYce,0.0);
   sigce_xy.initialize(nXce,nYce,0.0);
   kappae_x.initialize(nXce,nYcc,0.0);
   kappae_y.initialize(nXcc,nYce,0.0);
   //
   hy_cc.initialize(nXcc,nYcc,1.0);
   hy_cey.initialize(nXcc,nYce,1.0);
   hy_ce.initialize(nXce,nYcc,1.0);
   
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      Json::Value advScheme = Phys.get("advScheme",defValue);
      Json::Value geometry  = Phys.get("geometry",defValue);
      Json::Value gammaVal  = Phys.get("gammaC",defValue);
      
      //   get characteristic scales from input file
      // 
      Json::Value N0Val   = Phys.get("N0",defValue);
      Json::Value Ti0Val  = Phys.get("Ti0",defValue);
      Json::Value IscaleVal   = Phys.get("CurrScale_kA",defValue);
      Json::Value dtIscaleVal = Phys.get("dtCurrScale_ns",defValue);
      Json::Value epsilonRelVal = Phys.get("epsilonRel",defValue);
      //
      if(N0Val == defValue) {
         cout << "input ERROR: did not set N0 correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         N0 = N0Val.asDouble();
         if(procID==0) cout << "plasma density [1/cm^3] = " << N0 << endl;
      }
      //
      if(Ti0Val == defValue) {
         cout << "input ERROR: did not set Ti0 correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         Ti0 = Ti0Val.asDouble();
         if(procID==0) cout << "ion temperature [eV] = " << Ti0 << endl;
      }
      //
      if(IscaleVal == defValue) {
         cout << "input ERROR: did not set CurrScale_kA correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         Iscale_kA = IscaleVal.asDouble();
         if(procID==0) cout << "current scale [kA] = " << Iscale_kA << endl;
      }
      //
      if(dtIscaleVal == defValue) {
         cout << "input ERROR: did not set dtCurrScale_ns correctly" << endl;
	 exit (EXIT_FAILURE);
      } else {
         dtIscale_ns = dtIscaleVal.asDouble();
         if(procID==0) cout << "current rise time [ns] = " << dtIscale_ns << endl;
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

      //   set fundamental constants
      //
      //double kB  = 1.6022e-12; // Boltzmann constant [erg/eV]
      //double qe  = 1.6022e-19; // electron charge [C]       
      double me  = 9.1094e-28; // electron mass [g]
      double amu = 1.6605e-24; // atomic mass unit [g]
      double Mi  = 2.0*amu;    // ion mass [h] 
      double mu0 = pi*4.0e-7;  // permeability of free space [H/m]
      //double ep0 = 8.8542e-12; // permittivity of free space [F/m]
      meMi = me/Mi;

      //   calculate derived parameter scales
      //
      //N0  = 1.0e-14;   // density scale [1/cm^3]
      double Tescale = 1.0;      // electron temperature scale [eV]
      double Clog = 23.0 - 0.5*log(N0/Tescale/Tescale/Tescale);
      Clog = min(2.0,Clog);

      double wpescale  = 5.64e4*pow(N0,0.5)*1.0e-9; // ele plasma freq [rad/ns]
      taue0 = 3.44e4*pow(Tescale,1.5)/N0/1.0e-9; // collision time [ns] at Te=1eV
      double Lescale   = cvac/wpescale;    // ele inertial scale [cm]
      double Ve0 = 4.19*sqrt(Tescale)*1.0e-2; // Ve [cm/ns]
      sig0 = wpescale*wpescale/(4.0*pi)*taue0;
      kappae0 = 3.2*N0*Ve0*Ve0*taue0;


      if(procID==0) {
         cout << endl;
         cout << "derived scales:" << endl;
         cout << "ele plasma freq [rad/s] = " << wpescale*1.0e9 << endl; 
         cout << "ele collision time [s]  = " << taue0*1.0e-9 << endl; 
         cout << "ele inertial length [m] = " << Lescale << endl; 
      }


      //   calculate dimensionless parameters
      //
      double dyIscale = 2.0*pi*Xgrid.Xmin;
      B00 = mu0*(Iscale_kA*1.0e3)/(dyIscale/100.0)*1.0e4; // [Gauss]
      if(procID==0) {
	 cout << endl;
         cout << "simulation paramters:" << endl;
         cout << "conductivity  = " << sig0 << " nHz" << endl;
         cout << "heat conduct  = " << kappae0 << " 3.2*n0*Ve0^2*taue0" << endl;
         cout << "Iscale = " << Iscale_kA << " kA" << endl;      
         cout << "rise time = " << dtIscale_ns << " ns" << endl;      
      }


      if(gammaVal == defValue || geometry == defValue) {
         cout << "ERROR: geometry or gamma " << endl;
         cout << "not declared in input file" << endl;
         exit (EXIT_FAILURE);
      } 
      
      geometry0 = geometry.asString();
      if(geometry0=="CAR" || geometry0=="CYL") {
         if(procID==0) {
            cout << "geometry is " << geometry0 << endl;
         }
	 if(geometry0=="CYL") {
            for (auto i=0; i<nXcc; i++) {
               for (auto j=0; j<nYcc; j++) {
                  hy_cc(i,j) = Xgrid.Xcc.at(i);
                  if(i<nXce) hy_ce(i,j) = Xgrid.Xce.at(i);
               }
               for (auto j=0; j<nYce; j++) {
                  hy_cey(i,j) = Xgrid.Xcc.at(i);
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
      
   }
   else {
      cout << "value for key \"Physics\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
  
   //   get initial profiles for variables
   //
   const Json::Value Tevar = Phys.get("Te",defValue);
   if(Tevar.isObject()) { 
      Xgrid.setInitialProfile(Te,Tevar);
      if(procID==0) setXminBoundary(Te, 0.0, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(Te, 0.0, 1.0);   
      setZboundaryPeriodic(Te);
      Xgrid.communicate(Te);
   } else {
      cout << "value for Physics variable \"Te\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   //
   const Json::Value deltaTevar = Phys.get("deltaTe",defValue);
   if(deltaTevar.isObject()) { 
      Xgrid.setInitialProfile(deltaTe,deltaTevar);
      if(procID==0) setXminBoundary(deltaTe, 0.0, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(deltaTe, 0.0, 1.0);   
      setZboundaryPeriodic(deltaTe);
      Xgrid.communicate(deltaTe);
   } else {
      cout << "value for Physics variable \"deltaTe\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   //
   Te = Te*(1.0 + deltaTe);
   Ee = N0*Te/(gamma0-1.0);
   Eeold = Ee;

   //taue = taue0*pow(Te,1.5);
   //sig = sig0*pow(Te,1.5);
  
   // need inital flux calculation before add to output
   // and before first call to Physics.advance   
   computeFluxes(Xgrid, 1); // 1 for first order calculation   
  

   // add stuff to output files
   //
   dataFile.add(Jx, "Jx", 1);   // x-current density 
   dataFile.add(Jy, "Jy", 1);   // y-current density 
   dataFile.add(Jz, "Jz", 1);   // z-current density 
   dataFile.add(Ex, "Ex", 1);   // x-electric field
   dataFile.add(Ey, "Ey", 1);   // y-electric field
   dataFile.add(Ez, "Ez", 1);   // z-electric field
   dataFile.add(Bx, "Bx", 1);   // x-magnetic field
   dataFile.add(By, "By", 1);   // y-magnetic field 
   dataFile.add(Bz, "Bz", 1);   // z-magnetic field 
   dataFile.add(qex, "qex", 1);   // x-heat flux
   dataFile.add(qey, "qey", 1);   // y-heat flux
   dataFile.add(qex0, "qex0", 1);   // x-heat flux
   dataFile.add(qey0, "qey0", 1);   // y-heat flux
   dataFile.add(sig, "sig", 1);  // conductivity [1/ns]
   dataFile.add(kappae, "kappae", 1);  // conductivity [1/ns]
   dataFile.add(sigce_x, "sigce_x", 1);  // conductivity [1/ns]
   dataFile.add(sigce_y, "sigce_y", 1);  // conductivity [1/ns]
   dataFile.add(sigce_xy, "sigce_xy", 1);  // conductivity [1/ns]
   dataFile.add(SourceEe, "SourceEe", 1);   // ele heat source
   dataFile.add(etaJsq, "etaJsq", 1);   // joule heating
   dataFile.add(Stherm, "Stherm", 1);   // energy transfer to ions
   dataFile.add(Ex_cc, "Ex_cc", 1);   // z-electric field
   dataFile.add(Ey_cc, "Ey_cc", 1);   // z-electric field
   dataFile.add(Ez_cc, "Ez_cc", 1);   // z-electric field
   dataFile.add(Te, "Te", 1);   // ele temperature
   dataFile.add(Ee, "Ee", 1);    // electron energy density
   dataFile.add(taue, "taue", 1);  // electron collision time [ns]
   dataFile.add(gamma0,"gamma0",0); 
   //
   dataFile.add(Iscale_kA,"Iscale_kA",0);
   dataFile.add(N0,"N0",0);        // plasma density [1/cm^3]
   dataFile.add(meMi,"meMi",0);        // plasma density [1/cm^3]
   dataFile.add(Ti0, "Ti0", 0);    // ion temperature [eV]
   //
   dataFile.add(hy_cc,"hy_cc",0);
   dataFile.add(hy_ce,"hy_ce",0);

} // end Physics.initilize


void Physics::advance(const domainGrid& Xgrid, const double dt)
{
   const int nXcc = Xgrid.nXcc;
   const int nYcc = Xgrid.nZcc;
   const int nXg = Xgrid.nXg;
   const int nYg = Xgrid.nZg;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   //domainGrid* mesh = domainGrid::mesh;

   double thisdt;
   double Eemin = N0*Tthresh/(gamma0-1.0);
   double Eemax = N0*50.0/(gamma0-1.0);
   double dY = Xgrid.dZ; // Zgrid is actually Ygrid
   double dX = Xgrid.dX;
   const int Nsub=1;
   for (auto n=1; n<Nsub+1; n++) {
      thisdt = dt*n/Nsub;

      // update magnetic field and electron energy density from n to n+1
      //
      for (auto i=1; i<nXcc-1; i++) {
	 for (auto j=1; j<nYcc-1; j++) {

	    Bx(i,j) = Bxold(i,j) - cvac*thisdt*(Ez(i,j)-Ez(i,j-1))/dY/hy_ce(i,j);
	    By(i,j) = Byold(i,j) + cvac*thisdt*(Ez(i,j)-Ez(i-1,j))/dX;
	    Bz(i,j) = Bzold(i,j) - cvac*thisdt*(Ey(i,j)*hy_ce(i,j)-Ey(i-1,j)*hy_ce(i-1,j))/dX/hy_cc(i,j)
	    	                 + cvac*thisdt*(Ex(i,j)-Ex(i,j-1))/dY/hy_cc(i,j);
	    //
	    Ee(i,j) = Eeold(i,j) + thisdt*SourceEe(i,j) 
		                 - thisdt*(qex(i,j)*hy_ce(i,j)-qex(i-1,j)*hy_ce(i-1,j))/dX/hy_cc(i,j)
				 - thisdt*(qey(i,j)-qey(i,j-1))/dY/hy_cc(i,j);
	    if(Ee(i,j)<=Eemin) Ee(i,j) = Eemin;
	    if(Ee(i,j)>Eemax)  Ee(i,j) = Eemax;

	 }
      }
   
      // set magnetic field BCs at rmin
      //
      timeDomain* tmesh = timeDomain::tmesh;
      //cout << "thist = " << tmesh->tSim << endl;
      double thist = tmesh->tSim;
      B0 = B00*thist/dtIscale_ns; // By [Gauss] at lower boundary (Xmin-dX/2)
      if(thist>dtIscale_ns) B0=B00;

      if(procID==0) {
         setXminBoundary(Bx, 0.0, 1.0);   
         setXminBoundary(By, B0/hy_cc(nXg,1), 0.0);   
         setXminBoundary(Ee, 0.0, 1.0);   
         //setXminBoundary(Bz, 0.0, -1.0);   
      }
      if(procID==numProcs-1) {
         setXmaxBoundary(Ee, 0.0, 1.0);   
      }

      setZboundaryPeriodic(Bx);
      setZboundaryPeriodic(By);
      setZboundaryPeriodic(Bz);
      setZboundaryPeriodic(Ee);
      Xgrid.communicate(Bx);
      Xgrid.communicate(By);
      Xgrid.communicate(Bz);
      Xgrid.communicate(Ee);

      // compute fluxes using partialy updated fields at n+1/2 for 2nd order
      // or n+1 for 1st order
      computeFluxes(Xgrid, Nsub);


      // advance electric field from t^{n+1/2} to t^{n+3/2} 
      //
      double d0;
      for (auto i=nXg-1; i<nXcc-nXg; i++) {
	 for (auto j=nYg-1; j<nYcc-nYg; j++) {
            //d0 = 4.0*pi*sig0*pow(0.1,1.5)*thisdt/epsilonRel;	   
            //
            Ex(i,j) = Exold(i,j) + cvac*thisdt/epsilonRel*(Bz(i,j+1)-Bz(i,j))/dY/hy_cc(i,j);
            d0 = 4.0*pi*sigce_y(i,j)*thisdt/epsilonRel;	   
            Ex(i,j) /= (1.0 + d0); 
            //
            Ey(i,j) = Eyold(i,j) - cvac*thisdt/epsilonRel*(Bz(i+1,j)-Bz(i,j))/dX;
            d0 = 4.0*pi*sigce_x(i,j)*thisdt/epsilonRel;	   
            Ey(i,j) /= (1.0 + d0); 
            //
            Ez(i,j) = Ezold(i,j) + cvac*thisdt/epsilonRel*(By(i+1,j)*hy_cc(i+1,j)-By(i,j)*hy_cc(i,j))/dX/hy_ce(i,j)
                                 - cvac*thisdt/epsilonRel*(Bx(i,j+1)-Bx(i,j))/dY/hy_ce(i,j);
            d0 = 4.0*pi*sigce_xy(i,j)*thisdt/epsilonRel;	   
            Ez(i,j) /= (1.0 + d0); 

	    // update heat fluxes
	    //
	    d0 = cvac*cvac/epsilonRel*thisdt*N0/kappae_x(i,j);
	    qex(i,j) = (qexold(i,j) + d0*qex0(i,j))/(1.0 + d0); 
	    // 
	    d0 = cvac*cvac/epsilonRel*thisdt*N0/kappae_y(i,j);
	    qey(i,j) = (qeyold(i,j) + d0*qey0(i,j))/(1.0 + d0); 
	 }
      }

      if(procID==0) {
         //setXminBoundaryEz(Ez);   
         setXminBoundary(qex, 0.0, 1.0);   
      }
      if(procID==numProcs-1) {
         setXmaxBoundary(Ey, 0.0, 1.0);   
         setXmaxBoundary(Ez, 0.0, 1.0);   
         setXmaxBoundary(qex,0.0, 1.0);   
      }
      
      setZboundaryPeriodic(Ex);
      setZboundaryPeriodic(Ey);
      setZboundaryPeriodic(Ez);
      setZboundaryPeriodic(qex);
      setZboundaryPeriodic(qey);
      Xgrid.communicate(Ex);
      Xgrid.communicate(Ey);
      Xgrid.communicate(Ez);
      Xgrid.communicate(qex);
      Xgrid.communicate(qey);
   
   }

   // update old fields
   //
   Bxold = Bx;
   Byold = By;
   Bzold = Bz;
   Eeold = Ee;
   //
   Exold = Ex;
   Eyold = Ey;
   Ezold = Ez;
   //
   qexold = qex;
   qeyold = qey;
   
   // compute fluxes using fully updated fields at n+1
   //computeFluxes(Xgrid, 1);

} // end Physics.advance


void computeFluxes(const domainGrid& Xgrid, const int order)
{

   const int nXcc = Xgrid.nXcc;
   //const int nXce = Xgrid.nXce;
   const int nYcc = Xgrid.nZcc;
   //const int nYce = Xgrid.nZce;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   //matrix2D<double> Clog;
   //Clog.initialize(nXcc,nYcc,0.0);
  


   //  compute conductivity
   //
   Te = Ee/N0*(gamma0-1.0);
   //Clog = 23.0 - 0.5*log(N0/Te/Te/Te);
   //Clog = min(2.0,Clog);
   taue = taue0*pow(Te,1.5); //*10.0/Clog;
   sig  = sig0*pow(Te,1.5); //*10.0/Clog;
   kappae  = kappae0*pow(Te,2.5); //*10.0/Clog;

   /*
   double tauemin0, sigmin0;
   tauemin0 = taue0*pow(Ti0,1.5);
   sigmin0 = sig0*pow(Ti0,1.5);
   vector<double> tauemin, sigmin;
   tauemin.assign(nYcc,tauemin0);
   sigmin.assign(nYcc,sigmin0);
   if(procID==0) { // set insulator boundary values
      setXminBoundary(taue,tauemin);
      setXminBoundary(sig,sigmin);
   }
   */

   Xgrid.InterpToCellEdges(sigce_x,sig,sig,"C2",0);
   Xgrid.InterpToCellEdges(sigce_y,sig,sig,"C2",1);
   Xgrid.InterpToCellEdges(sigce_xy,sigce_x,sigce_x,"C2",1);
   Xgrid.communicate(sigce_x);
   Xgrid.communicate(sigce_y);
   Xgrid.communicate(sigce_xy);
   
   Xgrid.InterpToCellEdges(kappae_x,kappae,kappae,"C2",0);
   Xgrid.InterpToCellEdges(kappae_y,kappae,kappae,"C2",1);
   Xgrid.DDX(qex0,Te);
   Xgrid.DDZ(qey0,Te);
   Xgrid.communicate(qex0);
   Xgrid.communicate(qey0);
   qex0 *= kappae_x;
   qex0 *= -1.0;
   qey0 *= kappae_y/hy_cey; // need to divide by r
   qey0 *= -1.0;

   
   // compute eta*J^2 at cell center
   //
   Xgrid.InterpToCellCenter(Ex_cc,Ex);
   Xgrid.InterpToCellCenter(Ey_cc,Ey);
   Xgrid.InterpToCellCenter(Ez_cc,Ez);
   Xgrid.communicate(Ex_cc);
   Xgrid.communicate(Ey_cc);
   Xgrid.communicate(Ez_cc);
   setZboundaryPeriodic(Ey_cc);
   setZboundaryPeriodic(Ez_cc);


   etaJsq = sig*(Ex_cc*Ex_cc + Ey_cc*Ey_cc + Ez_cc*Ez_cc);
   etaJsq /= 1.6022e-12; // convert to correct units [eV/cm3/ns]
   Stherm = 2.0*meMi/taue/(gamma0-1.0)*N0*(Te - Ti0);
   SourceEe = etaJsq - Stherm;

   //Jx = sigce_y*Ex;
   //Jy = sigce_x*Ey;
   //Jz = sigce_xy*Ez;
   

}


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

   //domainGrid* mesh = domainGrid::mesh;
   //const int ishift = mesh->nXg;

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
   const int thisnY = var.size1();

   domainGrid* mesh = domainGrid::mesh;
   const int jshift = mesh->nZg;
   const int nYcc = mesh->nZcc;
   const int nYce = mesh->nZce;

   assert(jshift==2 || jshift==1);

   if(thisnY==nYcc) {

   if(jshift==2) {
      for (auto i=0; i<thisnX; i++) {
	 var(i,1) = var(i,thisnY-3);
	 var(i,0) = var(i,thisnY-4);

	 var(i,thisnY-2) = var(i,2);
	 var(i,thisnY-1) = var(i,3);
      }
   }
   else {
      for (auto i=0; i<thisnX; i++) {
	 var(i,0) = var(i,thisnY-2);
	 var(i,thisnY-1) = var(i,1);
      }
   }

   }

   if(thisnY==nYce) {

   if(jshift==2) {
      for (auto i=0; i<thisnX; i++) {
	 var(i,0) = var(i,thisnY-3);
	 var(i,thisnY-1) = var(i,2);
      }
   }
   else {
      for (auto i=0; i<thisnX; i++) {
	 //var(i,0) = var(i,thisnZ-2);
	 //var(i,thisnZ-1) = var(i,1);
      }
   }

   }

}

void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid, const int verbose)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   /*
   matrix2D<double> Cchar_x(Cs), Cchar_z(Cs);
   double Cmax_x, Cmax_z, Cmax;
   Cchar_x += abs(Vx);
   Cchar_z += abs(Vz);
   Cmax_x = max(Cchar_x);
   Cmax_z = max(Cchar_z);
   Cmax = max(Cmax_x,Cmax_z);
   //cout << "Cmax = " << Cmax << endl;
   //cout << "abs(V) = " << max(abs(V))<< endl;
   */
   matrix2D<double> nuSourceEe;
   nuSourceEe.initialize(Xgrid.nXcc,Xgrid.nZcc,0.0);
   nuSourceEe = 1.0+abs(SourceEe)/(gamma0-1.0)/(N0*Te);


   const double dX = Xgrid.dX;
   const double dY = Xgrid.dZ;

   double dtCFL_light = dX*dY/(dX+dY)*sqrt(epsilonRel)/cvac;
   double dtmax = dtCFL_light;

   double dtSourceEe = 1.0/max(nuSourceEe);

   dtSim = min(dtmax/tDom.dtFrac,dtSourceEe/tDom.dtFrac);
   dtSim = min(dtSim,tDom.dtOut);


   double sigmax = max(sig);
   if(procID==0 && verbose) {
      cout << "4*pi*sig*dt/epsilonRel = " << 4.0*pi*sigmax*dtSim/epsilonRel << endl;
      cout << "dtCFL_light = " << dtCFL_light << endl;
      cout << "dtSourceEe  = " << dtSourceEe << endl;
      cout << "dtSim = " << dtSim << endl;
   }

}



