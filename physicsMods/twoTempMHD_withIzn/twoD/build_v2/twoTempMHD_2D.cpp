/***
 * 
 * physics module for 2D two Temp MHD dpf rundown
 * and lift-off with ionization physics included
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
 * rhoScale = Mn*Nscale [kg/m^3], Mn=Amass*amu [kg]
 * Vscale   = sqrt(Pscale/rhoScale) [m/s]
 * tscale   = Xscale/Vscale [s]
 * Bscale   = sqrt(mu0*Pscale) [tesla]
 * Jscale   = Bscale/Xscale/mu0 [Amps/m^2]
 * Escale   = Vscale*Bscale [V/m]
 *
 * Governing Equations:
 * dN/dt  + d(Mx)/dx = 0
 * dMx/dt + d(Mx*Vx + P)/dx = -Jz*By
 * dEi/dt + d((Ei+Pi)*Vx)/dx = NeUdotE + Qie 
 * dEe/dt + d((Ee+Pe)*Vx)/dx = Jz*Ez  - Qie - NeUdotE + SEe
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
 * P = Pn + Pi + Pe
 * Ta = Pn/Nn = Pi/Ne
 * Te = Pe/Ne
 * eta = eta0/Te^1.5
 * NeUdotE = -Ux*(Jz*By + dPe/dx) (Ex from force balance)
 * Qie = 3.0*me/Mi*Ne/tauei*(Te-Ta)
 * Qne = 3.0*me/Mn*Ne/tauen*(Te-Ta)
 * taue = taue0*Te^1.5/Ne
 * Se = Ne*Nn*kiz
 * SEe = -Uizn*Se
 *
 * Dimensionless parameters:
 * delta0 = (V0/cvac)^2
 * Li0 = Li0/r0, Li0 = cvac/wpi0
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

#include "json/json.h"
#include "vectorMath.h"
#include "matrix2D.h"
#include "domainGrid.h"
#include "timeDomain.h"
#include "HDF5dataFile.h"

#include "Physics.h"


using namespace std;


domainGrid m_Xgrid;
double m_dX, m_dZ, m_nXce, m_nZce, m_nXcc, m_nZcc;
int m_nXg, m_nZg;

string advScheme0, advSchemeHall0;  // advection differencing scheme
string XlowBC, XhiBC;  // boundary condition strings
string geometry0;   // CAR or CYL
bool useEleScaleCorrections = false;
bool useIonScaleCorrections = false;
double gamma0;      // adiabatic coefficient
double eta0;        // resistivity coefficient
double taue0;       // ele collision time coefficient
double taui0;       // ion collision time coefficient
double Ve0sq;       // square of electron therm speed norm
double mM;          // me/Mi
double tauiMax, taueMax;  // maximum collision times for viscosity
double delta0;      // relaxation const (V0/cvac)^2
double Li0;         // normalized ion skin depth
double Le0sq;       // normalized electron skin depth squared
double B0, B00;     // boundary value of magnetic field
int Nsub;           // time-solver subcycle steps
matrix2D<double> N, Mx, Mz, Ei, Ee, By, Ez, Ex, Jz, Jx;   // time-evolving variables
matrix2D<double> eta, Cs, Cspeed_x, Cspeed_z, Vx, Vz, P, Pe, Pi, Te, Ti;
matrix2D<double> Jx0, Jz0, Jx0stagTime, Jz0stagTime;
matrix2D<double> J0x, J0x_x, J0x_z, J0z, J0z_x, J0z_z;
matrix2D<double> Jzcc, Ezcc, Jxcc, Excc, Eidealz, Eidealx;
matrix2D<double> eta_x, Ne_x, N_x, eta_z, Ne_z, N_z, Qvisc;
matrix2D<double> qex, qexold, qex0, kappae_x;  // ele heat flux
matrix2D<double> qez, qezold, qez0, kappae_z;
matrix2D<double> chie_perp, kappae, divqe;     // ele heat flux
matrix2D<double> Nold, Mxold, Mzold, Eiold, Eeold, Byold, Ezold, Exold, Jzold, Jxold;
matrix2D<double> FluxRatio_x, FluxLim_x, FluxR_x, FluxL_x;
matrix2D<double> FluxRatio_z, FluxLim_z, FluxR_z, FluxL_z;  
matrix2D<double> FluxN_x, FluxMx_x, FluxMz_x, FluxEi_x, FluxEe_x;
matrix2D<double> FluxN_z, FluxMx_z, FluxMz_z, FluxEi_z, FluxEe_z;
matrix2D<double> Qie, NeUdotE, JdotE, JcrossBx, JcrossBz, taue, taui, nue_spi, nue_vac;
matrix2D<double> deltaN;

// terms for finite ion and electron inertial length corrections
//
matrix2D<double> Vex, Vez, Vhallx, Vhallz, Vhallx_z, Vhallz_x;
matrix2D<double> Ehallx, Ehallz, Egradx, Egradz; 
matrix2D<double> Pe_eff, curlJ0, curlVe;
matrix2D<double> eJxFlux_x, eJzFlux_x, eJxFlux_z, eJzFlux_z;
matrix2D<double> divJxstress, divJzstress;


// ion viscousity terms
//
matrix2D<double> Pii_xxold, Pii_xx0, Pii_xx;
matrix2D<double> Pii_xzold, Pii_xz0, Pii_xz;
matrix2D<double> Pii_zzold, Pii_zz0, Pii_zz;
matrix2D<double> Pii_yy;
matrix2D<double> divPii_x, divPii_z, divUidotPii;
matrix2D<double> etaVis_ion, etaVis_ion_x, etaVis_ion_z, UidotPii_x, UidotPii_z;


// ele viscousity terms
//
matrix2D<double> Pie_xxold, Pie_xx0, Pie_xx;
matrix2D<double> Pie_xzold, Pie_xz0, Pie_xz;
matrix2D<double> Pie_zzold, Pie_zz0, Pie_zz;
matrix2D<double> Pie_yy;
matrix2D<double> divPiex, divPiez, divPiex_z, divPiez_x, divUedotPie;
matrix2D<double> etaVis_ele, etaVis_ele_x, etaVis_ele_z, UedotPie_x, UedotPie_z;


// new stuff for ionization
//
const double Uizn = 13.6;
const double g0 = 2, g1=1;
double Zmin = 1.0;
matrix2D<double> Zbar, Neold, Se, SEe, Ne, Nn, nue_neu, nue_izn;
matrix2D<double> FluxNe_x, FluxNe_z;

double Nscale, Xscale, Amass, Iscale, dtIscale;
double Pscale, Tscale, Bscale, Jscale, Ezscale, Vscale, tscale, Mi, Mn;
double wcescale, epsilonRel, meRel;

matrix2D<double> hy_cc, hy_ce; // y-dimension lame coefficient

// Set lower threshold values for N, T, and S
double Nthresh=1.0e-4, Tthresh=2.5e-2, Ethresh, Pthresh;

// Set conditions for vaccum resistivity model
double NvacC=0.01, NvacP = 4;

void initializeMemberVars(const domainGrid& );
void computeFluxes(const domainGrid&, const int);
void computeHeatFluxes(const domainGrid&);
void computeIonViscStressTensor(const domainGrid&);
void computeEleViscStressTensor(const domainGrid&);
void computeStrainTensor( matrix2D<double>&, 
                          matrix2D<double>&,
                          matrix2D<double>&,
                          matrix2D<double>&,
                    const matrix2D<double>&,
                    const matrix2D<double>&,
                    const domainGrid&        );
void computeStrainTensorStag( matrix2D<double>&, 
                              matrix2D<double>&,
                              matrix2D<double>&,
                              matrix2D<double>&,
                        const matrix2D<double>&,
                        const matrix2D<double>& );
void computeJ0(matrix2D<double>&, matrix2D<double>&, 
         const matrix2D<double>&, const domainGrid& );
void setJ0( const domainGrid& );
void updatePhysicalVars(const domainGrid& );
void updateJccAndEcc( const domainGrid& );
void updateCollisionTerms(const domainGrid& );
void computeDivOfOhmsLawStressTensor( const domainGrid& );
void computeDivOfElectronViscosity( const domainGrid& );
void computeIdealEatEdges( const domainGrid&, 
                           const matrix2D<double>&,
                           const matrix2D<double>& );
void computeHallEatEdges( const domainGrid& );
void computeGradEatEdges( const domainGrid& );
void computeJouleHeatingTerms(const domainGrid& );
void setNeutralInteractionRates(const domainGrid&,     const matrix2D<double>&,
                                const matrix2D<double>&, const matrix2D<double>& );
void advanceSemiImplicitVars(const domainGrid&, const int, const double, const double);
void advanceExplicitVars(const int, const double, const double);
void advanceElectronViscFluxes(const int, const double, const double);

void setXminExtrap(matrix2D<double>&);

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

   m_Xgrid = Xgrid;
   m_nXg   = Xgrid.nXg;
   m_nZg   = Xgrid.nZg;
   m_nXcc  = Xgrid.nXcc;
   m_nXce  = Xgrid.nXce;
   m_nZcc  = Xgrid.nZcc;
   m_nZce  = Xgrid.nZce;
   m_dX    = Xgrid.dX;
   m_dZ    = Xgrid.dZ;


   //domainGrid* mesh = domainGrid::mesh;
   const int nXg  = Xgrid.nXg;
   const int nXcc = Xgrid.Xcc.size();
   const int nXce = Xgrid.Xce.size();
   
   initializeMemberVars(Xgrid);

   parseInputFile(Xgrid,root);

   // define remaining vars from those specified in input file
   //
   Ne = Zbar*N;
   Nn = N-Ne; 
   Mx  = N*Vx;
   Mz  = N*Vz;
   Pe = Ne*Te;
   Pi = N*Ti;
   Ei = 0.5*(Mx*Mx + Mz*Mz)/N + Pi/(gamma0-1.0);
   Ee = Pe/(gamma0-1.0);

   // initialize current density from magnetic field
   //
   computeJ0(Jx0,Jz0,By,Xgrid);
   Jx = Jx0;
   Jz = Jz0;
   Jx0stagTime = Jx;
   Jz0stagTime = Jz;
   
   Xgrid.DDX(Jzcc,hy_cc*By);
   Jzcc = Jzcc/hy_cc; 
   Xgrid.InterpToCellCenter(Jxcc,Jx0);
   Xgrid.communicate(Jxcc);
   Xgrid.communicate(Jzcc);


   // need a complete RHS eval before writing at t=0
   // and before first call to Physics.advance   

   updatePhysicalVars(Xgrid);
   
   updateCollisionTerms(Xgrid);  
   
   computeHeatFluxes(Xgrid);
   qex = qex0;
   qez = qez0;

   updateJccAndEcc( Xgrid );
   
   if(tauiMax>0.0) {
      computeIonViscStressTensor(Xgrid);  
      Pii_xx = Pii_xx0;  
      Pii_xz = Pii_xz0;  
      Pii_zz = Pii_zz0;  
   }
   if(taueMax>0.0) {
      computeEleViscStressTensor(Xgrid);  
      Pie_xx = Pie_xx0;  
      Pie_xz = Pie_xz0;  
      Pie_zz = Pie_zz0;  
      computeDivOfElectronViscosity(Xgrid);
   }

   computeIdealEatEdges(Xgrid,abs(Vx),By);
   
   if(useIonScaleCorrections) {
      computeHallEatEdges(Xgrid);
      computeGradEatEdges(Xgrid);
      computeDivOfOhmsLawStressTensor(Xgrid);
   }

   Ez = eta_x*Jz + Eidealz; // initialize electric field from Ohm's law 
   Ex = eta_z*Jx + Eidealx; // initialize electric field from Ohm's law 

   computeJouleHeatingTerms(Xgrid);
   
   computeFluxes(Xgrid, 1); // 1 for first order calculation   

   // dont forget to initialize all the old values
   //
   Nold  = N;
   Mxold  = Mx;  
   Mzold  = Mz;  
   Eiold = Ei;
   Eeold = Ee;
   Byold  = By;
   Neold = Ne;
   //
   Ezold = Ez;
   Exold = Ex;
   Jzold = Jz;
   Jxold = Jx;
   //
   qexold = qex;
   qezold = qez;
   //
   Pii_xxold = Pii_xx;
   Pii_xzold = Pii_xz;
   Pii_zzold = Pii_zz;
   //
   Pie_xxold = Pie_xx;
   Pie_xzold = Pie_xz;
   Pie_zzold = Pie_zz;

   // add selected members to output files
   //
   addMembersToDataFile( dataFile );

} // end Physics.initilize


void Physics::advance(const domainGrid& Xgrid, const double a_dt)
{
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg = Xgrid.nXg;
   const int nZg = Xgrid.nZg;
   const double dX = Xgrid.dX;
   const double dZ = Xgrid.dZ;
   
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
      thist = tmesh->tSim + thisdt; // time at end of step
      dt_stag = a_dt; 
      if(thist==thisdt && n==1) dt_stag = thisdt;

      advanceExplicitVars(n, thisdt, thist);

      updatePhysicalVars(Xgrid);
      
      updateCollisionTerms(Xgrid);  
     
      if(n==2) { // for du-Fort and Frenkel, only update at stage=2
         computeHeatFluxes(Xgrid);  
         if(tauiMax>0.0) computeIonViscStressTensor(Xgrid);  
         if(taueMax>0.0) computeEleViscStressTensor(Xgrid);  
      }
      
      computeIdealEatEdges(Xgrid,abs(Vx),By);
      
      if(useIonScaleCorrections) {
         computeHallEatEdges(Xgrid);
         computeGradEatEdges(Xgrid);
         computeDivOfOhmsLawStressTensor(Xgrid);
      }
      
      if(taueMax>0.0) {
         advanceElectronViscFluxes(n,a_dt,thist);
         //advanceElectronViscFluxes(n,thisdt,thist);
         computeDivOfElectronViscosity(Xgrid);
      }

      advanceSemiImplicitVars(Xgrid,n,dt_stag,thist);

      updateJccAndEcc( Xgrid );
   
      computeJouleHeatingTerms(Xgrid);
      
      computeFluxes(Xgrid, Nsub);

   } // finish subcycle steps

} // end Physics.advance


void advanceExplicitVars(const int stage, const double thisdt, const double thist)
{
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
      
   for (auto i=m_nXg; i<m_nXcc-m_nXg; i++) {
      for (auto j=m_nZg; j<m_nZcc-m_nZg; j++) {

	 N(i,j)  = Nold(i,j)  - thisdt*(FluxN_x(i,j) -FluxN_x(i-1,j))/hy_cc(i,j)/m_dX
	                      - thisdt*(FluxN_z(i,j) -FluxN_z(i,j-1))/m_dZ;
         Mx(i,j) = Mxold(i,j) - thisdt*(FluxMx_x(i,j)-FluxMx_x(i-1,j))/hy_cc(i,j)/m_dX
	                      - thisdt*(FluxMx_z(i,j)-FluxMx_z(i,j-1))/m_dZ
	                      - thisdt*divPii_x(i,j)
                              + thisdt*JcrossBx(i,j);
	 if(geometry0=="CYL") Mx(i,j) = Mx(i,j) + thisdt*P(i,j)/hy_cc(i,j);
         Mz(i,j) = Mzold(i,j) - thisdt*(FluxMz_x(i,j)-FluxMz_x(i-1,j))/hy_cc(i,j)/m_dX
                              - thisdt*(FluxMz_z(i,j)-FluxMz_z(i,j-1))/m_dZ
	                      - thisdt*divPii_z(i,j)
                              + thisdt*JcrossBz(i,j);
         Ei(i,j) = Eiold(i,j) - thisdt*(FluxEi_x(i,j)-FluxEi_x(i-1,j))/hy_cc(i,j)/m_dX
                              - thisdt*(FluxEi_z(i,j)-FluxEi_z(i,j-1))/m_dZ
	                      - thisdt*divUidotPii(i,j)
      		              + thisdt*(NeUdotE(i,j) + Qie(i,j));
         Ee(i,j) = Eeold(i,j) - thisdt*(FluxEe_x(i,j)-FluxEe_x(i-1,j))/hy_cc(i,j)/m_dX
                              - thisdt*(FluxEe_z(i,j)-FluxEe_z(i,j-1))/m_dZ
      		              - thisdt*divqe(i,j)
      		              + thisdt*(JdotE(i,j) - Qie(i,j) - NeUdotE(i,j) + SEe(i,j));
	 By(i,j) = Byold(i,j) + thisdt*(Ez(i,j)-Ez(i-1,j))/m_dX
                              - thisdt*(Ex(i,j)-Ex(i,j-1))/m_dZ;
	 
         Ne(i,j) = Neold(i,j) - thisdt*(FluxNe_x(i,j)-FluxNe_x(i-1,j))/hy_cc(i,j)/m_dX
                              - thisdt*(FluxNe_z(i,j)-FluxNe_z(i,j-1))/m_dZ
                              + thisdt*Se(i,j);
	 
         // check thresholds
         //
         if(N(i,j)<=Nthresh) N(i,j) = Nthresh;
	 if(N(i,j)!=N(i,j)) {
            cout << "bout to go bad: N(i,j) = " << N(i,j) << endl;
            cout << "thist = " << thist << endl;
            cout << "X(i) = " << m_Xgrid.Xcc.at(i) << endl;
            cout << "Z(j) = " << m_Xgrid.Zcc.at(j) << endl;
	    exit (EXIT_FAILURE);
	 }
         if(Ne(i,j)<N(i,j)*Zmin) Ne(i,j) = N(i,j)*Zmin;
	 if(Ne(i,j)>N(i,j))      Ne(i,j) = N(i,j);
	 
         Ethresh = N(i,j)/Nthresh*Pthresh/(gamma0-1.0) 
                 + 0.5*(Mx(i,j)*Mx(i,j)+Mz(i,j)*Mz(i,j))/N(i,j);
	 if(Ei(i,j)<=Ethresh) Ei(i,j) = Ethresh;
	 Ethresh = Pthresh/(gamma0-1.0)*Ne(i,j)/Nthresh;
	 if(Ee(i,j)<=Ethresh) Ee(i,j) = Ethresh;

      }
   }
      
   //  apply boundary conditions and communicate
   //
   B0 = thist/dtIscale*B00;
   if(B0>B00) B0 = B00;
   //if(thist>=0.03) B0 = B00*pow(0.03/thist,2);

   if(procID==0) {
      m_Xgrid.setXminBoundary(N, 0.0, 1.0);   
      m_Xgrid.setXminBoundary(Ne, 0.0, 1.0);   
      m_Xgrid.setXminBoundary(Mx, 0.0, -1.0);   
      m_Xgrid.setXminBoundary(Mz, 0.0, 0.0);   
      m_Xgrid.setXminBoundary(Ei, 0.0, 1.0);
      m_Xgrid.setXminBoundary(Ee, 0.0, 1.0);
      //if(N(0,j)<Nthresh || N(1,j)<Nthresh) {
      //   N(0,j) = Nthresh;
      //   N(1,j) = Nthresh;
      //}
         
      if(XlowBC.compare("axis") == 0 || XlowBC.compare("symmetry") == 0) {
         m_Xgrid.setXminBoundary(By, 0.0, -1.0);   
      } 
      else if(XlowBC.compare("insulator") == 0) {
         m_Xgrid.setXminBoundary_J(By,B0,0.0);   
      } 
      else {
         cout << "low MagFieldBC not defined !!!! " << endl;
      }
   }
   
   if(procID==numProcs-1) {
      m_Xgrid.setXmaxBoundary(N, 0.0, 1.0);   
      m_Xgrid.setXmaxBoundary(Ne, 0.0, 1.0);   
      m_Xgrid.setXmaxBoundary(Mx, 0.0, -1.0);   
      m_Xgrid.setXmaxBoundary(Mz, 0.0, 0.0);   
      m_Xgrid.setXmaxBoundary(Ei, 0.0, 1.0);   
      m_Xgrid.setXmaxBoundary(Ee, 0.0, 1.0); 
         
      if(XhiBC.compare("conductor") ==0 ) {
         m_Xgrid.setXmaxBoundary_J(By,0,1);   
      } 
      else if(XhiBC.compare("insulator") == 0) {
         //cout << "hy_cc.at(nXcc-nXg) = " << hy_cc.at(nXcc-nXg) << endl;	 
         m_Xgrid.setXmaxBoundary_J(By,B0,0); 
      } 
      else {
         cout << "Xhi MagFieldBC not defined !!!! " << endl;
      }
   }      
      
   m_Xgrid.setZboundaryPeriodic(N); 
   m_Xgrid.setZboundaryPeriodic(Ne); 
   m_Xgrid.setZboundaryPeriodic(Mx); 
   m_Xgrid.setZboundaryPeriodic(Mz); 
   m_Xgrid.setZboundaryPeriodic(Ei); 
   m_Xgrid.setZboundaryPeriodic(Ee); 
   m_Xgrid.setZboundaryPeriodic(By); 
   m_Xgrid.communicate(N);
   m_Xgrid.communicate(Ne);
   m_Xgrid.communicate(Mx);
   m_Xgrid.communicate(Mz);
   m_Xgrid.communicate(Ei);
   m_Xgrid.communicate(Ee);
   m_Xgrid.communicate(By);
   
   if(stage==2) {
      // update old values for explicit vars
      //
      Nold = N;
      Mxold = Mx;
      Mzold = Mz;
      Eiold = Ei;
      Eeold = Ee;
      Byold = By;
      Neold = Ne;
   }

}

void advanceElectronViscFluxes(const int stage, const double thisdt, const double thist)
{
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   //if(!procID) cout << "dt_stag = " << thisdt << endl;
   //if(!procID) cout << "thist   = " << thist << endl;
   //if(!procID) cout << "stage   = " << stage << endl;
   
   
   double d0Pie;
   for (auto j=m_nZg-1; j<m_nZcc-m_nZg; j++) {
      for (auto i=m_nXg-1; i<m_nXcc-m_nXg; i++) {
         d0Pie = thisdt/delta0*Ne_x(i,j)/etaVis_ele_x(i,j)*mM;
         Pie_xx(i,j) = (Pie_xxold(i,j) + d0Pie*Pie_xx0(i,j))/(1.0 + d0Pie);
         //
         d0Pie = thisdt/delta0*Ne(i,j)/etaVis_ele(i,j)*mM;
         Pie_xz(i,j) = (Pie_xzold(i,j) + d0Pie*Pie_xz0(i,j))/(1.0 + d0Pie);
         //
         d0Pie = thisdt/delta0*Ne_z(i,j)/etaVis_ele_z(i,j)*mM;
         Pie_zz(i,j) = (Pie_zzold(i,j) + d0Pie*Pie_zz0(i,j))/(1.0 + d0Pie);
      }
   }
   if(procID==0) {
      m_Xgrid.setXminBoundary(Pie_xz, 0.0, -1.0);   
   }
   if(procID==numProcs-1) {
      m_Xgrid.setXmaxBoundary(Pie_xz, 0.0, -1.0);   
   }  
   m_Xgrid.setZboundaryPeriodic(Pie_xx); 
   m_Xgrid.setZboundaryPeriodic(Pie_xz); 
   m_Xgrid.setZboundaryPeriodic(Pie_zz); 
   m_Xgrid.communicate(Pie_xx); 
   m_Xgrid.communicate(Pie_xz); 
   m_Xgrid.communicate(Pie_zz); 
   if(stage==2) {
      Pie_xxold = Pie_xx;
      Pie_xzold = Pie_xz;
      Pie_zzold = Pie_zz;
   }
   
   //Pie_xx = Pie_xx0;
   //Pie_zz = Pie_zz0;
   //Pie_xz = Pie_xz0;

}

void advanceSemiImplicitVars(const domainGrid& Xgrid, const int stage, 
                             const double thisdt, const double thist)
{
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg = Xgrid.nXg;
   const int nZg = Xgrid.nZg;
   
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
      computeJ0(Jx0stagTime,Jz0stagTime,By,Xgrid);
   }

   // advance electric field and current density
   //
   double d0, e0_x, e0_z, E0z, E0x;
   double d0qex, d0qez, d0Pii;
   d0 = delta0/thisdt;	   
   for (auto j=nZg; j<nZcc-nZg; j++) {
      for (auto i=nXg-1; i<nXcc-nXg; i++) {
         E0z  = Eidealz(i,j) + Ehallz(i,j) + Egradz(i,j);
         E0z += divJzstress(i,j)/Ne_x(i,j);
         E0z -= Li0*divPiez_x(i,j)/Ne_x(i,j);
         e0_x = Le0sq/thisdt/Ne_x(i,j);	   
         //
         Ez(i,j)  = Jz0stagTime(i,j) + d0*Ezold(i,j)
                  - (e0_x*Jzold(i,j) - E0z)/(e0_x + eta_x(i,j));
         Ez(i,j) /= d0 + 1.0/(e0_x + eta_x(i,j)); 
         //
         Jz(i,j)  = e0_x*Jzold(i,j) + Ez(i,j) - E0z;
         Jz(i,j) /= e0_x + eta_x(i,j);
         //
         d0qex = thisdt/delta0*Ne_x(i,j)/kappae_x(i,j);
         qex(i,j) = (qexold(i,j) + d0qex*qex0(i,j))/(1.0 + d0qex);
         //
         d0Pii = thisdt/delta0*N_x(i,j)/etaVis_ion_x(i,j);
         Pii_xx(i,j) = (Pii_xxold(i,j) + d0Pii*Pii_xx0(i,j))/(1.0 + d0Pii);
         //
         d0Pii = thisdt/delta0*N(i,j)/etaVis_ion(i,j);
         Pii_xz(i,j) = (Pii_xzold(i,j) + d0Pii*Pii_xz0(i,j))/(1.0 + d0Pii);
      }
   }
   //Pii_xz = Pii_xz0; // this guy cell-center...need BC's

   for (auto j=nZg-1; j<nZcc-nZg; j++) {
      for (auto i=nXg; i<nXcc-nXg; i++) {
         E0x  = Eidealx(i,j) + Ehallx(i,j) + Egradx(i,j);
         E0x += divJxstress(i,j)/Ne_z(i,j);
         E0x -= divPiex_z(i,j)/Ne_z(i,j);
         e0_z = Le0sq/thisdt/Ne_z(i,j);	   
         //
         Ex(i,j)  = Jx0stagTime(i,j) + d0*Exold(i,j)
                  - (e0_z*Jxold(i,j) - E0x)/(e0_z + eta_z(i,j));
         Ex(i,j) /= d0 + 1.0/(e0_z + eta_z(i,j)); 
         //
         Jx(i,j)  = e0_z*Jxold(i,j) + Ex(i,j) - E0x;
         Jx(i,j) /= e0_z + eta_z(i,j);
         //
         d0qez = thisdt/delta0*Ne_z(i,j)/kappae_z(i,j);
         qez(i,j)  = (qezold(i,j) + d0qez*qez0(i,j))/(1.0 + d0qez);
         //
         d0Pii = thisdt/delta0*N_z(i,j)/etaVis_ion_z(i,j);
         Pii_zz(i,j) = (Pii_zzold(i,j) + d0Pii*Pii_zz0(i,j))/(1.0 + d0Pii);
      }
   }
   //cout << "Xce(nXg-1) =" << Xgrid.Xce.at(nXg-1) << endl;
   //cout << "Xce(nXg-1) =" << Xgrid.Xce.at(nXg-1) << endl;
   //cout << "Xce(nXcc-Xg-1) =" << Xgrid.Xce.at(nXcc-nXg-1) << endl;
   if(procID==0) {
      Xgrid.setXminBoundary(Pii_xz, 0.0, -1.0);   
   }
   if(procID==numProcs-1) {
      for (auto j=0; j<Ez.size1(); j++) {
         Ez(nXcc-nXg,j) = Ez(nXcc-nXg-1,j);
      }
      Xgrid.setXminBoundary(Pii_xz, 0.0, -1.0);   
   }  
   Xgrid.setZboundaryPeriodic(Ez); 
   Xgrid.setZboundaryPeriodic(Jz); 
   Xgrid.setZboundaryPeriodic(qez); 
   Xgrid.setZboundaryPeriodic(Ex); 
   Xgrid.setZboundaryPeriodic(Jx); 
   Xgrid.setZboundaryPeriodic(qex); 
   Xgrid.setZboundaryPeriodic(Pii_xx); 
   Xgrid.setZboundaryPeriodic(Pii_xz); 
   Xgrid.setZboundaryPeriodic(Pii_zz); 
   Xgrid.communicate(Ez);
   Xgrid.communicate(Jz);
   Xgrid.communicate(qez);
   Xgrid.communicate(Ex);
   Xgrid.communicate(Jx);
   Xgrid.communicate(qex);
   Xgrid.communicate(Pii_xx);
   Xgrid.communicate(Pii_xz);
   Xgrid.communicate(Pii_zz);

   if(stage==1) {
      Ezold = Ez;
      Exold = Ex;
      Jzold = Jz;
      Jxold = Jx;
      qexold = qex;
      qezold = qez;
      Pii_xxold = Pii_xx;
      Pii_xzold = Pii_xz;
      Pii_zzold = Pii_zz;
   }

}


void computeFluxes(const domainGrid& Xgrid, const int order)
{ // order not being used right now

   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXce = Xgrid.nXce;
   const int nXg = Xgrid.nXg;
   const int nZg = Xgrid.nZg;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   matrix2D<double> FluxNcc_x, FluxMxcc_x, FluxMzcc_x, FluxEicc_x, FluxEecc_x;
   matrix2D<double> FluxNcc_z, FluxMxcc_z, FluxMzcc_z, FluxEicc_z, FluxEecc_z;
   matrix2D<double> FluxNecc_x, FluxNecc_z; 
   FluxNcc_x.initialize(nXcc,nZcc,0.0);
   FluxMxcc_x.initialize(nXcc,nZcc,0.0);
   FluxMzcc_x.initialize(nXcc,nZcc,0.0);
   FluxEicc_x.initialize(nXcc,nZcc,0.0);
   FluxEecc_x.initialize(nXcc,nZcc,0.0);
   FluxNecc_x.initialize(nXcc,nZcc,0.0);
   //
   FluxNcc_z.initialize(nXcc,nZcc,0.0);
   FluxMxcc_z.initialize(nXcc,nZcc,0.0);
   FluxMzcc_z.initialize(nXcc,nZcc,0.0);
   FluxEicc_z.initialize(nXcc,nZcc,0.0);
   FluxEecc_z.initialize(nXcc,nZcc,0.0);
   FluxNecc_z.initialize(nXcc,nZcc,0.0);


   // compute x-fluxes at cell center
   //
   FluxNcc_x  = hy_cc*Mx;
   FluxMxcc_x = hy_cc*(Mx*Vx + P);
   FluxMzcc_x = hy_cc*Mz*Vx;
   FluxEicc_x = hy_cc*(Ei + Pi)*Vx;
   FluxEecc_x = hy_cc*(Ee + Pe)*Vx;
   FluxNecc_x = Zbar*FluxNcc_x;
   
   // compute z-fluxes at cell center
   //
   FluxNcc_z  = Mz;
   FluxMxcc_z = Mx*Vz;
   FluxMzcc_z = Mz*Vz + P;
   FluxEicc_z = (Ei + Pi)*Vz;
   FluxEecc_z = (Ee + Pe)*Vz;
   FluxNecc_z = Zbar*FluxNcc_z;
   
   // compute advective flux using
   // specified scheme from input file
   //
   //Nsub = 2;
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxNcc_x,Cspeed_x, hy_cc*N,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxMx_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMxcc_x,Cspeed_x,hy_cc*Mx,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxMz_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMzcc_x,Cspeed_x,hy_cc*Mz,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxEi_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxEicc_x,Cspeed_x,hy_cc*Ei,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxEe_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxEecc_x,Cspeed_x,hy_cc*Ee,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxNe_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxNecc_x,Cspeed_x,hy_cc*Ne,"minmod",0,Nsub);
      //
      Xgrid.computeFluxTVD(FluxN_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxNcc_z,Cspeed_z, N,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxMx_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxMxcc_z,Cspeed_z,Mx,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxMz_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxMzcc_z,Cspeed_z,Mz,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxEi_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxEicc_z,Cspeed_z,Ei,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxEe_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxEecc_z,Cspeed_z,Ee,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxNe_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           FluxNecc_z,Cspeed_z,Ne,"minmod",1,Nsub);
   }
   else {
      Xgrid.computeFluxTVDsimple(FluxN_x, FluxL_x,FluxR_x,FluxNcc_x, Cspeed_x,hy_cc*N, advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxMx_x,FluxL_x,FluxR_x,FluxMxcc_x,Cspeed_x,hy_cc*Mx,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxMz_x,FluxL_x,FluxR_x,FluxMzcc_x,Cspeed_x,hy_cc*Mz,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxEi_x,FluxL_x,FluxR_x,FluxEicc_x,Cspeed_x,hy_cc*Ei,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxEe_x,FluxL_x,FluxR_x,FluxEecc_x,Cspeed_x,hy_cc*Ee,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxNe_x,FluxL_x,FluxR_x,FluxNecc_x,Cspeed_x,hy_cc*Ne,advScheme0,0);
      // 
      Xgrid.computeFluxTVDsimple(FluxN_z, FluxL_z,FluxR_z,FluxNcc_z, Cspeed_z,N, advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxMx_z,FluxL_z,FluxR_z,FluxMxcc_z,Cspeed_z,Mx,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxMz_z,FluxL_z,FluxR_z,FluxMzcc_z,Cspeed_z,Mz,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxEi_z,FluxL_z,FluxR_z,FluxEicc_z,Cspeed_z,Ei,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxEe_z,FluxL_z,FluxR_z,FluxEecc_z,Cspeed_z,Ee,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxNe_z,FluxL_z,FluxR_z,FluxNecc_z,Cspeed_z,Ne,advScheme0,1);
   } 

   vector<double> P0;
   P0.assign(nZcc,0.0);
   if(procID==0) {
      Xgrid.setXminBoundary(FluxN_x, 0.0, 0.0);
      //cout << "hy_ce.at(1)" << hy_ce.at(1) << endl;
      for (auto j=0; j<nZcc; j++) {
         P0.at(j) = hy_ce(nXg-1,j)*(P(nXg,j)+P(nXg-1,j))/2.0;
      }   
      Xgrid.setXminBoundary(FluxMx_x, P0);   
      Xgrid.setXminBoundary(FluxMz_x, 0.0, 0.0);   
      Xgrid.setXminBoundary(FluxEi_x, 0.0, 0.0);
      Xgrid.setXminBoundary(FluxEe_x, 0.0, 0.0);
      Xgrid.setXminBoundary(FluxNe_x, 0.0, 0.0);
   }
   if(procID==numProcs-1) {
      Xgrid.setXmaxBoundary(FluxN_x, 0.0, 0.0);
      Xgrid.setXmaxBoundary(FluxEi_x, 0.0, 0.0);
      Xgrid.setXmaxBoundary(FluxEe_x, 0.0, 0.0);
      const int thisnX = P.size0();
      //cout << "hy_ce.at(thisnX-3)" << hy_ce.at(thisnX-3) << endl;   
      //double P0 = hy_ce.at(thisnX-nXg-1)*(P.at(thisnX-nXg-1)+P.at(thisnX-nXg))/2.0;
      for (auto j=0; j<nZcc; j++) {
         P0.at(j) = hy_ce(thisnX-nXg-1,j)*(P(thisnX-nXg-1,j)+P(thisnX-nXg,j))/2.0;
      }   
      Xgrid.setXmaxBoundary(FluxMx_x, P0);   
      Xgrid.setXmaxBoundary(FluxMz_x, 0.0, 0.0);   
      Xgrid.setXmaxBoundary(FluxNe_x, 0.0, 0.0);
   }   
   Xgrid.communicate(FluxN_x);   
   Xgrid.communicate(FluxMx_x);   
   Xgrid.communicate(FluxMz_x);   
   Xgrid.communicate(FluxEi_x);   
   Xgrid.communicate(FluxEe_x);   
   Xgrid.communicate(FluxNe_x);   
   
}

void updatePhysicalVars( const domainGrid&  Xgrid )
{
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   Vx  = Mx/N;
   Vz  = Mz/N;
   Zbar = Ne/N;
   Nn = N-Ne;
   if(min(N)<0.0)    cout << " N IS LESS THAN ZERO " << endl;
   if(min(Zbar)<Zmin*0.99) {
      cout << " Zbar IS LESS THAN Zmin " << endl;
      cout << " min(Zbar) = " << min(Zbar) << endl;
   }
   if(max(Zbar)>1.0)  {
      cout << " Zbar IS LARGER THAN ONE " << endl;
      cout << " max(Zbar) = " << max(Zbar) << endl;
   }

   Pi  = (Ei - 0.5*(Vx*Mx + Vz*Mz))*(gamma0-1.0);
   Pe  = Ee*(gamma0-1.0);
   Pe_eff = Pe - mM*Pi;
   P = Pe + Pi; 
   Ti  = Pi/N;
   Te  = Pe/Ne;
   if(min(Te)<0.0) cout << " Te IS LESS THAN ZERO on proc " << procID << endl;
   if(min(Ti)<0.0) cout << " Ti IS LESS THAN ZERO on proc " << procID << endl;
   Cs = sqrt(gamma0*P/N + By*By/(N+1.0e-4));  // 1.0e-4 here for time step
   Cspeed_x  = abs(Vx) + Cs;                  // adv flux jacobian
   Cspeed_z  = abs(Vz) + Cs;                  // adv flux jacobian
   
   // interp electron density to edges
   //
   Xgrid.InterpToCellEdges(Ne_x,Ne,Ne,"C2",0);
   Xgrid.InterpToCellEdges(Ne_z,Ne,Ne,"C2",1);
   Xgrid.communicate(Ne_x);
   
   // interp density to edges
   //
   Xgrid.InterpToCellEdges(N_x,N,N,"C2",0);
   Xgrid.InterpToCellEdges(N_z,N,N,"C2",1);
   Xgrid.communicate(N_x);

   // compute Hall and electron velocity at cell center
   //
   Vhallx_z = -Li0*Jx/Ne_z;
   Vhallz_x = -Li0*Jz/Ne_x;
   Vhallx = -Li0*Jxcc/Ne;
   Vhallz = -Li0*Jzcc/Ne;
   Vex = Vx + Vhallx;
   Vez = Vz + Vhallz;

}

void updateJccAndEcc( const domainGrid&  Xgrid )
{

   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   // compute curlB at cell edges and at cell center
   //
   computeJ0(Jx0,Jz0,By,Xgrid);
   bool useJ0forJcc = false;
   if(useJ0forJcc) {   
      Xgrid.InterpToCellCenter(Jzcc,Jz0);
      Xgrid.InterpToCellCenter(Jxcc,Jx0);
   } else {
      Xgrid.InterpToCellCenter(Jzcc,Jz);
      Xgrid.InterpToCellCenter(Jxcc,Jx);
   }
   Xgrid.communicate(Jzcc);
   Xgrid.communicate(Jxcc);
   
   
   // compute electric field at cell center
   ///
   Xgrid.InterpToCellCenter(Ezcc,Ez);
   Xgrid.InterpToCellCenter(Excc,Ex);
   Xgrid.communicate(Ezcc);
   
}


void updateCollisionTerms( const domainGrid&  Xgrid )
{
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg  = Xgrid.nXg;
   const int nZg  = Xgrid.nZg;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   matrix2D<double> nue, Te_eV;
   nue.initialize(nXcc,nZcc,0.0);
   Te_eV.initialize(nXcc,nZcc,0.0);
   

   //  set collision times and resistivity
   //
   Te_eV = Te*Tscale; 
   setNeutralInteractionRates(Xgrid,Te_eV,Ne,Nn);

   nue_spi = pow(Tscale,1.5)/taue0*Ne*( 1.0/pow(Te_eV,1.5) );
   nue_vac = pow(Tscale,1.5)/taue0*Ne*( 0.01*pow(NvacC/N,NvacP) );
   nue = nue_spi + nue_neu + nue_vac;
   
   double nueR_min0 = pow(Tscale,1.5)/taue0*(1.0/pow(0.1,1.5));
   vector<double> nueR_min;
   nueR_min.assign(nZcc,nueR_min0);
   const int thisnX = nue.size0();
   if(procID==0 && XlowBC.compare("insulator") == 0) { // set resistivity at insulator boundary
      for (auto j=0; j<nZcc; j++) {
         if(nue(1,j) < nueR_min0) {
            Xgrid.setXminBoundary(nue, nueR_min);
            break;
         }
      }   
   }
   if(procID==numProcs-1 && XhiBC.compare("insulator") == 0) { // set resistivity at insulator boundary
      for (auto j=0; j<nZcc; j++) {
         if(nue(thisnX-nXg+1,j) < nueR_min0) {
            Xgrid.setXmaxBoundary(nue, nueR_min);
            break;
         }
      }   
   }
   
   eta = eta0*taue0*nue/Ne;   
   taue = 1.0/nue; 
   taui = taui0*pow(Ti,1.5)/N; 
   Qie   = 2.0/(gamma0-1.0)*mM*(nue_spi + nue_neu)*Ne*(Te-Ti); // use taue_spi here for time-step

   Xgrid.InterpToCellEdges(eta_x,eta,eta,"C2",0);
   Xgrid.InterpToCellEdges(eta_z,eta,eta,"C2",1);
   Xgrid.communicate(eta_x);


   // set heat flux coefficients
   //

}

void computeJ0( matrix2D<double>&  a_Jx0_z,
                matrix2D<double>&  a_Jz0_x,
          const matrix2D<double>&  a_By_cc,
          const domainGrid&        Xgrid )
{
   const int nXce = Xgrid.nXce;
   const int nZce = Xgrid.nZce;
   assert(a_Jx0_z.size1()==nZce);
   assert(a_Jz0_x.size0()==nXce);
   
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg =  Xgrid.nXg;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
      
   Xgrid.DDZ(a_Jx0_z,a_By_cc);
   a_Jx0_z = -1.0*a_Jx0_z;
   
   Xgrid.DDX(a_Jz0_x,hy_cc*a_By_cc);
   a_Jz0_x = a_Jz0_x/hy_ce;
   
   vector<double> J00;
   J00.assign(nZcc,0.0);
   if(procID==0 && geometry0=="CYL" && XlowBC.compare("axis") == 0) {
      for (auto j=0; j<nZcc; j++) {
         J00.at(j) = 2.0*a_By_cc(nXg,j)/hy_cc(nXg,j);
      }
      Xgrid.setXminBoundary(a_Jz0_x,J00);
   } 
   Xgrid.communicate(a_Jz0_x);
   Xgrid.communicate(a_Jx0_z);

}

void setJ0( const domainGrid&  Xgrid )
{
   // compute each component of J0=curlB on each edge
   //

   const int nXce = Xgrid.nXce;
   const int nZce = Xgrid.nZce;
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg =  Xgrid.nXg;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
    
   // compute J0x on z-edge and J0z on x-edge
   // and set appropriate BC's
   //
   Xgrid.DDZ(J0x_z,By);
   J0x_z *= -1.0;
   
   Xgrid.DDX(J0z_x,hy_cc*By);
   J0z_x = J0z_x/hy_ce;
   
   vector<double> J00;
   J00.assign(nZcc,0.0);
   if(procID==0 && geometry0=="CYL" && XlowBC.compare("axis") == 0) {
      for (auto j=0; j<nZcc; j++) {
         J00.at(j) = 2.0*By(nXg,j)/hy_cc(nXg,j);
      }
      Xgrid.setXminBoundary(J0z_x,J00);
   } 
   Xgrid.setZboundaryPeriodic(J0x_z);
   Xgrid.communicate(J0z_x);

   // interp edge values to center
   //
   Xgrid.InterpToCellCenter(J0x,J0x_z);
   Xgrid.InterpToCellCenter(J0z,J0z_x);
   Xgrid.setZboundaryPeriodic(J0x);
   Xgrid.communicate(J0z);

   // interp center values to opposite edges
   //
   Xgrid.InterpToCellEdges(J0x_x,J0x,J0x,"C2",0);
   Xgrid.InterpToCellEdges(J0z_z,J0z,J0z,"C2",1);
   
}

void computeHeatFluxes( const domainGrid&  Xgrid )
{
   // compute heat fluxes at cell edges
   // qex0 = -kappae*dTe/dx
   // qez0 = -kappae*dTe/dz
   // kappae = kappae0*Pe*taue
   //

   const int nXce = Xgrid.nXce;
   const int nZce = Xgrid.nZce;
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg =  Xgrid.nXg;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
    

   // compute kappae at cell center
   //
   matrix2D<double> xe, xe2, xe4;
   xe.initialize(nXcc,nZcc,0.0);
   xe2.initialize(nXcc,nZcc,0.0);
   xe4.initialize(nXcc,nZcc,0.0);
   
   // compute magnetization coefficient
   //
   xe = wcescale*By*taue*tscale; // Omega_ce*taue 
   xe2 = xe*xe;
   xe4 = xe2*xe2;
   chie_perp = (4.664*xe2 + 11.92)/(xe4 + 14.79*xe2 + 3.770); 

   // compute kappae at cell edges
   //
   kappae = Ve0sq*Pe*taue*chie_perp;
   Xgrid.InterpToCellEdges(kappae_x,kappae,kappae,"C2",0);
   Xgrid.InterpToCellEdges(kappae_z,kappae,kappae,"C2",1);

   // compute gradient of Te
   //
   Xgrid.DDX(qex0,Te);
   Xgrid.DDZ(qez0,Te);
   Xgrid.communicate(qex0);
   Xgrid.setZboundaryPeriodic(qez0);

   // multily kappae by grad(Te) and negate
   //
   qex0 *= kappae_x;   
   qex0 *= -1.0;   
   qez0 *= kappae_z;   
   qez0 *= -1.0;   

}

void computeIonViscStressTensor( const domainGrid&  Xgrid )
{
   // compute ion stress tensor for momentum and energy equations
   //

   const int nXce = Xgrid.nXce;
   const int nZce = Xgrid.nZce;
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg =  Xgrid.nXg;
   
   // set ion viscousity coefficient
   //
   etaVis_ion = 0.96*Pi*min(taui,tauiMax);
   Xgrid.InterpToCellEdges(etaVis_ion_x,etaVis_ion,etaVis_ion,"C2",0);
   Xgrid.InterpToCellEdges(etaVis_ion_z,etaVis_ion,etaVis_ion,"C2",1);
   Xgrid.communicate(etaVis_ion_x);
   
   // compute ion strain tensor
   //
   computeStrainTensor( Pii_xx0, Pii_zz0, Pii_yy, Pii_xz0, Vx, Vz, Xgrid );
  
   // multiply ion strain tensor by -viscosity to get ion stress tensor
   //
   Pii_xx0 *= -etaVis_ion_x;
   Pii_zz0 *= -etaVis_ion_z;
   Pii_yy  *= -etaVis_ion;
   Pii_xz0 *= -etaVis_ion;

}   

void computeEleViscStressTensor( const domainGrid&  Xgrid )
{
   // compute ele stress tensor for Ohms law
   //

   const int nXce = Xgrid.nXce;
   const int nZce = Xgrid.nZce;
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg =  Xgrid.nXg;

   // set electron viscosity coefficient
   //
   etaVis_ele = 0.73*Pe*min(taue,taueMax);
   //etaVis_ele += 0.73*Pe*1.0e-3;
   Xgrid.InterpToCellEdges(etaVis_ele_x,etaVis_ele,etaVis_ele,"C2",0);
   Xgrid.InterpToCellEdges(etaVis_ele_z,etaVis_ele,etaVis_ele,"C2",1);
   Xgrid.communicate(etaVis_ele_x);

   // compute ele strain tensor
   //
   //computeStrainTensorStag( Pie_xx0, Pie_zz0, Pie_yy, Pie_xz0, Vex_z, Vez_x);
   computeStrainTensorStag( Pie_xx0, Pie_zz0, Pie_yy, Pie_xz0, Vhallx_z, Vhallz_x);
  
   // multiply ele strain tensor by -viscosity to get ele stress tensor
   //
   Pie_xx0 *= -etaVis_ele_x;
   Pie_zz0 *= -etaVis_ele_z;
   Pie_yy  *= -etaVis_ele;
   Pie_xz0 *= -etaVis_ele;

}

void computeStrainTensor( matrix2D<double>&  a_Pi_xx,
                          matrix2D<double>&  a_Pi_zz,
                          matrix2D<double>&  a_Pi_yy,
                          matrix2D<double>&  a_Pi_xz,
                    const matrix2D<double>&  a_Vx,
                    const matrix2D<double>&  a_Vz,
                    const domainGrid&        Xgrid )
{
   // compute strain tensor for viscosity
   //
   
   const int nXce = Xgrid.nXce;
   const int nZce = Xgrid.nZce;
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg =  Xgrid.nXg;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
    

   // compute Vx derivatives
   //
   matrix2D<double> dVxdx_x, dVxdz;
   dVxdx_x.initialize(nXce,nZcc,0.0); 
   dVxdz.initialize(nXcc,nZcc,0.0); 
   //
   Xgrid.DDX(dVxdx_x,a_Vx);
   Xgrid.DDZ(dVxdz,a_Vx);


   // compute Vz derivatives
   //
   matrix2D<double> dVzdz_z, dVzdx;
   dVzdz_z.initialize(nXcc,nZce,0.0); 
   dVzdx.initialize(nXcc,nZcc,0.0); 
   //
   Xgrid.DDZ(dVzdz_z,a_Vz);
   Xgrid.DDX(dVzdx,a_Vz);
   Xgrid.communicate(dVzdx);
   

   // compute divergence of velocity
   //
   matrix2D<double> divUz, divUx, divU, divUx_x, divU_x, divU_z;
   divUz.initialize(nXcc,nZcc,0.0); 
   divUx.initialize(nXcc,nZcc,0.0); 
   divU.initialize(nXcc,nZcc,0.0); 
   divUx_x.initialize(nXce,nZcc,0.0); 
   divU_x.initialize(nXce,nZcc,0.0); 
   divU_z.initialize(nXcc,nZce,0.0); 

   // div at cell center
   //
   Xgrid.DDZ(divUz,a_Vz);
   Xgrid.setZboundaryPeriodic(divUz);
   Xgrid.DDX(divUx,a_Vx*hy_cc);
   divUx /= hy_cc;
   Xgrid.communicate(divUx);
   divU = divUx + divUz;

   // div at z-edge
   //
   Xgrid.InterpToCellEdges(divU_z,divUx,divUx,"C2",1);
   divU_z += dVzdz_z;
   
   // div at x-edge
   //
   Xgrid.InterpToCellEdges(divU_x,divUz,divUz,"C2",0);
   Xgrid.DDX(divUx_x,hy_cc*a_Vx);
   divUx_x /= hy_ce;
   vector<double> divUx0;
   divUx0.assign(nZcc,0.0);
   if(procID==0 && geometry0=="CYL" && XlowBC.compare("axis") == 0) {
      for (auto j=0; j<nZcc; j++) {
         divUx0.at(j) = 2.0*a_Vx(nXg,j)/hy_cc(nXg,j);
      }
      Xgrid.setXminBoundary(divUx_x,divUx0);
   } 
   Xgrid.communicate(divUx_x);
   divU_x += divUx_x;
   

   // set components of stress tensor on cell edges
   // for momentum equation
   //
   a_Pi_xx = (2.0*dVxdx_x - 2.0/3.0*divU_x);
   a_Pi_zz = (2.0*dVzdz_z - 2.0/3.0*divU_z);
   a_Pi_yy = (2.0*a_Vx/hy_cc  - 2.0/3.0*divU);
   a_Pi_xz = (dVxdz + dVzdx);

}

void computeStrainTensorStag( matrix2D<double>&  a_Pi_xx,
                              matrix2D<double>&  a_Pi_zz,
                              matrix2D<double>&  a_Pi_yy,
                              matrix2D<double>&  a_Pi_xz,
                        const matrix2D<double>&  a_Vx,
                        const matrix2D<double>&  a_Vz )
{
   // compute strain tensor for viscosity
   // with input velocities defined on cell-edges
   //
   assert(a_Vx.size1()==m_nZce);
   assert(a_Vz.size0()==m_nXce);

   const int nXce = m_Xgrid.nXce;
   const int nZce = m_Xgrid.nZce;
   const int nXcc = m_Xgrid.nXcc;
   const int nZcc = m_Xgrid.nZcc;
   const int nXg =  m_Xgrid.nXg;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
    
   // compute Vx derivatives
   //
   matrix2D<double> dVxdx_x, dVxdz;
   dVxdx_x.initialize(nXce,nZcc,0.0); 
   dVxdz.initialize(nXcc,nZcc,0.0); 
   //
   //m_Xgrid.DDX(dVxdx_x,a_Vx);
   //m_Xgrid.communicate(dVxdx_x);
   //m_Xgrid.DDZ(dVxdz,a_Vx);
   //m_Xgrid.setZboundaryPeriodic(dVxdz);


   // compute Vz derivatives
   //
   matrix2D<double> dVzdz_z, dVzdx;
   dVzdz_z.initialize(nXce,nZce,0.0); 
   dVzdx.initialize(nXcc,nZcc,0.0); 
   //
   m_Xgrid.DDZ(dVzdz_z,a_Vz);
   m_Xgrid.DDX(dVzdx,a_Vz); 
   if(procID==0) m_Xgrid.setXminBoundary(dVzdx,0.0,-1.0);
   if(procID==numProcs-1) m_Xgrid.setXmaxBoundary(dVzdx,0.0,-1.0);
   m_Xgrid.communicate(dVzdx);
   m_Xgrid.setZboundaryPeriodic(dVzdz_z);
   
   /*
   // compute divergence of velocity
   //
   matrix2D<double> divUz, divUx, divU, divUx_x, divU_x, divU_z;
   divUz.initialize(nXcc,nZcc,0.0); 
   divUx.initialize(nXcc,nZcc,0.0); 
   divU.initialize(nXcc,nZcc,0.0); 
   divUx_x.initialize(nXce,nZcc,0.0); 
   divU_x.initialize(nXce,nZcc,0.0); 
   divU_z.initialize(nXcc,nZce,0.0); 

   // div at cell center
   //
   m_Xgrid.DDZ(divUz,a_Vz);
   m_Xgrid.setZboundaryPeriodic(divUz);
   m_Xgrid.DDX(divUx,a_Vx*hy_cc);
   divUx /= hy_cc;
   m_Xgrid.communicate(divUx);
   divU = divUx + divUz;

   // div at z-edge
   //
   m_Xgrid.InterpToCellEdges(divU_z,divUx,divUx,"C2",1);
   divU_z += dVzdz_z;
   
   // div at x-edge
   //
   m_Xgrid.InterpToCellEdges(divU_x,divUz,divUz,"C2",0);
   m_Xgrid.DDX(divUx_x,hy_cc*a_Vx);
   divUx_x /= hy_ce;
   vector<double> divUx0;
   divUx0.assign(nZcc,0.0);
   if(procID==0 && geometry0=="CYL" && XlowBC.compare("axis") == 0) {
      for (auto j=0; j<nZcc; j++) {
         divUx0.at(j) = 2.0*a_Vx(nXg,j)/hy_cc(nXg,j);
      }
      m_Xgrid.setXminBoundary(divUx_x,divUx0);
   } 
   m_Xgrid.communicate(divUx_x);
   divU_x += divUx_x;
   */

   // set components of stress tensor on cell edges
   // for momentum equation
   //
   //a_Pi_xx = (2.0*dVxdx_x - 2.0/3.0*divU_x);
   //a_Pi_zz = (2.0*dVzdz_z - 2.0/3.0*divU_z);
   //a_Pi_yy = (2.0*a_Vx/hy_cc  - 2.0/3.0*divU);
   a_Pi_xz = dVxdz + dVzdx;

}

void computeIdealEatEdges( const domainGrid&      Xgrid, 
                           const matrix2D<double>&  a_Cspeed,
                           const matrix2D<double>&  a_By )
{
   // Eideal = -VxB
   //

   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   matrix2D<double> VxBy, VzBy;
   VxBy.initialize(nXcc,nZcc,0.0);
   VzBy.initialize(nXcc,nZcc,0.0);
   VxBy = Vx*a_By;
   VzBy = Vz*a_By;
   if(advScheme0=="TVD") {
      Xgrid.computeFluxTVD(Eidealz,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           VxBy,abs(Vx),a_By,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(Eidealx,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           VzBy,abs(Vz),a_By,"minmod",1,Nsub);
   }
   else {
      Xgrid.computeFluxTVDsimple(Eidealz,FluxL_x,FluxR_x,VxBy,abs(Vx),a_By,advScheme0,0);
      Xgrid.computeFluxTVDsimple(Eidealx,FluxL_z,FluxR_z,VzBy,abs(Vz),a_By,advScheme0,1);
   }
   if(procID==0) {
      Xgrid.setXminBoundary(Eidealz, 0.0, 0.0);
      Xgrid.setXminBoundary(Eidealx, 0.0, 1.0);
   }
   if(procID==numProcs-1) {
      Xgrid.setXmaxBoundary(Eidealz, 0.0, 0.0);
      Xgrid.setXmaxBoundary(Eidealx, 0.0, 1.0);
   }
   Eidealz *= -1.0;
   Xgrid.communicate(Eidealz);
   Xgrid.communicate(Eidealx);

}

void computeHallEatEdges( const domainGrid&  Xgrid )
{
   // Ehall = Li0/Ne*J0xOmega 
   //
   
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXce = Xgrid.nXce;
   const int nZce = Xgrid.nZce;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);


   // get Hall velocity on opposite edges for upwinding
   //
   matrix2D<double> Vhallx_x, Vhallz_z;
   Vhallx_x.initialize(nXce,nZcc,0.0);
   Vhallz_z.initialize(nXcc,nZce,0.0);
   
   Xgrid.InterpToCellEdges(Vhallx_x,Vhallx,Vhallx,"C2",0);
   Xgrid.InterpToCellEdges(Vhallz_z,Vhallz,Vhallz,"C2",1);
   Xgrid.communicate(Vhallx_x);
   Xgrid.communicate(Vhallz_z);

   // upwind By to edges
   //
   Xgrid.InterpCellToEdges(Ehallx,By,Vhallz_z,advSchemeHall0,1);
   Xgrid.InterpCellToEdges(Ehallz,By,Vhallx_x,advSchemeHall0,0);
   Xgrid.communicate(Ehallx);
   Xgrid.communicate(Ehallz);
   Xgrid.setZboundaryPeriodic(Ehallx);
   Xgrid.setZboundaryPeriodic(Ehallz);

   // put it all together
   //
   Ehallx *= Vhallz_z;
   Ehallz *= -1.0*Vhallx_x;

}

void computeDivOfOhmsLawStressTensor( const domainGrid&  Xgrid )
{
   // dF/dt + div(U*F + F*Ue)
   // Ue = U - Li0*J/Ne*F
   // F = J*Le0sq
   //
   // Note that J lives on cell edges. The divergence of the stress
   // tensor is computed by first calculting at cell center in the
   // normal fashion and thn interpolating to the edges
   //
   
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXce = Xgrid.nXce;
   const int nZce = Xgrid.nZce;
   const int nXg = Xgrid.nXg;
   const int nZg = Xgrid.nZg;
   const double dX = Xgrid.dX;
   const double dZ = Xgrid.dZ; 

   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

  
   // compute fluxes at cell center
   //
   matrix2D<double> eJxFluxX_cc, eJxFluxZ_cc, eJzFluxX_cc, eJzFluxZ_cc;
   eJxFluxX_cc.initialize(nXcc,nZcc,0.0);
   eJxFluxZ_cc.initialize(nXcc,nZcc,0.0);
   eJzFluxX_cc.initialize(nXcc,nZcc,0.0);
   eJzFluxZ_cc.initialize(nXcc,nZcc,0.0);
   
   eJxFluxX_cc = Le0sq*(Vx*Jxcc + Jxcc*Vex)*hy_cc;
   eJxFluxZ_cc = Le0sq*(Vz*Jxcc + Jzcc*Vex);
   eJzFluxX_cc = Le0sq*(Vx*Jzcc + Jxcc*Vez)*hy_cc;
   eJzFluxZ_cc = Le0sq*(Vz*Jzcc + Jzcc*Vez);

   // upwind cell center fluxes to cell faces
   //
   const int Nsub = 2;
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(eJxFlux_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           eJxFluxX_cc,Cspeed_x, hy_cc*Le0sq*Jxcc,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(eJzFlux_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           eJzFluxX_cc,Cspeed_x, hy_cc*Le0sq*Jzcc,"minmod",0,Nsub);
      
      Xgrid.computeFluxTVD(eJxFlux_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           eJxFluxZ_cc,Cspeed_z, Le0sq*Jxcc,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(eJzFlux_z,FluxL_z,FluxR_z,FluxRatio_z,FluxLim_z,
                           eJzFluxZ_cc,Cspeed_z, Le0sq*Jzcc,"minmod",1,Nsub);
   }
   else {
      Xgrid.computeFluxTVDsimple(eJxFlux_x,FluxL_x,FluxR_x,eJxFluxX_cc,
                                 Cspeed_x,hy_cc*Le0sq*Jxcc,advScheme0,0);
      Xgrid.computeFluxTVDsimple(eJzFlux_x,FluxL_x,FluxR_x,eJzFluxX_cc,
                                 Cspeed_x,hy_cc*Le0sq*Jxcc,advScheme0,0);
        
      Xgrid.computeFluxTVDsimple(eJxFlux_z,FluxL_z,FluxR_z,eJxFluxZ_cc, 
                                 Cspeed_z,Le0sq*Jxcc,advScheme0,1);
      Xgrid.computeFluxTVDsimple(eJzFlux_z,FluxL_z,FluxR_z,eJzFluxZ_cc,
                                 Cspeed_z,Le0sq*Jzcc,advScheme0,1);
   }
   
   // set flux BCs
   //
   if(procID==0) {
      Xgrid.setXminBoundary(eJxFlux_x, 0.0, 0.0);
      Xgrid.setXminBoundary(eJzFlux_x, 0.0, 0.0);
   }
   if(procID==numProcs-1) {
      //Xgrid.setXmaxBoundary(eJxFlux_x, 0.0, 0.0);
      //Xgrid.setXmaxBoundary(eJzFlux_x, 0.0, 0.0);
   }
 
   // compute divergence of fluxes
   //
   matrix2D<double> divJxstress_cc, divJzstress_cc;
   divJxstress_cc.initialize(nXcc,nZcc,0.0);
   divJzstress_cc.initialize(nXcc,nZcc,0.0);
   for (auto i=nXg; i<nXcc-nXg; i++) {
      for (auto j=nZg; j<nZcc-nZg; j++) {
         divJxstress_cc(i,j) = (eJxFlux_x(i,j) - eJxFlux_x(i-1,j))/hy_cc(i,j)/dX
                             + (eJxFlux_z(i,j) - eJxFlux_z(i,j-1))/dZ;
         divJzstress_cc(i,j) = (eJzFlux_x(i,j) - eJzFlux_x(i-1,j))/hy_cc(i,j)/dX
                             + (eJzFlux_z(i,j) - eJzFlux_z(i,j-1))/dZ;
      }
   }
   
   // apply BC's before interpretting to cell edges
   //
   if(procID==0) {
      Xgrid.setXminBoundary(divJxstress_cc, 0.0, 1.0);
      Xgrid.setXminBoundary(divJzstress_cc, 0.0, 1.0);
   }
   if(procID==numProcs-1) {
      Xgrid.setXmaxBoundary(divJxstress_cc, 0.0, 1.0);
      Xgrid.setXmaxBoundary(divJzstress_cc, 0.0, 1.0);
   }

   // interp divJstress to edges
   //
   Xgrid.InterpCellToEdges(divJxstress,divJxstress_cc,divJxstress,"C2",1);
   Xgrid.InterpCellToEdges(divJzstress,divJzstress_cc,divJzstress,"C2",0);
   Xgrid.communicate(divJxstress);
   Xgrid.communicate(divJzstress);
   Xgrid.setZboundaryPeriodic(divJxstress);
   Xgrid.setZboundaryPeriodic(divJzstress);
   
}

void computeGradEatEdges( const domainGrid&  Xgrid )
{
   // Egrad  = -Li0/Ne*grad(Pe_eff)
   // Pe_eff = Pe - me/Mi*Pi
   //  
 
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   matrix2D<double> gradPex, gradPez;
   gradPez.initialize(nXcc,nZcc,0.0);
   gradPex.initialize(nXcc,nZcc,0.0);
   
   Xgrid.DDX(gradPex,Pe_eff);
   Xgrid.DDZ(gradPez,Pe_eff);
   Xgrid.communicate(gradPex);   
   
   Xgrid.InterpToCellEdges(Egradz,gradPez,gradPez,"C2",0);
   Xgrid.InterpToCellEdges(Egradx,gradPex,gradPex,"C2",1);
   Xgrid.communicate(Egradz);
   
   Egradx *= -Li0/Ne_z;
   Egradz *= -Li0/Ne_x;

}

void computeJouleHeatingTerms( const domainGrid&  Xgrid )
{
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXce = Xgrid.nXce;
   const int nZce = Xgrid.nZce;
   const int nXg  = Xgrid.nXg;
   const int nZg  = Xgrid.nZg;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   matrix2D<double> dPedx, dPedz;
   dPedx.initialize(nXcc,nZcc,0.0);
   dPedz.initialize(nXcc,nZcc,0.0);


   // set momentum density source terms... move this
   //
   JcrossBx = -Jzcc*By;
   JcrossBz =  Jxcc*By;

  
   //  calculate work source terms for energy equations
   //
   if(useIonScaleCorrections) {
      NeUdotE = Ne*(Vz*(Ezcc - eta*Jzcc) + Vx*(Excc - eta*Jxcc));
      NeUdotE /= Li0;
   } 
   else {
      Xgrid.DDX(dPedx,Pe);
      Xgrid.DDZ(dPedz,Pe);
      Xgrid.communicate(dPedx);   
      NeUdotE = Vz*(Jxcc*By - dPedz) - Vx*(Jzcc*By + dPedx);
   }
   JdotE = Jzcc*Ezcc + Jxcc*Excc;
   

   ////////////
   //
   // compute divergence of heat flux
   //
   matrix2D<double> dummydiv;
   dummydiv.initialize(nXcc,nZcc,0.0); 
   
   Xgrid.DDX(dummydiv,hy_ce*qex);
   dummydiv /= hy_cc;
   Xgrid.communicate(dummydiv);
   Xgrid.DDZ(divqe,qez);
   divqe += dummydiv;
   
   if(tauiMax>0.0) {
     
      // interp Pii_xz to cell edges in x and z
      //
      matrix2D<double> Pii_xz_x, Pii_xz_z;
      Pii_xz_x.initialize(nXce,nZcc,0.0);
      Pii_xz_z.initialize(nXcc,nZce,0.0);
      Xgrid.InterpToCellEdges(Pii_xz_x,Pii_xz,Pii_xz,"C2",0);
      Xgrid.InterpToCellEdges(Pii_xz_z,Pii_xz,Pii_xz,"C2",1);

      ////////////
      //
      // compute divergence of stress tensor for x-momentum equation
      //
      Xgrid.DDX(divPii_x,Pii_xx*hy_ce);
      divPii_x /= hy_cc;
      Xgrid.DDZ(dummydiv,Pii_xz_z);
      divPii_x += dummydiv;
      if(geometry0=="CYL") divPii_x -= Pii_yy/hy_cc;
      Xgrid.communicate(divPii_x);
   
 
      ////////////
      //
      // compute divergence of stress tensor for z-momentum equation
      //
      Xgrid.DDX(divPii_z,Pii_xz_x*hy_ce);
      divPii_z /= hy_cc;
      Xgrid.DDZ(dummydiv,Pii_zz);
      divPii_z += dummydiv;
      Xgrid.communicate(divPii_z);
   

      ////////////
      //
      // interpolate velocity to edges and 
      // compute stress tensor source term
      // in energy equation
      //
      matrix2D<double> Vx_x, Vx_z, Vz_x, Vz_z;
      Vx_x.initialize(nXce,nZcc,0.0); 
      Vz_x.initialize(nXce,nZcc,0.0); 
      Vx_z.initialize(nXcc,nZce,0.0); 
      Vz_z.initialize(nXcc,nZce,0.0); 
      Xgrid.InterpToCellEdges(Vx_x,Vx,Vx,"C2",0);
      Xgrid.InterpToCellEdges(Vz_x,Vz,Vz,"C2",0);
      Xgrid.InterpToCellEdges(Vx_z,Vx,Vx,"C2",1);
      Xgrid.InterpToCellEdges(Vz_z,Vz,Vz,"C2",1);
   
   
      // set components of stress tensor for energy equations
      //
      UidotPii_x = Vx_x*Pii_xx   + Vz_x*Pii_xz_x; 
      UidotPii_z = Vx_z*Pii_xz_z + Vz_z*Pii_zz; 

   
      // take divergence for energy equation
      //
      Xgrid.DDX(divUidotPii,UidotPii_x*hy_ce);
      divUidotPii /= hy_cc;
      Xgrid.DDZ(dummydiv,UidotPii_z);
      divUidotPii += dummydiv;
      Xgrid.communicate(divUidotPii);

   }

} 

     
void computeDivOfElectronViscosity( const domainGrid&  Xgrid )
{
      
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   // interp Pie_xz to cell edges in x and z
   //
   matrix2D<double> dummy_onx, dummy_onz;
   dummy_onx.initialize(m_nXce,m_nZcc,0.0);
   dummy_onz.initialize(m_nXcc,m_nZce,0.0);

   // compute divergence of viscous stress tensor for Jx equation
   //
   m_Xgrid.DDX(divPiex,Pie_xx*hy_ce);
   divPiex /= hy_cc;
   if(geometry0=="CYL") divPiex -= Pie_yy/hy_cc;
   m_Xgrid.InterpToCellEdges(divPiex_z,divPiex,divPiex_z,"C2",1);
   m_Xgrid.DDZ(dummy_onz,Pie_xz);
   divPiex_z += dummy_onz;
   m_Xgrid.communicate(divPiex_z);
   m_Xgrid.setZboundaryPeriodic(divPiex_z); 
  
      
   // compute divergence of viscous stress tensor for Jz equation
   //
   m_Xgrid.DDX(divPiez_x,Pie_xz*hy_cc);
   divPiez_x /= hy_ce;
   vector<double> divPiexz0;
   divPiexz0.assign(m_nZcc,0.0);
   if(procID==0 && geometry0=="CYL" && XlowBC.compare("axis") == 0) {
      for (auto j=0; j<m_nZcc; j++) {
         divPiexz0.at(j) = 2.0*Pie_xz(m_nXg,j)/hy_cc(m_nXg,j);
      }
      m_Xgrid.setXminBoundary(divPiez_x,divPiexz0);
   } 
   m_Xgrid.communicate(divPiez_x);
    
   matrix2D<double> dummydiv;
   dummydiv.initialize(m_nXcc,m_nZcc,0.0); 
   m_Xgrid.DDZ(dummydiv,Pie_zz);
   m_Xgrid.setZboundaryPeriodic(dummydiv); 
   m_Xgrid.InterpToCellEdges(dummy_onx,dummydiv,dummy_onx,"C2",0);
   divPiez_x += dummy_onx;

   // 
   //divPiez_x *= Li0;
   //divPiex_z *= Li0;

}

void setNeutralInteractionRates( const domainGrid&     Xgrid, 
                                 const matrix2D<double>&  a_Te_eV, 
                                 const matrix2D<double>&  a_Ne,
                                 const matrix2D<double>&  a_Nn )
{
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   matrix2D<double> kizn, K3br, kmom, TeoU;
   kizn.initialize(nXcc,nZcc,0.0);
   TeoU.initialize(nXcc,nZcc,0.0);
   K3br.initialize(nXcc,nZcc,0.0);
   kmom.initialize(nXcc,nZcc,1.0e-7); // [cm^3/s]

   //if(!procID) cout << "Te_eV = " << a_Te_eV.at(10) << endl;
   //cout << "Ne    = " << a_Ne*Nscale/1.0e6 << endl;
   //cout << "Nn    = " << a_Nn*Nscale/1.0e6 << endl;

   // compute ionization and 3-body recombination rate constant
   //
   TeoU = (1.0/Uizn)*a_Te_eV;
   kizn = 1.0e-5*sqrt(TeoU)/Uizn/sqrt(Uizn)/(6.0+TeoU)*exp(-1.0/TeoU); // [cm^3/s]
   K3br = kizn*1.66e-22*g1/g0/pow(a_Te_eV,1.5)*exp(1.0/TeoU); // [cm^6/s]
 
   // set charge density source/sink rate
   //
   nue_izn = Nscale/1.0e6*a_Nn;
   nue_izn *= kizn; // - pow(a_Ne*Nscale/1.0e6,2)*K3br; // [1/s]
   nue_izn *= tscale; // normalized

   // compute momentum-exchange rate constant and frequency
   //
   //kmom = 1.0e-7; // [cm^3/s]
   nue_neu = Nscale/1.0e6*kmom*a_Nn; // [1/s] 
   nue_neu *= tscale;           // normalized
   nue_neu += nue_izn;         // add growth term

   // compute normalized source terms for charge density and electron energy density
   //
   Se = nue_izn*a_Ne;          // electron density source term
   SEe = -1.0*Uizn/Tscale*Se; // Electron energy density source term

}

void initializeMemberVars( const domainGrid& Xgrid )
{ 
   const int nXce = Xgrid.Xce.size();
   const int nXcc = Xgrid.Xcc.size();
   const int nZce = Xgrid.Zce.size();
   const int nZcc = Xgrid.Zcc.size();
   const int nXg = Xgrid.nXg;
   const int nZg = Xgrid.nZg;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   deltaN.initialize(nXcc,nZcc,0.0);
   
   N.initialize( nXcc,nZcc,0.0);
   Mx.initialize( nXcc,nZcc,0.0);
   Mz.initialize( nXcc,nZcc,0.0);
   Ei.initialize(nXcc,nZcc,0.0);
   Ee.initialize(nXcc,nZcc,0.0);
   By.initialize( nXcc,nZcc,0.0);
   //
   Nold.initialize( nXcc,nZcc,0.0);
   Mxold.initialize( nXcc,nZcc,0.0);
   Mzold.initialize( nXcc,nZcc,0.0);
   Eiold.initialize(nXcc,nZcc,0.0);
   Eeold.initialize(nXcc,nZcc,0.0);
   Byold.initialize( nXcc,nZcc,0.0);
   //
   P.initialize(  nXcc,nZcc,0.0);
   Pi.initialize( nXcc,nZcc,0.0);
   Pe.initialize( nXcc,nZcc,0.0);
   Ti.initialize( nXcc,nZcc,0.0);
   Te.initialize( nXcc,nZcc,0.0);
   eta.initialize(nXcc,nZcc,0.0);
   Cs.initialize( nXcc,nZcc,0.0);
   Vx.initialize( nXcc,nZcc,0.0);
   Vz.initialize( nXcc,nZcc,0.0);
   Jzcc.initialize(  nXcc,nZcc,0.0);
   Jxcc.initialize(  nXcc,nZcc,0.0);
   Ezcc.initialize(  nXcc,nZcc,0.0);
   Excc.initialize(  nXcc,nZcc,0.0);
   Qvisc.initialize( nXcc,nZcc,0.0);
   Cspeed_x.initialize(nXcc,nZcc,0.0);
   Cspeed_z.initialize(nXcc,nZcc,0.0);
   //
   Qie.initialize(nXcc,nZcc,0.0);
   NeUdotE.initialize(nXcc,nZcc,0.0);
   JdotE.initialize(nXcc,nZcc,0.0);
   JcrossBx.initialize(nXcc,nZcc,0.0);
   JcrossBz.initialize(nXcc,nZcc,0.0);
   taue.initialize(nXcc,nZcc,0.0);
   taui.initialize(nXcc,nZcc,0.0);
   nue_spi.initialize(nXcc,nZcc,0.0);
   nue_vac.initialize(nXcc,nZcc,0.0);
   
   // Ez and Jz are defined on x-edges
   //
   Ez.initialize(nXce,nZcc,0.0);
   Ezold.initialize(nXce,nZcc,0.0);
   Jz.initialize(nXce,nZcc,0.0);
   Jzold.initialize(nXce,nZcc,0.0);
   Jz0.initialize(nXce,nZcc,0.0);
   Jz0stagTime.initialize(nXce,nZcc,0.0);
   eta_x.initialize(nXce,nZcc,0.0);
   Ne_x.initialize(nXce,nZcc,0.0);
   N_x.initialize(nXce,nZcc,0.0);
   Eidealz.initialize(nXce,nZcc,0.0);
   
   // Ex and Jx are defined on z-edges
   //
   Ex.initialize(nXcc,nZce,0.0);
   Exold.initialize(nXcc,nZce,0.0);
   Jx.initialize(nXcc,nZce,0.0);
   Jxold.initialize(nXcc,nZce,0.0);
   Jx0.initialize(nXcc,nZce,0.0);
   Jx0stagTime.initialize(nXcc,nZce,0.0);
   eta_z.initialize(nXcc,nZce,0.0);
   Ne_z.initialize(nXcc,nZce,0.0);
   N_z.initialize(nXcc,nZce,0.0);
   Eidealx.initialize(nXcc,nZce,0.0);
 

 
   // curlB at cell center and both edges
   //
   J0x.initialize(nXcc,nZcc,0.0);
   J0z.initialize(nXcc,nZcc,0.0);
   J0x_x.initialize(nXce,nZcc,0.0);
   J0z_x.initialize(nXce,nZcc,0.0);
   J0x_z.initialize(nXcc,nZce,0.0);
   J0z_z.initialize(nXcc,nZce,0.0);
    


   //
   FluxRatio_x.initialize(nXce,nZcc,0.0);
   FluxLim_x.initialize(nXce,nZcc,0.0);
   FluxN_x.initialize(nXce,nZcc,0.0);
   FluxMx_x.initialize(nXce,nZcc,0.0);
   FluxMz_x.initialize(nXce,nZcc,0.0);
   FluxEe_x.initialize(nXce,nZcc,0.0);
   FluxEi_x.initialize(nXce,nZcc,0.0);
   FluxNe_x.initialize(nXce,nZcc,0.0);
   FluxR_x.initialize(nXce,nZcc,0.0);
   FluxL_x.initialize(nXce,nZcc,0.0);
   //
   hy_cc.initialize(nXcc,nZcc,1.0);
   hy_ce.initialize(nXce,nZcc,1.0);
   
   
   FluxRatio_z.initialize(nXcc,nZce,0.0);
   FluxLim_z.initialize(nXcc,nZce,0.0);
   FluxN_z.initialize(nXcc,nZce,0.0);
   FluxMx_z.initialize(nXcc,nZce,0.0);
   FluxMz_z.initialize(nXcc,nZce,0.0);
   FluxEe_z.initialize(nXcc,nZce,0.0);
   FluxEi_z.initialize(nXcc,nZce,0.0);
   FluxNe_z.initialize(nXcc,nZce,0.0);
   FluxR_z.initialize(nXcc,nZce,0.0);
   FluxL_z.initialize(nXcc,nZce,0.0);
   //

   // ionization terms
   //
   //Zbar.initialize(nXcc,nZcc,Zmin);
   Neold.initialize(nXcc,nZcc,0.0);
   Se.initialize(nXcc,nZcc,0.0);
   SEe.initialize(nXcc,nZcc,0.0);
   Ne.initialize(nXcc,nZcc,0.0);
   Nn.initialize(nXcc,nZcc,0.0);
   nue_neu.initialize(nXcc,nZcc,0.0);
   nue_izn.initialize(nXcc,nZcc,0.0);

   // finite ion and electron inertial length terms
   //
   Ehallx.initialize(nXcc,nZce,0.0); 
   Ehallz.initialize(nXce,nZcc,0.0); 
   Egradx.initialize(nXcc,nZce,0.0); 
   Egradz.initialize(nXce,nZcc,0.0); 
   Vhallx.initialize(nXcc,nZcc,0.0); 
   Vhallz.initialize(nXcc,nZcc,0.0); 
   Vhallx_z.initialize(nXcc,nZce,0.0);
   Vhallz_x.initialize(nXce,nZcc,0.0);
   Vex.initialize(nXcc,nZcc,0.0); 
   Vez.initialize(nXcc,nZcc,0.0); 
   divJxstress.initialize(nXcc,nZce,0.0); 
   divJzstress.initialize(nXce,nZcc,0.0); 
   Pe_eff.initialize(nXcc,nZcc,0.0); 
   curlJ0.initialize(nXcc,nZcc,0.0); 
   curlVe.initialize(nXcc,nZcc,0.0); 
   eJxFlux_x.initialize(nXce,nZcc,0.0); 
   eJzFlux_x.initialize(nXce,nZcc,0.0); 
   eJxFlux_z.initialize(nXcc,nZce,0.0); 
   eJzFlux_z.initialize(nXcc,nZce,0.0); 

   // electron heat flux
   //
   qex.initialize(nXce,nZcc,0.0);
   qexold.initialize(nXce,nZcc,0.0);
   qex0.initialize(nXce,nZcc,0.0);
   kappae_x.initialize(nXce,nZcc,0.0);
   qez.initialize(nXcc,nZce,0.0);
   qezold.initialize(nXcc,nZce,0.0);
   qez0.initialize(nXcc,nZce,0.0);
   kappae_z.initialize(nXcc,nZce,0.0);
   //
   divqe.initialize(nXcc,nZcc,0.0);
   kappae.initialize(nXcc,nZcc,0.0);
   chie_perp.initialize(nXcc,nZcc,0.0);

   
   // ion viscosity
   //
   Pii_xxold.initialize(nXce,nZcc,0.0);
   Pii_xx0.initialize(nXce,nZcc,0.0);
   Pii_xx.initialize(nXce,nZcc,0.0);
   Pii_xzold.initialize(nXcc,nZcc,0.0);
   Pii_xz0.initialize(nXcc,nZcc,0.0);
   Pii_xz.initialize(nXcc,nZcc,0.0);
   Pii_zzold.initialize(nXcc,nZce,0.0);
   Pii_zz0.initialize(nXcc,nZce,0.0);
   Pii_zz.initialize(nXcc,nZce,0.0);
   Pii_yy.initialize(nXcc,nZcc,0.0);

   divPii_x.initialize(nXcc,nZcc,0.0);
   divPii_z.initialize(nXcc,nZcc,0.0);
   divUidotPii.initialize(nXcc,nZcc,0.0);
   etaVis_ion.initialize(nXcc,nZcc,0.0);
   etaVis_ion_x.initialize(nXce,nZcc,0.0); 
   etaVis_ion_z.initialize(nXcc,nZce,0.0); 
   UidotPii_x.initialize(nXce,nZcc,0.0);
   UidotPii_z.initialize(nXcc,nZce,0.0);
   

   // ele viscosity
   //
   Pie_xxold.initialize(nXce,nZcc,0.0);
   Pie_xx0.initialize(nXce,nZcc,0.0);
   Pie_xx.initialize(nXce,nZcc,0.0);
   Pie_xzold.initialize(nXcc,nZcc,0.0);
   Pie_xz0.initialize(nXcc,nZcc,0.0);
   Pie_xz.initialize(nXcc,nZcc,0.0);
   Pie_zzold.initialize(nXcc,nZce,0.0);
   Pie_zz0.initialize(nXcc,nZce,0.0);
   Pie_zz.initialize(nXcc,nZce,0.0);
   Pie_yy.initialize(nXcc,nZcc,0.0);

   divPiex.initialize(nXcc,nZcc,0.0);
   divPiez.initialize(nXcc,nZcc,0.0);
   divUedotPie.initialize(nXcc,nZcc,0.0);
   etaVis_ele.initialize(nXcc,nZcc,0.0);
   etaVis_ele_x.initialize(nXce,nZcc,0.0); 
   etaVis_ele_z.initialize(nXcc,nZce,0.0); 
   UedotPie_x.initialize(nXce,nZcc,0.0);
   UedotPie_z.initialize(nXcc,nZce,0.0);
   divPiex_z.initialize(nXcc,nZce,0.0);
   divPiez_x.initialize(nXce,nZcc,0.0);


}

void addMembersToDataFile( HDF5dataFile&  dataFile )
{
   dataFile.add(N, "N", 1);      // density 
   dataFile.add(Mx, "Mx", 1);      // momentum density 
   dataFile.add(Mz, "Mz", 1);      // momentum density 
   dataFile.add(By, "By", 1);      // magnetic field
   dataFile.add(Ei, "Ei", 1);    // total ion energy
   dataFile.add(Ee, "Ee", 1);    // total ele energy
   dataFile.add(P, "P", 1);      // total pressure
   dataFile.add(Pi, "Pi", 1);    // ion pressure
   dataFile.add(Pe, "Pe", 1);    // ele pressure
   dataFile.add(Ti, "Ti", 1);    // ion temperature
   dataFile.add(Te, "Te", 1);    // ele temperature
   dataFile.add(eta, "eta", 1);  // resistivity
   dataFile.add(taue, "taue", 1);  // ele collision time
   dataFile.add(taui, "taui", 1);  // ion collision time
   dataFile.add(nue_spi, "nue_spi", 1);  // spitzer coll freq
   dataFile.add(nue_vac, "nue_vac", 1);  // vac coll freq
   
   dataFile.add(Ne,  "Ne", 1);       // Ne = Zbar*N   [cm^3/s]
   dataFile.add(Zbar,"Zbar", 1);     // Ne = Zbar*N   [cm^3/s]
   dataFile.add(nue_neu, "nue_neu", 1);   // e-n coll freq [Hz]
   dataFile.add(nue_izn, "nue_izn", 1);  // e-n izn - recom [Hz]
   
   // electron heat flux
   //
   dataFile.add(qex,"qex", 1); 
   dataFile.add(qez,"qez", 1);  
   dataFile.add(qex0,"qex0", 1); 
   dataFile.add(qez0,"qez0", 1);  
   dataFile.add(chie_perp,"chie_perp", 1);  
   dataFile.add(kappae,"kappae", 1);  
   dataFile.add(divqe,"divqe", 1);  
   
   // ion viscosity
   //
   dataFile.add(Pii_xx0,"Pii_xx0", 1); 
   dataFile.add(Pii_xx,"Pii_xx", 1); 
   dataFile.add(Pii_xz0,"Pii_xz0", 1); 
   dataFile.add(Pii_xz,"Pii_xz", 1); 
   dataFile.add(Pii_zz0,"Pii_zz0", 1); 
   dataFile.add(Pii_zz,"Pii_zz", 1); 
   dataFile.add(etaVis_ion,"etaVis_ion", 1); 
   dataFile.add(divPii_x,"divPii_x", 1); 
   dataFile.add(divPii_z,"divPii_z", 1); 
   dataFile.add(UidotPii_x,"UidotPii_z", 1); 
   dataFile.add(divUidotPii,"divUidotPii", 1); 
   
   // ele viscosity
   //
   dataFile.add(Pie_xx0,"Pie_xx0", 1); 
   dataFile.add(Pie_xx,"Pie_xx", 1); 
   dataFile.add(Pie_xz0,"Pie_xz0", 1); 
   dataFile.add(Pie_xz,"Pie_xz", 1); 
   dataFile.add(Pie_zz0,"Pie_zz0", 1); 
   dataFile.add(Pie_zz,"Pie_zz", 1); 
   dataFile.add(etaVis_ele,"etaVis_ele", 1); 
   dataFile.add(divPiex_z,"divPiex_z", 1); 
   dataFile.add(divPiez_x,"divPiez_x", 1); 
   dataFile.add(UedotPie_x,"UidotPie_z", 1); 
   dataFile.add(divUedotPie,"divUidotPie", 1); 


   //dataFile.add(Vx, "Vx", 1);      // velocity
   //dataFile.add(Vz, "Vz", 1);      // velocity
   dataFile.add(Jz, "Jz", 1);     // z-current density
   dataFile.add(Jx, "Jx", 1);     // x-current density
   dataFile.add(Jzcc, "Jzcc", 1); // z-current density at cell-center
   dataFile.add(Jxcc, "Jxcc", 1); // x-current density at cell-center
   dataFile.add(Jz0, "Jz0", 1);   // curl of B
   dataFile.add(Jx0, "Jx0", 1);   // curl of B
   dataFile.add(Ez, "Ez", 1);          // z-electric field
   dataFile.add(Ex, "Ex", 1);          // x-electric field
   dataFile.add(Ezcc, "Ezcc", 1);          // z-electric field
   dataFile.add(Excc, "Excc", 1);          // x-electric field
   dataFile.add(Ehallz, "Ehallz", 1);  // z-Hall electric field
   dataFile.add(Ehallx, "Ehallx", 1);  // x-Hall electric field
   dataFile.add(Egradz, "Egradz", 1);  // z-grad electric field
   dataFile.add(Egradx, "Egradx", 1);  // x-grad electric field
   dataFile.add(curlJ0, "curlJ0", 1);  // curlJ0 = -delp2(By)
   dataFile.add(curlVe, "curlVe", 1);  // curlVe 
   dataFile.add(Vhallx, "Vhallx", 1);  // Vhallx = - Li0/Ne*Jx
   dataFile.add(Vhallz, "Vhallz", 1);  // Vhallz = - Li0/Ne*Jz
   dataFile.add(Vex, "Vex", 1);    // Vex = Vx - Li0/Ne*Jx
   dataFile.add(Vez, "Vez", 1);    // Vez = Vz - Li0/Ne*Jz
   dataFile.add(eJxFlux_x, "eJxFlux_x", 1);
   dataFile.add(eJxFlux_z, "eJxFlux_z", 1);
   dataFile.add(eJzFlux_x, "eJzFlux_x", 1);
   dataFile.add(eJzFlux_z, "eJzFlux_z", 1);
   dataFile.add(divJxstress, "divJxstress", 1);  // x-comp of Ohm's law stress tensor
   dataFile.add(divJzstress, "divJzstress", 1);  // z-comp of Ohm's law stress tensor
   dataFile.add(Cs,"Cs",1);          // sound speed
  
   dataFile.add(JdotE,"JdotE",1);
   dataFile.add(NeUdotE,"NeUdotE",1);
   dataFile.add(Qie,"Qie",1);
   dataFile.add(SEe,"SEe",1);
  
   //dataFile.add(Cspeed_x,"Cspeed_x",1);  // max char speed
   dataFile.add(gamma0,"gamma0",0); 
   dataFile.add(delta0,"delta0",0); 
   dataFile.add(Li0,"Li0",0); 
   dataFile.add(Le0sq,"Le0sq",0); 
   //
   dataFile.add(FluxN_x, "FluxN_x", 1);  
   dataFile.add(FluxMx_x, "FluxMx_x", 1);  
   dataFile.add(FluxMz_x, "FluxMz_x", 1);  
   dataFile.add(FluxEi_x, "FluxEi_x", 1);  
   dataFile.add(FluxEe_x, "FluxEe_x", 1);  
   dataFile.add(FluxNe_x, "FluxNe_x", 1);  
   //
   //dataFile.add(FluxRatio_x, "FluxRatio_x", 1);  
   //dataFile.add(FluxLim_x, "FluxLim_x", 1);  
   //dataFile.add(FluxR_x, "FluxR_x", 1);
   //dataFile.add(FluxL_x, "FluxL_x", 1);
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
   dataFile.add(wcescale,"wcescale",0);
   dataFile.add(Mi,"Mi",0);
   dataFile.add(Mn,"Mn",0);
   //
   dataFile.add(hy_cc,"hy_cc",0);
   dataFile.add(hy_ce,"hy_ce",0);

}

void parseInputFile( const domainGrid& Xgrid, const Json::Value& a_root )
{ 
   const int nXce = Xgrid.nXce;
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg = Xgrid.nXg;
   const int nZg = Xgrid.nZg;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = a_root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      Json::Value advScheme = Phys.get("advScheme",defValue);
      Json::Value advSchemeHall = Phys.get("advSchemeHall",defValue);
      Json::Value useEleScale = Phys.get("useEleScaleCorrections",defValue);
      Json::Value useIonScale = Phys.get("useIonScaleCorrections",defValue);
      Json::Value geometry  = Phys.get("geometry",defValue);
      Json::Value gammaVal  = Phys.get("gammaC",defValue);
      Json::Value ZminVal   = Phys.get("Zmin",defValue);
      Json::Value NsubVal   = Phys.get("Nsub",defValue);
      Json::Value tauiMaxVal = Phys.get("tauiMax",defValue);
      Json::Value taueMaxVal = Phys.get("taueMax",defValue);
     
  
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
      Mn      = Amass*amu;                  // neutral mass [kg]
      Mi      = Mn - me;                    // ion mass [kg]
      Vscale  = pow(Pscale/Mn/Nscale,0.5);  // velocity scale [m/s]
      Bscale  = pow(mu0*Pscale,0.5);        // magnetic field scale [T]
      Jscale  = Bscale/Xscale/mu0;          // current density scale [A/m^2]
      Tscale  = Pscale/Nscale/qe;           // temperature scale [eV]
      Ezscale = Vscale*Bscale;              // electric field scale [V/m]
      tscale  = Xscale/Vscale;              // time scale [s]
      double etascale  = Xscale*Xscale*mu0/tscale; // resistivity scale [Ohm-m]
      mM = me/Mn;
      double wpescale = 5.64e4*pow(Nscale/1.0e6,0.5); // ele plasma freq [rad/s]
      double wpiscale = wpescale*pow(me/Mn,0.5);      // ion plasma freq [rad/s]
      wcescale = qe*Bscale/me;   // ele cyclotron freq [rad/s]
      double wciscale = qe*Bscale/Mn;   // ion cyclotron freq [rad/s]
      double tauescale = 3.44e5/10.0*pow(Tscale,1.5)/(Nscale/1.0e6); // collision time [s]
      double tauiscale = 2.09e7/10.0*pow(Tscale,1.5)/(Nscale/1.0e6)*sqrt(Mi/Mp); // collision time [s]
      //
      double Lescale  = cvac/wpescale;    // ele inertial scale [m]
      double Liscale  = cvac/wpiscale;    // ion inertial scale with neutral mass [m]
      
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
         eta0   = 1.03e-4/10.0/pow(Tscale,1.5)/(Xscale*Xscale*mu0/tscale); // norm res
      } else {
	 eta0 = etaVal.asDouble();
      } 
      // 
      taue0 = tauescale/tscale;
      taui0 = tauiscale/tscale;
      double Vte  = 4.19e5*sqrt(Tscale); // characteristic ele therm speed [m/s]
      Ve0sq = Vte*Vte/Vscale/Vscale;
      delta0 = pow(Vscale/cvac,2.0)*epsilonRel;
      Le0sq = pow(Lescale/Xscale,2.0)*meRel;
      Li0   = sqrt(Le0sq/mM/meRel); // no meRel here
      B00   = sqrt(2.0);  // see normalization
      if(procID==0) {
	 cout << endl;
         cout << "dimensionless parameters:" << endl;
         cout << "normalized resistivity = " << eta0 << endl;
	 if(etaVal != defValue) cout << "WARNING: USING eta0 FROM INPUT FILE !!!" << endl;
         cout << "taue/tscale = " << taue0 << endl;
         cout << "taui/tscale = " << taui0 << endl;
         cout << "wce*taue = " << wcescale*tauescale << endl;
         cout << "wci*taui = " << wciscale*tauiscale << endl;
         cout << "(Vte/V0)^2 = " << Ve0sq  << endl;
         cout << "(Li0/r0)   = " << Li0   << " norm ion inert length " << endl;      
         cout << "(Le0/r0)^2 = " << Le0sq << " norm ele inert length squared" << endl;      
         cout << "(V0/c)^2   = " << delta0   << " (Jz relaxation const)" << endl;      
      }

      if(advScheme == defValue || gammaVal == defValue ||
	 geometry == defValue || NsubVal == defValue || 
         ZminVal == defValue || advSchemeHall == defValue ) {
         cout << "ERROR: advScheme/Hall or gamma " << endl;
         cout << "or Nsub or geometry or Zmin" << endl;
         cout << "not declared in input file" << endl;
         exit (EXIT_FAILURE);
      } 
      
      advScheme0 = advScheme.asString();
      if(procID==0) {
         cout << "advection diff/interp scheme is " << advScheme0 << endl;
      }
      advSchemeHall0 = advSchemeHall.asString();
      if(procID==0) {
         cout << "Hall advection diff/interp scheme is " << advSchemeHall0 << endl;
      }
      useEleScaleCorrections = useEleScale.asBool();
      if(procID==0 && useEleScaleCorrections) {
         cout << "electron scale corrections being used" << endl;
      }
      useIonScaleCorrections = useIonScale.asBool();
      if(procID==0 && useIonScaleCorrections) {
         cout << "Ion scale corrections being used" << endl;
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

      Zmin = ZminVal.asDouble();
      if(procID==0) cout << "Zmin = " << Zmin << endl;
      //Zbar;      
      Zbar.initialize(nXcc,nZcc,Zmin);

      Nsub = NsubVal.asInt();
      if(procID==0) cout << "Nsub = " << Nsub << endl;
      if(gamma0 < 1.0) {
         printf("ERROR: Nsub must be int >= 1\n");
         exit (EXIT_FAILURE);
      }
      Pthresh = Nthresh*Tthresh/Tscale;
      Ethresh = Pthresh/(gamma0-1.0);

      tauiMax = tauiMaxVal.asDouble();
      if(procID==0) cout << "max taui for visc = " << tauiMax << endl;
      
      taueMax = taueMaxVal.asDouble();
      if(procID==0) cout << "max taue for visc = " << taueMax << endl;

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
   
   // initialize plasma density perturbation
   //
   const Json::Value deltaNvar = Phys.get("deltaN",defValue);
   if(deltaNvar.isObject()) { 
      Xgrid.setInitialProfile(deltaN,deltaNvar);
      if(procID==0) Xgrid.setXminBoundary(deltaN, 0.0, 1.0);   
      if(procID==numProcs-1) Xgrid.setXmaxBoundary(deltaN, 0.0, 1.0);   
      Xgrid.setZboundaryPeriodic(deltaN);
      Xgrid.communicate(deltaN);
   } else {
      cout << "value for Physics variable \"deltaN\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   N = N*(1.0+deltaN);

   // initialize plasma velocity
   //
   const Json::Value Vxvar = Phys.get("Vx",defValue);
   if(Vxvar.isObject()) { 
      Xgrid.setInitialProfile(Vx,Vxvar);
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
   const Json::Value Byvar = Phys.get("By",defValue);
   if(Byvar.isObject()) { 
      Xgrid.setInitialProfile(By,Byvar);
      //B0 = 0;
      //if(procID==0) Xgrid.setXminBoundary(By, 0.0, 0.0);   
      //if(procID==numProcs-1) Xgrid.setXmaxBoundary_J(By,B0,0);
      //Xgrid.communicate(By);
   } else {
      cout << "value for Physics variable \"By\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }

}

void setXminExtrap(matrix2D<double>& var)
{
   /*
   //var.at(1) = 2.0*var.at(2) - var.at(3);
   //var.at(0) = 2.0*var.at(1) - var.at(2);
   var.at(1) = 3.0*(var.at(2) - var.at(3)) + var.at(4);
   var.at(0) = 3.0*(var.at(1) - var.at(2)) + var.at(3);
   //var.at(0) = var.at(1);
   */
}


void Physics::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid, const int verbose)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   //matrix2D<double> Cchar;
   //Cchar.initialize(Cs.size0(),Cs.size1(),0.0);
   double Cmax_x, Cmax_z, nue_izn_max;
   double Vhall_max_x, Vhall_max_z;
   //Cchar = abs(Vx)+Cs;
   //Cmax = max(Cchar);
   Cmax_x = max(Cspeed_x);
   Cmax_z = max(Cspeed_z);
   Vhall_max_x = max(abs(Vhallx));
   Vhall_max_z = max(abs(Vhallz));
   nue_izn_max = max(abs(nue_izn));
   //cout << "Cmax_x = " << Cmax_x << endl;
   //cout << "abs(Vx) = " << max(abs(Vx))<< endl;
   //cout << "Cs = " << max(Cs)<< endl;

   const double dX = Xgrid.dX;
   const double dZ = Xgrid.dZ;

   double dtCFL_sound_x = dX/Cmax_x;
   double dtCFL_sound_z = dZ/Cmax_z;
   double dtCFL_sound = min(dtCFL_sound_x,dtCFL_sound_z);
   
   double dtCFL_hall_x = dX/Vhall_max_x;
   double dtCFL_hall_z = dZ/Vhall_max_z;
   double dtCFL_hall = min(dtCFL_hall_x,dtCFL_hall_z);

   double dt_izn = 1.0/nue_izn_max;

   double dtCFL_light = dX*sqrt(delta0);

   double dtmax = min(dtCFL_sound,dtCFL_light);
   dtmax = min(dtCFL_hall,dtmax);
   dtmax = min(dt_izn,dtmax);
   
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);
   //dtSim = 5.0e-5;
   if(procID==0 && verbose) {
      cout << "sigma_0*dt/delta = " << dtSim/delta0/eta0 << endl;
      cout << "dtCFL_sound = " << dtCFL_sound << endl;
      cout << "dtCFL_hall  = " << dtCFL_hall << endl;
      cout << "dt_izn      = " << dt_izn << endl;
      cout << "dtCFL_light = " << dtCFL_light << endl;
      cout << "dtSim = " << dtSim << endl;
   }
}

