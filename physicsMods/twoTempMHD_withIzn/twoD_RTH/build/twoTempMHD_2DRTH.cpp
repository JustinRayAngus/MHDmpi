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
double m_dY, m_nYce, m_nYcc;
int m_nXg, m_nZg, m_nYg;

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
double tauiMax=0.0;     // maximum collision times for viscosity
double taueMin=0.0, taueMax=0.0;  // min/max collision times for electron viscosity
double delta0;      // relaxation const (V0/cvac)^2
double Li0;         // normalized ion skin depth
double Le0sq;       // normalized electron skin depth squared
double B0, B00;     // boundary value of magnetic field
int Nsub;           // time-solver subcycle steps
matrix2D<double> N, Mx, My, Ei, Ee;   // time-evolving variables
matrix2D<double> Nold, Mxold, Myold, Eiold, Eeold;
matrix2D<double> eta, Cs, Cspeed_x, Cspeed_y, Vx, Vy, P, Pe, Pi, Te, Ti;
matrix2D<double> eta_x, Ne_x, N_x, eta_y, Ne_y, N_y, Qvisc;
//
matrix2D<double> qex, qexold, qex0, kappae_x;  // ele heat flux
matrix2D<double> qey, qeyold, qey0, kappae_y;
matrix2D<double> chie_perp, kappae, divqe;
//
matrix2D<double> qix, qixold, qix0, kappai_x;  // ion heat flux
matrix2D<double> qiy, qiyold, qiy0, kappai_y;
matrix2D<double> chii_perp, kappai, divqi;
//
matrix2D<double> Bx, Bxold, Bxcc;
matrix2D<double> By, Byold, Bycc;
matrix2D<double> Bz, Bzold;
matrix2D<double> Ex, Exold, Excc, Eidealx;
matrix2D<double> Ey, Eyold, Eycc, Eidealy;
matrix2D<double> Ez, Ezold, Ezcc, Eidealz;
matrix2D<double> Jx, Jxold, Jxcc, Jx0, Jx0stagTime;
matrix2D<double> Jy, Jyold, Jycc, Jy0, Jy0stagTime;
matrix2D<double> Jz, Jzold, Jzcc, Jz0, Jz0stagTime;
matrix2D<double> eta_xy, Ne_xy;
//
//
//
matrix2D<double> FluxRatio_x, FluxLim_x, FluxR_x, FluxL_x;
matrix2D<double> FluxRatio_y, FluxLim_y, FluxR_y, FluxL_y;  
matrix2D<double> FluxRatio_xy, FluxLim_xy, FluxR_xy, FluxL_xy, Vx_y, Vy_x;
matrix2D<double> FluxN_x, FluxMx_x, FluxMy_x, FluxEi_x, FluxEe_x;
matrix2D<double> FluxN_y, FluxMx_y, FluxMy_y, FluxEi_y, FluxEe_y;
matrix2D<double> Qie, Qiwall, NeUdotE, JdotE, JcrossBx, JcrossBy, taue, taui, nue_spi, nue_vac;
matrix2D<double> deltaN;

// terms for finite ion and electron inertial length corrections
//
matrix2D<double> Vex, Vey, Vez, Vhallx, Vhally, Vhallz, Vhallx_y, Vhally_x;
matrix2D<double> Ehallx, Ehallz, Egradx, Egradz; 
matrix2D<double> Pe_eff;
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
matrix2D<double> divPiex_z, divPiez_x;
matrix2D<double> etaVis_ele, etaVis_ele_x, etaVis_ele_z, etaVis_ele_xz;


// new stuff for ionization
//
const double Uizn = 13.6;
const double g0 = 2, g1=1;
double Zmin = 1.0;
matrix2D<double> Zbar, Neold, Se, SEe, Ne, Nn, nue_neu, nue_izn;
matrix2D<double> FluxNe_x, FluxNe_y;

double Nscale, Xscale, Amass, Iscale, dtIscale;
double Pscale, Tscale, Bscale, Jscale, Ezscale, Vscale, tscale, Mi, Mn;
double wcescale, epsilonRel, meRel;

matrix2D<double> hy_cc, hy_x, hy_y, hy_xy; // y-dimension lame coefficient

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
void computeJ0( matrix2D<double>&, matrix2D<double>&, matrix2D<double>&, 
          const matrix2D<double>&, 
          const matrix2D<double>&, 
          const matrix2D<double>& );
void updatePhysicalVars(const domainGrid& );
void updateJccAndEcc( const domainGrid& );
void updateCollisionTerms(const domainGrid& );
void computeDivOfOhmsLawStressTensor( const domainGrid& );
void computeDivOfElectronViscosity( const domainGrid& );
void computeIdealEatEdges();
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
   m_nXce  = Xgrid.nXce2;
   m_nZcc  = Xgrid.nZcc;
   m_nZce  = Xgrid.nZce2;
   m_dX    = Xgrid.dX;
   m_dZ    = Xgrid.dZ;
   //
   m_nYg   = Xgrid.nZg;
   m_dY    = Xgrid.dZ;
   m_nYcc  = Xgrid.nZcc;
   m_nYce  = Xgrid.nZce2;


   //domainGrid* mesh = domainGrid::mesh;
   const int nXg  = Xgrid.nXg;
   const int nXcc = Xgrid.Xcc.size();
   
   initializeMemberVars(Xgrid);

   parseInputFile(Xgrid,root);

   // define remaining vars from those specified in input file
   //
   Ne = Zbar*N;
   Nn = N-Ne; 
   Mx  = N*Vx;
   My  = N*Vy;
   Pe = Ne*Te;
   Pi = N*Ti;
   Ei = 0.5*(Mx*Mx + My*My)/N + Pi/(gamma0-1.0);
   Ee = Pe/(gamma0-1.0);

   // initialize current density from magnetic field
   //
   computeJ0(Jx0,Jy0,Jz0,Bx,By,Bz);
   Jx = Jx0;
   Jy = Jy0;
   Jz = Jz0;
   Jx0stagTime = Jx;
   Jy0stagTime = Jy;
   Jz0stagTime = Jz;
   
   Xgrid.DDX(Jzcc,hy_cc*Bycc);
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
   qey = qey0;
   qix = qix0;
   qiy = qiy0;

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

   computeIdealEatEdges();
   
   if(useIonScaleCorrections) {
      computeHallEatEdges(Xgrid);
      computeGradEatEdges(Xgrid);
      computeDivOfOhmsLawStressTensor(Xgrid);
   }

   //Ez = eta_x*Jz + Eidealz; // initialize electric field from Ohm's law 
   //Ex = eta_y*Jx + Eidealx; // initialize electric field from Ohm's law 

   computeJouleHeatingTerms(Xgrid);
   
   computeFluxes(Xgrid, 1); // 1 for first order calculation   

   // dont forget to initialize all the old values
   //
   Nold  = N;
   Neold = Ne;
   Mxold = Mx;  
   Myold = My;  
   Eiold = Ei;
   Eeold = Ee;
   //
   Bxold = Bx;
   Byold = By;
   Bzold = Bz;
   //
   Exold = Ex;
   Eyold = Ey;
   Ezold = Ez;
   Jxold = Jx;
   Jyold = Jy;
   Jzold = Jz;
   //
   qexold = qex;
   qeyold = qey;
   qixold = qix;
   qeyold = qey;
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
      
      computeIdealEatEdges();
      
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
      
   for (auto i=m_nXg; i<m_nXce-m_nXg; i++) {
      for (auto j=m_nYg; j<m_nYcc-m_nYg; j++) {

	 N(i,j)  = Nold(i,j)  - thisdt*(FluxN_x(i+1,j)  - FluxN_x(i,j) )/hy_cc(i,j)/m_dX
	                      - thisdt*(FluxN_y(i,j+1)  - FluxN_y(i,j) )/hy_cc(i,j)/m_dY;
         Mx(i,j) = Mxold(i,j) - thisdt*(FluxMx_x(i+1,j) - FluxMx_x(i,j))/hy_cc(i,j)/m_dX
	                      - thisdt*(FluxMx_y(i,j+1) - FluxMx_y(i,j))/hy_cc(i,j)/m_dY
	                      - thisdt*divPii_x(i,j)
                              + thisdt*JcrossBx(i,j);
	 if(geometry0=="CYL") Mx(i,j) = Mx(i,j) + thisdt*(P(i,j)+Vy(i,j)*My(i,j))/hy_cc(i,j);
         My(i,j) = Myold(i,j) - thisdt*(FluxMy_x(i+1,j) - FluxMy_x(i,j))/hy_cc(i,j)/m_dX
                              - thisdt*(FluxMy_y(i,j+1) - FluxMy_y(i,j))/hy_cc(i,j)/m_dY
	                      - thisdt*divPii_z(i,j)
                              + thisdt*JcrossBy(i,j);
	 if(geometry0=="CYL") My(i,j) = My(i,j) - thisdt*(Vy(i,j)*Mx(i,j))/hy_cc(i,j);
         Ei(i,j) = Eiold(i,j) - thisdt*(FluxEi_x(i+1,j) - FluxEi_x(i,j))/hy_cc(i,j)/m_dX
                              - thisdt*(FluxEi_y(i,j+1) - FluxEi_y(i,j))/hy_cc(i,j)/m_dY
	                      - thisdt*(divqi(i,j) + divUidotPii(i,j))
      		              + thisdt*(NeUdotE(i,j) + Qie(i,j) - Qiwall(i,j));
         Ee(i,j) = Eeold(i,j) - thisdt*(FluxEe_x(i+1,j) - FluxEe_x(i,j))/hy_cc(i,j)/m_dX
                              - thisdt*(FluxEe_y(i,j+1) - FluxEe_y(i,j))/hy_cc(i,j)/m_dY
      		              - thisdt*divqe(i,j)
      		              + thisdt*(JdotE(i,j) - Qie(i,j) - NeUdotE(i,j) + SEe(i,j));
	 By(i,j) = Byold(i,j) + thisdt*(Ez(i+1,j) - Ez(i,j))/m_dX;
	 
         Bx(i,j) = Bxold(i,j) - thisdt*(Ez(i,j+1) - Ez(i,j))/hy_x(i,j)/m_dY;
         
         Bz(i,j) = Bzold(i,j) - thisdt*(hy_x(i+1,j)*Ey(i+1,j) - hy_x(i,j)*Ey(i,j))/hy_cc(i,j)/m_dX
                              + thisdt*(Ex(i,j+1) - Ex(i,j))/hy_cc(i,j)/m_dY;
	 
         Ne(i,j) = Neold(i,j) - thisdt*(FluxNe_x(i+1,j) - FluxNe_x(i,j))/hy_cc(i,j)/m_dX
                              - thisdt*(FluxNe_y(i,j+1) - FluxNe_y(i,j))/hy_cc(i,j)/m_dY
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
                 + 0.5*(Mx(i,j)*Mx(i,j)+My(i,j)*My(i,j))/N(i,j);
	 if(Ei(i,j)<=Ethresh) Ei(i,j) = Ethresh;
	 //Ei(i,j) = Ethresh;
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
      m_Xgrid.setXminBoundary(My, 0.0, 0.0);   
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
      m_Xgrid.setXminBoundary(Bz, 0.0, 1.0);

   }
   
   if(procID==numProcs-1) {
      m_Xgrid.setXmaxBoundary(N, 0.0, 1.0);   
      m_Xgrid.setXmaxBoundary(Ne, 0.0, 1.0);   
      m_Xgrid.setXmaxBoundary(Mx, 0.0, -1.0);   
      m_Xgrid.setXmaxBoundary(My, 0.0, 0.0);   
      m_Xgrid.setXmaxBoundary(Ei, 0.0, 1.0);   
      m_Xgrid.setXmaxBoundary(Ee, 0.0, 1.0); 
         
      if(XhiBC.compare("conductor")==0 ) {
         m_Xgrid.setXmaxBoundary_J(By,0,1);   
      } 
      else if(XhiBC.compare("insulator") == 0) {
         //cout << "hy_cc.at(nXcc-nXg) = " << hy_cc.at(nXcc-nXg) << endl;	 
         m_Xgrid.setXmaxBoundary_J(By,B0,0); 
      } 
      else {
         cout << "Xhi MagFieldBC not defined !!!! " << endl;
      }
      m_Xgrid.setXmaxFluxBC(Bx, 0.0, 0.0);
      m_Xgrid.setXmaxBoundary(Bz, 0.0, 1.0);
   }      
      
   m_Xgrid.setZboundaryPeriodic(N); 
   m_Xgrid.setZboundaryPeriodic(Ne); 
   m_Xgrid.setZboundaryPeriodic(Mx); 
   m_Xgrid.setZboundaryPeriodic(My); 
   m_Xgrid.setZboundaryPeriodic(Ei); 
   m_Xgrid.setZboundaryPeriodic(Ee); 
   m_Xgrid.setZboundaryPeriodic(Bx); 
   m_Xgrid.setZboundaryPeriodic(By); 
   m_Xgrid.setZboundaryPeriodic(Bz); 
   m_Xgrid.communicate(N);
   m_Xgrid.communicate(Ne);
   m_Xgrid.communicate(Mx);
   m_Xgrid.communicate(My);
   m_Xgrid.communicate(Ei);
   m_Xgrid.communicate(Ee);
   m_Xgrid.communicate(Bx);
   m_Xgrid.communicate(By);
   m_Xgrid.communicate(Bz);
   
   if(stage==2) {
      // update old values for explicit vars
      //
      Nold = N;
      Neold = Ne;
      Mxold = Mx;
      Myold = My;
      Eiold = Ei;
      Eeold = Ee;
      Bxold = Bx;
      Byold = By;
      Bzold = Bz;
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
   for (auto j=m_nZg; j<m_nZce-m_nZg; j++) {
      for (auto i=m_nXg; i<m_nXce-m_nXg; i++) {
         d0Pie = thisdt/delta0*Ne_xy(i,j)/etaVis_ele_xz(i,j)*mM;
         Pie_xx(i,j) = (Pie_xxold(i,j) + d0Pie*Pie_xx0(i,j))/(1.0 + d0Pie);
         Pie_zz(i,j) = (Pie_zzold(i,j) + d0Pie*Pie_zz0(i,j))/(1.0 + d0Pie);
      }
   }
   for (auto j=m_nZg; j<m_nZcc-m_nZg; j++) {
      for (auto i=m_nXg; i<m_nXcc-m_nXg; i++) {
         d0Pie = thisdt/delta0*Ne(i,j)/etaVis_ele(i,j)*mM;
         Pie_xz(i,j) = (Pie_xzold(i,j) + d0Pie*Pie_xz0(i,j))/(1.0 + d0Pie);
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
      computeJ0(Jx0stagTime,Jy0stagTime,Jz0stagTime,Bx,By,Bz);
   }

   // advance electric field and current density
   //
   double d0, e0_xy, e0_y, e0_x, E0x, E0y, E0z;
   double d0qex, d0qey, d0Pii;
   double d0qix, d0qiy;
   d0 = delta0/thisdt;	

   // update vars stag in x and y
   //
   for (auto j=m_nYg; j<m_nYce-m_nYg; j++) {
      for (auto i=m_nXg; i<m_nXce-m_nXg; i++) {
         E0z  = Eidealz(i,j); // + Ehallz(i,j) + Egradz(i,j);
         //E0z += divJzstress(i,j)/Ne_xy(i,j);
         //E0z -= Li0*divPiez_x(i,j)/Ne_xy(i,j);
         e0_xy = Le0sq/thisdt/Ne_xy(i,j);	   
         //
         Ez(i,j)  = Jz0stagTime(i,j) + d0*Ezold(i,j)
                  - (e0_xy*Jzold(i,j) - E0z)/(e0_xy + eta_xy(i,j));
         Ez(i,j) /= d0 + 1.0/(e0_xy + eta_xy(i,j)); 
         //
         Jz(i,j)  = e0_xy*Jzold(i,j) + Ez(i,j) - E0z;
         Jz(i,j) /= e0_xy + eta_xy(i,j);
         
         /*
         d0Pii = thisdt/delta0*N_x(i,j)/etaVis_ion_x(i,j);
         Pii_xx(i,j) = (Pii_xxold(i,j) + d0Pii*Pii_xx0(i,j))/(1.0 + d0Pii);
         //
         if(i<m_nXcc-m_nXg && j<m_nYcc-m_nYg) {
            d0Pii = thisdt/delta0*N(i,j)/etaVis_ion(i,j);
            Pii_xz(i,j) = (Pii_xzold(i,j) + d0Pii*Pii_xz0(i,j))/(1.0 + d0Pii);
         }
         */
      }
   }
 
   // update vars stag in y
   //
   for (auto j=m_nYg; j<m_nYce-m_nYg; j++) {
      for (auto i=m_nXg; i<m_nXcc-m_nXg; i++) {
         E0x  = Eidealx(i,j); // + Ehallx(i,j) + Egradx(i,j);
         e0_y = Le0sq/thisdt/Ne_y(i,j);	   
         //
         Ex(i,j)  = Jx0stagTime(i,j) + d0*Exold(i,j)
                  - (e0_y*Jxold(i,j) - E0x)/(e0_y + eta_y(i,j));
         Ex(i,j) /= d0 + 1.0/(e0_y + eta_y(i,j)); 
         //
         Jx(i,j)  = e0_y*Jxold(i,j) + Ex(i,j) - E0x;
         Jx(i,j) /= e0_y + eta_y(i,j);
         //
         d0qey = thisdt/delta0*Ne_y(i,j)/kappae_y(i,j);
         qey(i,j)  = (qeyold(i,j) + d0qey*qey0(i,j))/(1.0 + d0qey);
         //
         d0qiy = thisdt/delta0*N_y(i,j)/kappai_y(i,j);
         qiy(i,j) = (qiyold(i,j) + d0qiy*qiy0(i,j))/(1.0 + d0qiy);
         //
         //d0Pii = thisdt/delta0*N_y(i,j)/etaVis_ion_z(i,j);
         //Pii_zz(i,j) = (Pii_zzold(i,j) + d0Pii*Pii_zz0(i,j))/(1.0 + d0Pii);
      }
   }
   
   // update vars stag in x
   //
   for (auto j=m_nYg; j<m_nYcc-m_nYg; j++) {
      for (auto i=m_nXg; i<m_nXce-m_nXg; i++) {
         E0y  = Eidealy(i,j); // + Ehally(i,j) + Egrady(i,j);
         e0_x = Le0sq/thisdt/Ne_x(i,j);	   
         //
         Ey(i,j)  = Jy0stagTime(i,j) + d0*Eyold(i,j)
                  - (e0_x*Jyold(i,j) - E0y)/(e0_x + eta_x(i,j));
         Ey(i,j) /= d0 + 1.0/(e0_x + eta_x(i,j)); 
         //
         Jy(i,j)  = e0_x*Jyold(i,j) + Ey(i,j) - E0y;
         Jy(i,j) /= e0_x + eta_x(i,j);
         //
         d0qex = thisdt/delta0*Ne_x(i,j)/kappae_x(i,j);
         qex(i,j)  = (qexold(i,j) + d0qex*qex0(i,j))/(1.0 + d0qex);
         //
         d0qix = thisdt/delta0*N_x(i,j)/kappai_x(i,j);
         qix(i,j)  = (qixold(i,j) + d0qix*qix0(i,j))/(1.0 + d0qix);
         //
         //d0qix = thisdt/delta0*N_x(i,j)/kappa_x(i,j);
         //qix(i,j)  = (qixold(i,j) + d0qix*qix0(i,j))/(1.0 + d0qix);
      }
   }
   
   if(procID==0) {
      Xgrid.setXminBoundary(Pii_xz, 0.0, -1.0);   
      Xgrid.setXminBoundary(Jz, 0.0, 1.0);   
      Xgrid.setXminBoundary(Jx, 0.0, -1.0);   
   }
   if(procID==numProcs-1) {
      Xgrid.setXmaxBoundary(Pii_xz, 0.0, -1.0);   
      Xgrid.setXmaxBoundaryExtrap(Jz);   
      Xgrid.setXmaxBoundaryExtrap(Jx);   
   }  
   Xgrid.setZboundaryPeriodic(Ey); 
   Xgrid.setZboundaryPeriodic(Jy); 
   Xgrid.setZboundaryPeriodic(Ez); 
   Xgrid.setZboundaryPeriodic(Jz); 
   Xgrid.setZboundaryPeriodic(Ex); 
   Xgrid.setZboundaryPeriodic(Jx); 
   Xgrid.setZboundaryPeriodic(qex); 
   Xgrid.setZboundaryPeriodic(qey); 
   Xgrid.setZboundaryPeriodic(qix); 
   //Xgrid.setZboundaryPeriodic(Pii_xx); 
   //Xgrid.setZboundaryPeriodic(Pii_xz); 
   //Xgrid.setZboundaryPeriodic(Pii_zz); 
   Xgrid.communicate(Ey);
   Xgrid.communicate(Jy);
   Xgrid.communicate(Ez);
   Xgrid.communicate(Jz);
   Xgrid.communicate(Ex);
   Xgrid.communicate(Jx);
   Xgrid.communicate(qex);
   Xgrid.communicate(qey);
   Xgrid.communicate(qix);
   //Xgrid.communicate(Pii_xx);
   //Xgrid.communicate(Pii_xz);
   //Xgrid.communicate(Pii_zz);

   if(stage==1) {
      Exold = Ex;
      Eyold = Ey;
      Ezold = Ez;
      Jxold = Jx;
      Jyold = Jy;
      Jzold = Jz;
      qexold = qex;
      qeyold = qey;
      qixold = qix;
      //Pii_xxold = Pii_xx;
      //Pii_xzold = Pii_xz;
      //Pii_zzold = Pii_zz;
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
   
   matrix2D<double> FluxNcc_x, FluxMxcc_x, FluxMycc_x, FluxEicc_x, FluxEecc_x;
   matrix2D<double> FluxNcc_y, FluxMxcc_y, FluxMycc_y, FluxEicc_y, FluxEecc_y;
   matrix2D<double> FluxNecc_x, FluxNecc_y; 
   FluxNcc_x.initialize( nXcc,nZcc,0.0);
   FluxMxcc_x.initialize(nXcc,nZcc,0.0);
   FluxMycc_x.initialize(nXcc,nZcc,0.0);
   FluxEicc_x.initialize(nXcc,nZcc,0.0);
   FluxEecc_x.initialize(nXcc,nZcc,0.0);
   FluxNecc_x.initialize(nXcc,nZcc,0.0);
   //
   FluxNcc_y.initialize( nXcc,nZcc,0.0);
   FluxMxcc_y.initialize(nXcc,nZcc,0.0);
   FluxMycc_y.initialize(nXcc,nZcc,0.0);
   FluxEicc_y.initialize(nXcc,nZcc,0.0);
   FluxEecc_y.initialize(nXcc,nZcc,0.0);
   FluxNecc_y.initialize(nXcc,nZcc,0.0);


   // compute x-fluxes at cell center
   //
   FluxNcc_x  = hy_cc*Mx;
   FluxMxcc_x = hy_cc*(Mx*Vx + P);
   FluxMycc_x = hy_cc*My*Vx;
   FluxEicc_x = hy_cc*(Ei + Pi)*Vx;
   FluxEecc_x = hy_cc*(Ee + Pe)*Vx;
   FluxNecc_x = Zbar*FluxNcc_x;
   
   // compute z-fluxes at cell center
   //
   FluxNcc_y  = My;
   FluxMxcc_y = Mx*Vy;
   FluxMycc_y = My*Vy + P;
   FluxEicc_y = (Ei + Pi)*Vy;
   FluxEecc_y = (Ee + Pe)*Vy;
   FluxNecc_y = Zbar*FluxNcc_y;
   
   // compute advective flux using
   // specified scheme from input file
   //
   //Nsub = 2;
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxNcc_x,Cspeed_x, hy_cc*N,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxMx_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMxcc_x,Cspeed_x,hy_cc*Mx,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxMy_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxMycc_x,Cspeed_x,hy_cc*My,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxEi_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxEicc_x,Cspeed_x,hy_cc*Ei,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxEe_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxEecc_x,Cspeed_x,hy_cc*Ee,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(FluxNe_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           FluxNecc_x,Cspeed_x,hy_cc*Ne,"minmod",0,Nsub);
      //
      Xgrid.computeFluxTVD(FluxN_y,FluxL_y, FluxR_y,FluxRatio_y,FluxLim_y,
                           FluxNcc_y,Cspeed_y, N,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxMx_y,FluxL_y,FluxR_y,FluxRatio_y,FluxLim_y,
                           FluxMxcc_y,Cspeed_y,Mx,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxMy_y,FluxL_y,FluxR_y,FluxRatio_y,FluxLim_y,
                           FluxMycc_y,Cspeed_y,My,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxEi_y,FluxL_y,FluxR_y,FluxRatio_y,FluxLim_y,
                           FluxEicc_y,Cspeed_y,Ei,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxEe_y,FluxL_y,FluxR_y,FluxRatio_y,FluxLim_y,
                           FluxEecc_y,Cspeed_y,Ee,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(FluxNe_y,FluxL_y,FluxR_y,FluxRatio_y,FluxLim_y,
                           FluxNecc_y,Cspeed_y,Ne,"minmod",1,Nsub);
   }
   else {
      Xgrid.computeFluxTVDsimple(FluxN_x, FluxL_x,FluxR_x,FluxNcc_x, Cspeed_x,hy_cc*N, advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxMx_x,FluxL_x,FluxR_x,FluxMxcc_x,Cspeed_x,hy_cc*Mx,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxMy_x,FluxL_x,FluxR_x,FluxMycc_x,Cspeed_x,hy_cc*My,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxEi_x,FluxL_x,FluxR_x,FluxEicc_x,Cspeed_x,hy_cc*Ei,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxEe_x,FluxL_x,FluxR_x,FluxEecc_x,Cspeed_x,hy_cc*Ee,advScheme0,0);
      Xgrid.computeFluxTVDsimple(FluxNe_x,FluxL_x,FluxR_x,FluxNecc_x,Cspeed_x,hy_cc*Ne,advScheme0,0);
      // 
      Xgrid.computeFluxTVDsimple(FluxN_y, FluxL_y,FluxR_y,FluxNcc_y, Cspeed_y,N, advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxMx_y,FluxL_y,FluxR_y,FluxMxcc_y,Cspeed_y,Mx,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxMy_y,FluxL_y,FluxR_y,FluxMycc_y,Cspeed_y,My,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxEi_y,FluxL_y,FluxR_y,FluxEicc_y,Cspeed_y,Ei,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxEe_y,FluxL_y,FluxR_y,FluxEecc_y,Cspeed_y,Ee,advScheme0,1);
      Xgrid.computeFluxTVDsimple(FluxNe_y,FluxL_y,FluxR_y,FluxNecc_y,Cspeed_y,Ne,advScheme0,1);
   } 

   vector<double> P0;
   P0.assign(nZcc,0.0);
   if(procID==0) {
      Xgrid.setXminFluxBC(FluxN_x, 0.0, 0.0);
      for (auto j=0; j<nZcc; j++) {
         P0.at(j) = hy_x(nXg,j)*(P(nXg,j)+P(nXg-1,j))/2.0;
      }   
      Xgrid.setXminFluxBC(FluxMx_x, P0);   
      Xgrid.setXminFluxBC(FluxMy_x, 0.0, 0.0);   
      Xgrid.setXminFluxBC(FluxEi_x, 0.0, 0.0);
      Xgrid.setXminFluxBC(FluxEe_x, 0.0, 0.0);
      Xgrid.setXminFluxBC(FluxNe_x, 0.0, 0.0);
   }
   if(procID==numProcs-1) {
      Xgrid.setXmaxFluxBC(FluxN_x, 0.0, 0.0);
      Xgrid.setXmaxFluxBC(FluxEi_x, 0.0, 0.0);
      Xgrid.setXmaxFluxBC(FluxEe_x, 0.0, 0.0);
      const int thisnX = P.size0();
      for (auto j=0; j<nZcc; j++) {
         P0.at(j) = hy_x(thisnX-nXg,j)*(P(thisnX-nXg-1,j)+P(thisnX-nXg,j))/2.0;
      }   
      Xgrid.setXmaxFluxBC(FluxMx_x, P0);   
      Xgrid.setXmaxFluxBC(FluxMy_x, 0.0, 0.0);   
      Xgrid.setXmaxFluxBC(FluxNe_x, 0.0, 0.0);
   }   
   //Xgrid.communicate(FluxN_x);   
   //Xgrid.communicate(FluxMx_x);   
   //Xgrid.communicate(FluxMy_x);   
   //Xgrid.communicate(FluxEi_x);   
   //Xgrid.communicate(FluxEe_x);   
   //Xgrid.communicate(FluxNe_x);   
   
}

void updatePhysicalVars( const domainGrid&  Xgrid )
{
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   Vx  = Mx/N;
   Vy  = My/N;
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

   Pi  = (Ei - 0.5*(Vx*Mx + Vy*My))*(gamma0-1.0);
   Pe  = Ee*(gamma0-1.0);
   Pe_eff = Pe - mM*Pi;
   P = Pe + Pi; 
   Ti  = Pi/N;
   Te  = Pe/Ne;
   if(min(Te)<0.0) cout << " Te IS LESS THAN ZERO on proc " << procID << endl;
   if(min(Ti)<0.0) cout << " Ti IS LESS THAN ZERO on proc " << procID << endl;
   if(!procID) m_Xgrid.setXminBoundary(Ti,Tthresh/Tscale,0.0);

   m_Xgrid.InterpToCellCenter(Bycc,By);
   m_Xgrid.InterpToCellCenter(Bxcc,Bx);

   Cs = sqrt(gamma0*P/N + (Bxcc*Bxcc + Bycc*Bycc + Bz*Bz)/(N+1.0e-4));  // 1.0e-4 here for time step
   Cspeed_x  = abs(Vx) + Cs;                  // adv flux jacobian
   Cspeed_y  = abs(Vy) + Cs;                  // adv flux jacobian
   
   // interp electron density to edges
   //
   Xgrid.InterpToCellEdges(Ne_x,Ne,Ne,"C2",0);
   Xgrid.communicate(Ne_x);
   Xgrid.InterpToCellEdges(Ne_y,Ne,Ne,"C2",1);
   Xgrid.InterpToCellEdges(Ne_xy,Ne_x,Ne_x,"C2",1);   
   Xgrid.setZboundaryPeriodic(Ne_y);
   Xgrid.setZboundaryPeriodic(Ne_xy);

   // interp density to edges
   //
   Xgrid.InterpToCellEdges(N_x,N,N,"C2",0);
   Xgrid.InterpToCellEdges(N_y,N,N,"C2",1);
   Xgrid.communicate(N_x);
   Xgrid.setZboundaryPeriodic(N_y);

   // compute Hall and electron velocity at cell center
   //
   //Vhallx_y = -Li0*Jx/Ne_y;
   //Vhally_x = -Li0*Jy/Ne_x;
   Vhallx = -Li0*Jxcc/Ne;
   Vhally = -Li0*Jycc/Ne;
   Vhallz = -Li0*Jzcc/Ne;
   Vex = Vx + Vhallx;
   Vey = Vy + Vhally;
   Vez = Vy + Vhallz;

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
   computeJ0(Jx0,Jy0,Jz0,Bx,By,Bz);
   bool useJ0forJcc = false;
   if(useJ0forJcc) {   
      Xgrid.InterpToCellCenter(Jzcc,Jz0);
      Xgrid.InterpToCellCenter(Jxcc,Jx0);
      Xgrid.InterpToCellCenter(Jzcc,Jz0);
   } else {
      Xgrid.InterpToCellCenter(Jzcc,Jz);
      Xgrid.InterpToCellCenter(Jxcc,Jx);
      Xgrid.InterpToCellCenter(Jzcc,Jz);
   }
   Xgrid.setZboundaryPeriodic(Jxcc);
   Xgrid.communicate(Jycc);
   Xgrid.communicate(Jzcc);
   Xgrid.setZboundaryPeriodic(Jzcc);
   
   
   // compute electric field at cell center
   ///
   Xgrid.InterpToCellCenter(Excc,Ex);
   Xgrid.InterpToCellCenter(Eycc,Ey);
   Xgrid.InterpToCellCenter(Ezcc,Ez);
   Xgrid.setZboundaryPeriodic(Excc);
   Xgrid.communicate(Eycc);
   Xgrid.communicate(Ezcc);
   Xgrid.setZboundaryPeriodic(Ezcc);
   
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
   //Qiwall = max(Qie,0.0);
   //Qiwall = 0.01*(nue_spi + nue_neu)*N*(Ti-Tthresh/Tscale);

   Xgrid.InterpToCellEdges(eta_x,eta,eta,"C2",0);
   Xgrid.communicate(eta_x);
   Xgrid.InterpToCellEdges(eta_y,eta,eta,"C2",1);
   Xgrid.InterpToCellEdges(eta_xy,eta_x,eta_x,"C2",1);
   Xgrid.setZboundaryPeriodic(eta_y);
   Xgrid.setZboundaryPeriodic(eta_xy);


   // set heat flux coefficients
   //

}

void computeJ0( matrix2D<double>&  a_Jx0,
                matrix2D<double>&  a_Jy0,
                matrix2D<double>&  a_Jz0,
          const matrix2D<double>&  a_Bx,
          const matrix2D<double>&  a_By,
          const matrix2D<double>&  a_Bz )
{
   assert(a_Jx0.size1()==m_nZce);
   assert(a_Jy0.size0()==m_nXce);
   assert(a_Jz0.size0()==m_nXce);
   assert(a_Jz0.size1()==m_nYce);
   assert(a_Bx.size0()==m_nXce);
   assert(a_By.size1()==m_nYce);
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
      
  
   // comptue Jz0 = ( d(r*Bth)/dr - d(Br)/dtheta )/r
   //
   matrix2D<double> dBrdth;
   dBrdth.initialize(m_nXce,m_nYce,0.0);
   m_Xgrid.DDZ(dBrdth,a_Bx);
   m_Xgrid.setZboundaryPeriodic(dBrdth);
 
   m_Xgrid.DDX(a_Jz0,hy_y*a_By);
   m_Xgrid.communicate(a_Jz0);
   a_Jz0 -= dBrdth;
   a_Jz0 /= hy_xy;

    
   // comptue Jx0 = d(Bz)/dtheta/r
   //
   m_Xgrid.DDZ(a_Jx0, a_Bz);
   m_Xgrid.setZboundaryPeriodic(a_Jx0);
   a_Jx0 /= hy_y;

   
   // compute Jy0 = -d(Bz)/dr
   //
   m_Xgrid.DDX(a_Jy0, a_Bz);
   m_Xgrid.communicate(a_Jy0);
   a_Jy0 *= -1.0;
   

}

void computeHeatFluxes( const domainGrid&  Xgrid )
{
   // compute heat fluxes at cell edges
   // qex0 = -kappae*dTe/dx
   // qey0 = -kappae*dTe/dy
   // kappae = kappae0*Pe*taue
   //

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
   xe = wcescale*Bycc*taue*tscale; // Omega_ce*taue 
   xe2 = xe*xe;
   xe4 = xe2*xe2;
   chie_perp = (4.664*xe2 + 11.92)/(xe4 + 14.79*xe2 + 3.770); 

   // compute kappae at cell edges
   //
   //kappae = Ve0sq*Pe*taue*chie_perp;
   kappae = Ve0sq*Pe*taue*11.92/3.77;
   Xgrid.InterpToCellEdges(kappae_x,kappae,kappae,"C2",0);
   Xgrid.InterpToCellEdges(kappae_y,kappae,kappae,"C2",1);

   // compute gradient of Te
   //
   Xgrid.DDX(qex0,Te);
   Xgrid.DDZ(qey0,Te);
   Xgrid.communicate(qex0);
   Xgrid.setZboundaryPeriodic(qey0);

   // multily kappae by grad(Te) and negate
   //
   qex0 *= kappae_x;   
   qex0 *= -1.0;   
   qey0 *= kappae_y;   
   qey0 *= -1.0;   

   ///////////////////////////////////////////////////////
   
   // compute kappai at cell edges
   //
   kappai = mM*Ve0sq*Pi*taui*3.9;
   Xgrid.InterpToCellEdges(kappai_x,kappai,kappai,"C2",0);
   Xgrid.InterpToCellEdges(kappai_y,kappai,kappai,"C2",1);

   // compute gradient of Ti
   //
   Xgrid.DDX(qix0,Ti);
   Xgrid.DDZ(qiy0,Ti);
   Xgrid.communicate(qix0);
   Xgrid.setZboundaryPeriodic(qiy0);

   // multily kappae by grad(Te) and negate
   //
   qix0 *= kappai_x;   
   qix0 *= -1.0;   
   qiy0 *= kappai_y;   
   qiy0 *= -1.0;   


}

void computeIonViscStressTensor( const domainGrid&  Xgrid )
{
   // compute ion stress tensor for momentum and energy equations
   //
   
   // set ion viscousity coefficient
   //
   etaVis_ion = 0.96*Pi*min(taui,tauiMax);
   Xgrid.InterpToCellEdges(etaVis_ion_x,etaVis_ion,etaVis_ion,"C2",0);
   Xgrid.InterpToCellEdges(etaVis_ion_z,etaVis_ion,etaVis_ion,"C2",1);
   Xgrid.communicate(etaVis_ion_x);
   
   // compute ion strain tensor
   //
   computeStrainTensor( Pii_xx0, Pii_zz0, Pii_yy, Pii_xz0, Vx, Vy, Xgrid );
  
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

   // set electron viscosity coefficient
   //
   matrix2D<double> taueVis0(taue);
   taueVis0 = min(taueVis0,taueMax);
   taueVis0 = max(taueVis0,taueMin);
   etaVis_ele = 0.73*Pe*taueVis0;
   Xgrid.InterpToCellEdges(etaVis_ele_x,etaVis_ele,etaVis_ele_x,"C2",0);
   Xgrid.communicate(etaVis_ele_x);
   Xgrid.InterpToCellEdges(etaVis_ele_z,etaVis_ele,etaVis_ele_z,"C2",1);
   Xgrid.InterpToCellEdges(etaVis_ele_xz,etaVis_ele_x,etaVis_ele_xz,"C2",1);

   // compute ele strain tensor
   //
   //computeStrainTensorStag( Pie_xx0, Pie_zz0, Pie_yy, Pie_xz0, Vex_z, Vez_x);
   computeStrainTensorStag( Pie_xx0, Pie_zz0, Pie_yy, Pie_xz0, Vhallx_y, Vhally_x);
  
   // multiply ele strain tensor by -viscosity to get ele stress tensor
   //
   Pie_xx0 *= -etaVis_ele_xz;
   Pie_zz0 *= -etaVis_ele_xz;
   Pie_yy  *= -etaVis_ele_z;
   Pie_xz0 *= -etaVis_ele;

}

void computeStrainTensor( matrix2D<double>&  a_Pi_xx,
                          matrix2D<double>&  a_Pi_zz,
                          matrix2D<double>&  a_Pi_yy,
                          matrix2D<double>&  a_Pi_xz,
                    const matrix2D<double>&  a_Vx,
                    const matrix2D<double>&  a_Vy,
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
   dVxdx_x.initialize(m_nXce,nZcc,0.0); 
   dVxdz.initialize(nXcc,nZcc,0.0); 
   //
   Xgrid.DDX(dVxdx_x,a_Vx);
   Xgrid.DDZ(dVxdz,a_Vx);


   // compute Vy derivatives
   //
   matrix2D<double> dVydz_z, dVydx;
   dVydz_z.initialize(nXcc,m_nZce,0.0); 
   dVydx.initialize(nXcc,nZcc,0.0); 
   //
   Xgrid.DDZ(dVydz_z,a_Vy);
   Xgrid.DDX(dVydx,a_Vy);
   Xgrid.communicate(dVydx);
   

   // compute divergence of velocity
   //
   matrix2D<double> divUz, divUx, divU, divUx_x, divU_x, divU_z;
   divUz.initialize(nXcc,nZcc,0.0); 
   divUx.initialize(nXcc,nZcc,0.0); 
   divU.initialize(nXcc,nZcc,0.0); 
   divUx_x.initialize(m_nXce,nZcc,0.0); 
   divU_x.initialize(m_nXce,nZcc,0.0); 
   divU_z.initialize(nXcc,m_nZce,0.0); 

   // div at cell center
   //
   Xgrid.DDZ(divUz,a_Vy);
   Xgrid.setZboundaryPeriodic(divUz);
   Xgrid.DDX(divUx,a_Vx*hy_cc);
   divUx /= hy_cc;
   Xgrid.communicate(divUx);
   divU = divUx + divUz;

   // div at z-edge
   //
   Xgrid.InterpToCellEdges(divU_z,divUx,divUx,"C2",1);
   divU_z += dVydz_z;
   
   // div at x-edge
   //
   Xgrid.InterpToCellEdges(divU_x,divUz,divUz,"C2",0);
   Xgrid.DDX(divUx_x,hy_cc*a_Vx);
   divUx_x /= hy_x;
   vector<double> divUx0;
   divUx0.assign(nZcc,0.0);
   if(procID==0 && geometry0=="CYL" && XlowBC.compare("axis") == 0) {
      for (auto j=0; j<nZcc; j++) {
         divUx0.at(j) = 2.0*a_Vx(nXg,j)/hy_cc(nXg,j);
      }
      Xgrid.setXminFluxBC(divUx_x,divUx0);
   } 
   Xgrid.communicate(divUx_x);
   divU_x += divUx_x;
   

   // set components of stress tensor on cell edges
   // for momentum equation
   //
   a_Pi_xx = (2.0*dVxdx_x - 2.0/3.0*divU_x);
   a_Pi_zz = (2.0*dVydz_z - 2.0/3.0*divU_z);
   a_Pi_yy = (2.0*a_Vx/hy_cc  - 2.0/3.0*divU);
   a_Pi_xz = (dVxdz + dVydx);

}

void computeStrainTensorStag( matrix2D<double>&  a_Pi_xx,
                              matrix2D<double>&  a_Pi_zz,
                              matrix2D<double>&  a_Pi_yy,
                              matrix2D<double>&  a_Pi_xz,
                        const matrix2D<double>&  a_Vx,
                        const matrix2D<double>&  a_Vy )
{
   // compute strain tensor for viscosity
   // with input velocities defined on cell-edges
   //
   assert(a_Pi_xx.size0()==m_nXce && a_Pi_xx.size1()==m_nZce);
   assert(a_Pi_zz.size0()==m_nXce && a_Pi_zz.size1()==m_nZce);
   assert(a_Vx.size1()==m_nZce);
   assert(a_Vy.size0()==m_nXce);

   const int nXcc = m_Xgrid.nXcc;
   const int nZcc = m_Xgrid.nZcc;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
    
   // compute Vx derivatives
   //
   matrix2D<double> dVxdx_xz, dVxdz;
   dVxdx_xz.initialize(m_nXce,m_nZce,0.0); 
   dVxdz.initialize(nXcc,nZcc,0.0); 
   //
   m_Xgrid.DDX(dVxdx_xz,a_Vx);
   m_Xgrid.communicate(dVxdx_xz);
   m_Xgrid.DDZ(dVxdz,a_Vx);
   m_Xgrid.setZboundaryPeriodic(dVxdz);


   // compute Vy derivatives
   //
   matrix2D<double> dVydz_xz, dVydx;
   dVydz_xz.initialize(m_nXce,m_nZce,0.0); 
   dVydx.initialize(nXcc,nZcc,0.0); 
   //
   m_Xgrid.DDZ(dVydz_xz,a_Vy);
   m_Xgrid.DDX(dVydx,a_Vy); 
   m_Xgrid.communicate(dVydx);
   m_Xgrid.setZboundaryPeriodic(dVydz_xz);
   

   // compute divU on nodes
   //
   matrix2D<double> divUx_xz, divU_xz;
   divUx_xz.initialize(m_nXce,m_nZce,0.0); 
   divU_xz.initialize( m_nXce,m_nZce,0.0); 
   m_Xgrid.DDZ(divU_xz,a_Vy);
   m_Xgrid.setZboundaryPeriodic(divU_xz);
   m_Xgrid.DDX(divUx_xz,hy_y*a_Vx);
   divUx_xz /= hy_xy;
   vector<double> divUx0;
   divUx0.assign(m_nZce,0.0);
   if(procID==0 && geometry0=="CYL" && XlowBC.compare("axis") == 0) {
      for (auto j=0; j<m_nZce; j++) {
         divUx0.at(j) = 2.0*a_Vx(m_nXg,j)/hy_y(m_nXg,j);
      }
      m_Xgrid.setXminFluxBC(divUx_xz,divUx0);
   } 
   m_Xgrid.communicate(divUx_xz);
   divU_xz += divUx_xz;
   

   // compute divU on z-edge
   //
   matrix2D<double> divU_z;
   divU_z.initialize(nXcc,m_nZce,0.0); 
   m_Xgrid.InterpNodesToEdges(divU_z,divU_xz,1);
   

   // set components of stress tensor on cell edges
   // for momentum equation
   //
   a_Pi_xx = (2.0*dVxdx_xz  - 2.0/3.0*divU_xz);
   a_Pi_zz = (2.0*dVydz_xz  - 2.0/3.0*divU_xz);
   a_Pi_yy = (2.0*a_Vx/hy_y - 2.0/3.0*divU_z);
   a_Pi_xz = dVxdz + dVydx;

}

void computeIdealEatEdges() 
{
   // Eideal = -VxB
   //
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   matrix2D<double> VxBy_y, VyBx_x, VyBz, VxBz, VyBx_xy;
   VxBy_y.initialize(m_nXcc,m_nYce,0.0);
   VyBx_x.initialize(m_nXce,m_nYcc,0.0);
   VyBx_xy.initialize(m_nXce,m_nYce,0.0);
   VyBz.initialize(m_nXcc,m_nYcc,0.0);
   VxBz.initialize(m_nXcc,m_nYcc,0.0);

   // compute VxBy on y-faces and VyBx on x-faces
   //
   m_Xgrid.InterpToCellEdges(Vx_y,Vx,Vx,"C2",1);
   m_Xgrid.setZboundaryPeriodic(Vx_y);
   m_Xgrid.InterpToCellEdges(Vy_x,Vy,Vy,"C2",0);
   m_Xgrid.communicate(Vy_x);
   VxBy_y = Vx_y*By;
   VyBx_x = Vy_x*Bx;

   // compute VxBz and VyBz at cell center
   //
   VyBz = Vy*Bz;
   VxBz = Vx*Bz;

   if(advScheme0=="TVD") {
      m_Xgrid.computeFluxTVD( Eidealz,FluxL_xy,FluxR_xy,FluxRatio_xy,FluxLim_xy,
                              VxBy_y,abs(Vx_y),By,"minmod",0,Nsub);
      m_Xgrid.computeFluxTVD( VyBx_xy,FluxL_xy,FluxR_xy,FluxRatio_xy,FluxLim_xy,
                              VyBx_x,abs(Vy_x),Bx,"minmod",1,Nsub);
      //
      m_Xgrid.computeFluxTVD( Eidealx,FluxL_y,FluxR_y,FluxRatio_y,FluxLim_y,
                              VyBz,abs(Vy),Bz,"minmod",1,Nsub);
      m_Xgrid.computeFluxTVD( Eidealy,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                              hy_cc*VxBz,abs(Vx),hy_cc*Bz,"minmod",0,Nsub);
   }
   else {
      m_Xgrid.computeFluxTVDsimple(Eidealz,FluxL_xy,FluxR_xy,VxBy_y,abs(Vx_y),By,advScheme0,0);
      m_Xgrid.computeFluxTVDsimple(VyBx_xy,FluxL_xy,FluxR_xy,VyBx_x,abs(Vy_x),Bx,advScheme0,1);
      //
      m_Xgrid.computeFluxTVDsimple(Eidealx,FluxL_y,FluxR_y,VyBz,abs(Vy),Bz,advScheme0,1);
      m_Xgrid.computeFluxTVDsimple(Eidealy,FluxL_x,FluxR_x,VxBz,abs(Vx),Bz,advScheme0,0);
   }
   Eidealz = VyBx_xy - Eidealz;
   Eidealx *= -1.0;
   Eidealy /= hy_x;
   if(procID==0) {
      m_Xgrid.setXminFluxBC(Eidealz, 0.0, 0.0);
      m_Xgrid.setXminFluxBC(Eidealy, 0.0, 0.0);
      //m_Xgrid.setXminBoundary(Eidealx, 0.0, 1.0);
   }
   if(procID==numProcs-1) {
      m_Xgrid.setXmaxFluxBC(Eidealz, 0.0, 0.0);
      m_Xgrid.setXmaxFluxBC(Eidealy, 0.0, 0.0);
      //m_Xgrid.setXmaxBoundary(Eidealx, 0.0, 1.0);
   }
   m_Xgrid.communicate(Eidealx);
   m_Xgrid.communicate(Eidealy);
   m_Xgrid.communicate(Eidealz);
   m_Xgrid.setZboundaryPeriodic(Eidealx);
   m_Xgrid.setZboundaryPeriodic(Eidealy);
   m_Xgrid.setZboundaryPeriodic(Eidealz);

}

void computeHallEatEdges( const domainGrid&  Xgrid )
{
   // Ehall = Li0/Ne*J0xOmega 
   //
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);


   // get Hall velocity on opposite edges for upwinding
   //
   matrix2D<double> Vhallx_x, Vhally_y;
   Vhallx_x.initialize(m_nXce,m_nYcc,0.0);
   Vhally_y.initialize(m_nXcc,m_nYce,0.0);
   
   Xgrid.InterpToCellEdges(Vhallx_x,Vhallx,Vhallx,"C2",0);
   Xgrid.InterpToCellEdges(Vhally_y,Vhally,Vhally,"C2",1);
   Xgrid.communicate(Vhallx_x);
   Xgrid.communicate(Vhally_y);

   // upwind By to edges
   //
   Xgrid.InterpCellToEdges(Ehallx,Bycc,Vhally_y,advSchemeHall0,1);
   Xgrid.InterpCellToEdges(Ehallz,Bycc,Vhallx_x,advSchemeHall0,0);
   Xgrid.communicate(Ehallx);
   Xgrid.communicate(Ehallz);
   Xgrid.setZboundaryPeriodic(Ehallx);
   Xgrid.setZboundaryPeriodic(Ehallz);

   // put it all together
   //
   Ehallx *= Vhally_y;
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
   eJxFluxZ_cc = Le0sq*(Vy*Jxcc + Jzcc*Vex);
   eJzFluxX_cc = Le0sq*(Vx*Jzcc + Jxcc*Vez)*hy_cc;
   eJzFluxZ_cc = Le0sq*(Vy*Jzcc + Jzcc*Vez);

   // upwind cell center fluxes to cell faces
   //
   const int Nsub = 2;
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(eJxFlux_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           eJxFluxX_cc,Cspeed_x, hy_cc*Le0sq*Jxcc,"minmod",0,Nsub);
      Xgrid.computeFluxTVD(eJzFlux_x,FluxL_x,FluxR_x,FluxRatio_x,FluxLim_x,
                           eJzFluxX_cc,Cspeed_x, hy_cc*Le0sq*Jzcc,"minmod",0,Nsub);
      
      Xgrid.computeFluxTVD(eJxFlux_z,FluxL_y,FluxR_y,FluxRatio_y,FluxLim_y,
                           eJxFluxZ_cc,Cspeed_y, Le0sq*Jxcc,"minmod",1,Nsub);
      Xgrid.computeFluxTVD(eJzFlux_z,FluxL_y,FluxR_y,FluxRatio_y,FluxLim_y,
                           eJzFluxZ_cc,Cspeed_y, Le0sq*Jzcc,"minmod",1,Nsub);
   }
   else {
      Xgrid.computeFluxTVDsimple(eJxFlux_x,FluxL_x,FluxR_x,eJxFluxX_cc,
                                 Cspeed_x,hy_cc*Le0sq*Jxcc,advScheme0,0);
      Xgrid.computeFluxTVDsimple(eJzFlux_x,FluxL_x,FluxR_x,eJzFluxX_cc,
                                 Cspeed_x,hy_cc*Le0sq*Jxcc,advScheme0,0);
        
      Xgrid.computeFluxTVDsimple(eJxFlux_z,FluxL_y,FluxR_y,eJxFluxZ_cc, 
                                 Cspeed_y,Le0sq*Jxcc,advScheme0,1);
      Xgrid.computeFluxTVDsimple(eJzFlux_z,FluxL_y,FluxR_y,eJzFluxZ_cc,
                                 Cspeed_y,Le0sq*Jzcc,advScheme0,1);
   }
   
   // set flux BCs
   //
   if(procID==0) {
      Xgrid.setXminFluxBC(eJxFlux_x, 0.0, 0.0);
      Xgrid.setXminFluxBC(eJzFlux_x, 0.0, 0.0);
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
         divJxstress_cc(i,j) = (eJxFlux_x(i+1,j) - eJxFlux_x(i,j))/hy_cc(i,j)/dX
                             + (eJxFlux_z(i,j+1) - eJxFlux_z(i,j))/dZ;
         divJzstress_cc(i,j) = (eJzFlux_x(i+1,j) - eJzFlux_x(i,j))/hy_cc(i,j)/dX
                             + (eJzFlux_z(i,j+1) - eJzFlux_z(i,j))/dZ;
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
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   matrix2D<double> gradPex, gradPez;
   gradPez.initialize(m_nXcc,m_nZcc,0.0);
   gradPex.initialize(m_nXcc,m_nZcc,0.0);
   
   Xgrid.DDX(gradPex,Pe_eff);
   Xgrid.DDZ(gradPez,Pe_eff);
   Xgrid.communicate(gradPex);   
   
   Xgrid.InterpToCellEdges(Egradz,gradPez,gradPez,"C2",0);
   Xgrid.InterpToCellEdges(Egradx,gradPex,gradPex,"C2",1);
   Xgrid.communicate(Egradz);
   
   Egradx *= -Li0/Ne_y;
   Egradz *= -Li0/Ne_x;

}

void computeJouleHeatingTerms( const domainGrid&  Xgrid )
{
   const int nXcc = Xgrid.nXcc;
   const int nZcc = Xgrid.nZcc;
   const int nXg  = Xgrid.nXg;
   const int nZg  = Xgrid.nZg;

   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   matrix2D<double> dPedx, dPedy;
   dPedx.initialize(nXcc,nZcc,0.0);
   dPedy.initialize(nXcc,nZcc,0.0);


   // set momentum density source terms... move this
   //
   JcrossBx =  Jycc*Bz   - Jzcc*Bycc;
   JcrossBy =  Jzcc*Bxcc - Jxcc*Bz;

  
   //  calculate work source terms for energy equations
   //
   if(useIonScaleCorrections) {
      NeUdotE = Ne*(Vy*(Eycc - eta*Jycc) + Vx*(Excc - eta*Jxcc));
      NeUdotE /= Li0;
   } 
   else {
      Xgrid.DDX(dPedx,Pe);
      Xgrid.communicate(dPedx);   
      Xgrid.DDZ(dPedy,Pe);
      dPedy /= hy_cc;
      Xgrid.setZboundaryPeriodic(dPedy);   
      NeUdotE = Vy*(Jzcc*Bxcc - Jxcc*Bz - dPedy) + Vx*(Jycc*Bz - Jzcc*Bycc - dPedx);
   }
   JdotE = Jxcc*Excc + Jycc*Eycc + Jzcc*Ezcc;
   

   ////////////
   //
   // compute divergence of heat flux
   //
   matrix2D<double> dummydiv;
   dummydiv.initialize(nXcc,nZcc,0.0); 
   
   Xgrid.DDX(dummydiv,hy_x*qex);
   Xgrid.DDZ(divqe,qey);
   divqe += dummydiv;
   divqe /= hy_cc;
   
   Xgrid.DDX(dummydiv,hy_x*qix);
   Xgrid.DDZ(divqi,qiy);
   divqi += dummydiv;
   divqi /= hy_cc;
   
   
   if(tauiMax>0.0) {
     
      // interp Pii_xz to cell edges in x and z
      //
      matrix2D<double> Pii_xz_x, Pii_xz_z;
      Pii_xz_x.initialize(m_nXce,nZcc,0.0);
      Pii_xz_z.initialize(nXcc,m_nZce,0.0);
      Xgrid.InterpToCellEdges(Pii_xz_x,Pii_xz,Pii_xz,"C2",0);
      Xgrid.InterpToCellEdges(Pii_xz_z,Pii_xz,Pii_xz,"C2",1);

      ////////////
      //
      // compute divergence of stress tensor for x-momentum equation
      //
      Xgrid.DDX(divPii_x,Pii_xx*hy_x);
      divPii_x /= hy_cc;
      Xgrid.DDZ(dummydiv,Pii_xz_z);
      divPii_x += dummydiv;
      if(geometry0=="CYL") divPii_x -= Pii_yy/hy_cc;
      Xgrid.communicate(divPii_x);
   
 
      ////////////
      //
      // compute divergence of stress tensor for z-momentum equation
      //
      Xgrid.DDX(divPii_z,Pii_xz_x*hy_x);
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
      matrix2D<double> Vx_x, Vx_z, Vy_x, Vy_z;
      Vx_x.initialize(m_nXce,nZcc,0.0); 
      Vy_x.initialize(m_nXce,nZcc,0.0); 
      Vx_z.initialize(nXcc,m_nZce,0.0); 
      Vy_z.initialize(nXcc,m_nZce,0.0); 
      Xgrid.InterpToCellEdges(Vx_x,Vx,Vx,"C2",0);
      Xgrid.InterpToCellEdges(Vy_x,Vy,Vy,"C2",0);
      Xgrid.InterpToCellEdges(Vx_z,Vx,Vx,"C2",1);
      Xgrid.InterpToCellEdges(Vy_z,Vy,Vy,"C2",1);
   
      // set components of stress tensor for energy equations
      //
      UidotPii_x = Vx_x*Pii_xx   + Vy_x*Pii_xz_x; 
      UidotPii_z = Vx_z*Pii_xz_z + Vy_z*Pii_zz; 

      // take divergence for energy equation
      //
      Xgrid.DDX(divUidotPii,UidotPii_x*hy_x);
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
   m_Xgrid.DDX(divPiex_z,Pie_xx*hy_xy);
   divPiex_z /= hy_y;
   if(geometry0=="CYL") divPiex_z -= Pie_yy/hy_y;
   m_Xgrid.DDZ(dummy_onz,Pie_xz);
   divPiex_z += dummy_onz;
   m_Xgrid.communicate(divPiex_z);
   m_Xgrid.setZboundaryPeriodic(divPiex_z); 
  
   // compute divergence of viscous stress tensor for Jz equation
   //
   m_Xgrid.DDX(divPiez_x,Pie_xz*hy_cc);
   divPiez_x /= hy_x;
   vector<double> divPiexz0;
   divPiexz0.assign(m_nZcc,0.0);
   if(procID==0 && geometry0=="CYL" && XlowBC.compare("axis") == 0) {
      for (auto j=0; j<m_nZcc; j++) {
         divPiexz0.at(j) = 2.0*Pie_xz(m_nXg,j)/hy_cc(m_nXg,j);
      }
      m_Xgrid.setXminFluxBC(divPiez_x,divPiexz0);
   } 
   m_Xgrid.communicate(divPiez_x);
   
   m_Xgrid.DDZ(dummy_onx,Pie_zz);
   m_Xgrid.setZboundaryPeriodic(dummy_onx); 
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
   const int nXcc = Xgrid.Xcc.size();
   const int nZcc = Xgrid.Zcc.size();
   const int nXg = Xgrid.nXg;
   const int nZg = Xgrid.nZg;
   
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);
   
   deltaN.initialize(nXcc,nZcc,0.0);
   
   N.initialize( nXcc,nZcc,0.0);
   Mx.initialize( nXcc,nZcc,0.0);
   My.initialize( nXcc,nZcc,0.0);
   Ei.initialize(nXcc,nZcc,0.0);
   Ee.initialize(nXcc,nZcc,0.0);
   //
   Nold.initialize( nXcc,nZcc,0.0);
   Mxold.initialize( nXcc,nZcc,0.0);
   Myold.initialize( nXcc,nZcc,0.0);
   Eiold.initialize(nXcc,nZcc,0.0);
   Eeold.initialize(nXcc,nZcc,0.0);
   
   //
   // By stag in y
   //
   By.initialize(   m_nXcc,m_nYce,0.0);
   Byold.initialize(m_nXcc,m_nYce,0.0);
   Bycc.initialize( m_nXcc,m_nYcc,0.0);
   //
   // Bx stag in x
   //
   Bx.initialize(   m_nXce,m_nYcc,0.0);
   Bxold.initialize(m_nXce,m_nYcc,0.0);
   Bxcc.initialize( m_nXcc,m_nYcc,0.0);
   //
   // Bz cell center
   //
   Bz.initialize(   m_nXcc,m_nYcc,0.0);
   Bzold.initialize(m_nXcc,m_nYcc,0.0);
   //
   // Ez and Jz stag in x and y
   //
   Ez.initialize(   m_nXce,m_nYce,0.0);
   Jz.initialize(   m_nXce,m_nYce,0.0);
   Ezold.initialize(m_nXce,m_nYce,0.0);
   Jzold.initialize(m_nXce,m_nYce,0.0);
   Ezcc.initialize( m_nXcc,m_nYcc,0.0);
   Jzcc.initialize( m_nXcc,m_nYcc,0.0);
   Jz0.initialize(  m_nXce,m_nYce,0.0);
   Jz0stagTime.initialize( m_nXce,m_nYce,0.0);
   eta_xy.initialize(  m_nXce,m_nYce,0.0);
   Ne_xy.initialize(   m_nXce,m_nYce,0.0);
   Eidealz.initialize( m_nXce,m_nYce,0.0);
   //
   FluxRatio_xy.initialize(m_nXce,m_nYce,0.0);
   FluxLim_xy.initialize(  m_nXce,m_nZce,0.0);
   FluxR_xy.initialize(    m_nXce,m_nZce,0.0);
   FluxL_xy.initialize(    m_nXce,m_nYce,0.0);
   Vx_y.initialize(        m_nXcc,m_nYce,0.0);
   Vy_x.initialize(        m_nXce,m_nYcc,0.0);
   //
   // Ex and Jx stag in y
   //
   Ex.initialize(   m_nXcc,m_nYce,0.0);
   Jx.initialize(   m_nXcc,m_nYce,0.0);
   Exold.initialize(m_nXcc,m_nYce,0.0);
   Jxold.initialize(m_nXcc,m_nYce,0.0);
   Excc.initialize( m_nXcc,m_nYcc,0.0);
   Jxcc.initialize( m_nXcc,m_nYcc,0.0);
   Jx0.initialize(  m_nXcc,m_nYce,0.0);
   Jx0stagTime.initialize( m_nXcc,m_nYce,0.0);
   eta_y.initialize(m_nXcc,m_nYce,0.0);
   Ne_y.initialize( m_nXcc,m_nYce,0.0);
   Eidealx.initialize(m_nXcc,m_nYce,0.0); 
   //
   // Ey and Jy stag in x
   //
   Ey.initialize(   m_nXce,m_nYcc,0.0);
   Jy.initialize(   m_nXce,m_nYcc,0.0);
   Eyold.initialize(m_nXce,m_nYcc,0.0);
   Jyold.initialize(m_nXce,m_nYcc,0.0);
   Eycc.initialize( m_nXcc,m_nYcc,0.0);
   Jycc.initialize( m_nXcc,m_nYcc,0.0);
   Jy0.initialize(  m_nXce,m_nYcc,0.0);
   Jy0stagTime.initialize( m_nXce,m_nYcc,0.0);
   eta_x.initialize(m_nXce,m_nYcc,0.0);
   Ne_x.initialize( m_nXce,m_nYcc,0.0);
   Eidealy.initialize(m_nXce,m_nYcc,0.0); 
   //
   //
   //

   P.initialize(  nXcc,nZcc,0.0);
   Pi.initialize( nXcc,nZcc,0.0);
   Pe.initialize( nXcc,nZcc,0.0);
   Ti.initialize( nXcc,nZcc,0.0);
   Te.initialize( nXcc,nZcc,0.0);
   eta.initialize(nXcc,nZcc,0.0);
   Cs.initialize( nXcc,nZcc,0.0);
   Vx.initialize( nXcc,nZcc,0.0);
   Vy.initialize( nXcc,nZcc,0.0);
   Qvisc.initialize( nXcc,nZcc,0.0);
   Cspeed_x.initialize(nXcc,nZcc,0.0);
   Cspeed_y.initialize(nXcc,nZcc,0.0);
   //
   Qie.initialize(nXcc,nZcc,0.0);
   Qiwall.initialize(nXcc,nZcc,0.0);
   NeUdotE.initialize(nXcc,nZcc,0.0);
   JdotE.initialize(nXcc,nZcc,0.0);
   JcrossBx.initialize(nXcc,nZcc,0.0);
   JcrossBy.initialize(nXcc,nZcc,0.0);
   taue.initialize(nXcc,nZcc,0.0);
   taui.initialize(nXcc,nZcc,0.0);
   nue_spi.initialize(nXcc,nZcc,0.0);
   nue_vac.initialize(nXcc,nZcc,0.0);
   
   N_x.initialize(m_nXce,nZcc,0.0);
   N_y.initialize(nXcc,m_nZce,0.0);
 

   //
   FluxN_x.initialize(m_nXce,nZcc,0.0);
   FluxMx_x.initialize(m_nXce,nZcc,0.0);
   FluxMy_x.initialize(m_nXce,nZcc,0.0);
   FluxEe_x.initialize(m_nXce,nZcc,0.0);
   FluxEi_x.initialize(m_nXce,nZcc,0.0);
   FluxNe_x.initialize(m_nXce,nZcc,0.0);
   FluxRatio_x.initialize(m_nXce,nZcc,0.0);
   FluxLim_x.initialize(m_nXce,nZcc,0.0);
   FluxR_x.initialize(m_nXce,nZcc,0.0);
   FluxL_x.initialize(m_nXce,nZcc,0.0);
   //
   hy_cc.initialize(m_nXcc,m_nYcc,1.0);
   hy_x.initialize( m_nXce,m_nYcc,1.0);
   hy_y.initialize( m_nXcc,m_nYce,1.0);
   hy_xy.initialize(m_nXce,m_nYce,1.0);
   

   FluxN_y.initialize( m_nXcc,m_nYce,0.0);
   FluxNe_y.initialize(m_nXcc,m_nYce,0.0);
   FluxMx_y.initialize(m_nXcc,m_nYce,0.0);
   FluxMy_y.initialize(m_nXcc,m_nYce,0.0);
   FluxEe_y.initialize(m_nXcc,m_nYce,0.0);
   FluxEi_y.initialize(m_nXcc,m_nYce,0.0);
   //
   FluxRatio_y.initialize(m_nXcc,m_nYce,0.0);
   FluxLim_y.initialize(m_nXcc,m_nYce,0.0);
   FluxR_y.initialize(m_nXcc,m_nYce,0.0);
   FluxL_y.initialize(m_nXcc,m_nYce,0.0);

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
   Ehallx.initialize(nXcc,m_nZce,0.0); 
   Ehallz.initialize(m_nXce,nZcc,0.0); 
   Egradx.initialize(nXcc,m_nZce,0.0); 
   Egradz.initialize(m_nXce,nZcc,0.0); 
   Vhallx.initialize(m_nXcc,m_nYcc,0.0); 
   Vhally.initialize(m_nXcc,m_nYcc,0.0); 
   Vhallz.initialize(m_nXcc,m_nYcc,0.0); 
   //Vhallx_y.initialize(m_nXcc,m_nYce,0.0);
   //Vhally_x.initialize(m_nXce,m_nYcc,0.0);
   Vex.initialize(nXcc,nZcc,0.0); 
   Vey.initialize(nXcc,nZcc,0.0); 
   Vez.initialize(nXcc,nZcc,0.0); 
   divJxstress.initialize(nXcc,m_nZce,0.0); 
   divJzstress.initialize(m_nXce,nZcc,0.0); 
   Pe_eff.initialize(nXcc,nZcc,0.0); 
   eJxFlux_x.initialize(m_nXce,nZcc,0.0); 
   eJzFlux_x.initialize(m_nXce,nZcc,0.0); 
   eJxFlux_z.initialize(nXcc,m_nZce,0.0); 
   eJzFlux_z.initialize(nXcc,m_nZce,0.0); 

   // electron heat flux
   //
   qex.initialize(m_nXce,m_nYcc,0.0);
   qexold.initialize(m_nXce,m_nYcc,0.0);
   qex0.initialize(m_nXce,m_nYcc,0.0);
   kappae_x.initialize(m_nXce,m_nYcc,0.0);
   qey.initialize(   m_nXcc,m_nYce,0.0);
   qeyold.initialize(m_nXcc,m_nYce,0.0);
   qey0.initialize(  m_nXcc,m_nYce,0.0);
   kappae_y.initialize(m_nXcc,m_nYce,0.0);
   //
   divqe.initialize(m_nXcc,m_nYcc,0.0);
   kappae.initialize(m_nXcc,m_nYcc,0.0);
   chie_perp.initialize(m_nXcc,m_nYcc,0.0);
   
   // ion heat flux
   //
   qix.initialize(m_nXce,m_nYcc,0.0);
   qixold.initialize(m_nXce,m_nYcc,0.0);
   qix0.initialize(m_nXce,m_nYcc,0.0);
   kappai_x.initialize(m_nXce,m_nYcc,0.0);
   qiy.initialize(   m_nXcc,m_nYce,0.0);
   qiyold.initialize(m_nXcc,m_nYce,0.0);
   qiy0.initialize(  m_nXcc,m_nYce,0.0);
   kappai_y.initialize(m_nXcc,m_nYce,0.0);
   //
   divqi.initialize(m_nXcc,m_nYcc,0.0);
   kappai.initialize(m_nXcc,m_nYcc,0.0);

   
   // ion viscosity
   //
   Pii_xxold.initialize(m_nXce,nZcc,0.0);
   Pii_xx0.initialize(m_nXce,nZcc,0.0);
   Pii_xx.initialize(m_nXce,nZcc,0.0);
   Pii_xzold.initialize(nXcc,nZcc,0.0);
   Pii_xz0.initialize(nXcc,nZcc,0.0);
   Pii_xz.initialize(nXcc,nZcc,0.0);
   Pii_zzold.initialize(nXcc,m_nZce,0.0);
   Pii_zz0.initialize(nXcc,m_nZce,0.0);
   Pii_zz.initialize(nXcc,m_nZce,0.0);
   Pii_yy.initialize(nXcc,nZcc,0.0);

   divPii_x.initialize(nXcc,nZcc,0.0);
   divPii_z.initialize(nXcc,nZcc,0.0);
   divUidotPii.initialize(nXcc,nZcc,0.0);
   etaVis_ion.initialize(nXcc,nZcc,0.0);
   etaVis_ion_x.initialize(m_nXce,nZcc,0.0); 
   etaVis_ion_z.initialize(nXcc,m_nZce,0.0); 
   UidotPii_x.initialize(m_nXce,nZcc,0.0);
   UidotPii_z.initialize(nXcc,m_nZce,0.0);
   

   // ele viscosity
   //
   Pie_xxold.initialize(m_nXce,m_nZce,0.0);
   Pie_xx0.initialize(  m_nXce,m_nZce,0.0);
   Pie_xx.initialize(   m_nXce,m_nZce,0.0);
    
   Pie_xzold.initialize(nXcc,nZcc,0.0);
   Pie_xz0.initialize(  nXcc,nZcc,0.0);
   Pie_xz.initialize(   nXcc,nZcc,0.0);
   
   Pie_zzold.initialize(m_nXce,m_nZce,0.0);
   Pie_zz0.initialize(  m_nXce,m_nZce,0.0);
   Pie_zz.initialize(   m_nXce,m_nZce,0.0);
   Pie_yy.initialize(   nXcc,m_nZce,0.0);

   etaVis_ele.initialize(nXcc,nZcc,0.0);
   etaVis_ele_x.initialize(m_nXce,m_nZcc,0.0); 
   etaVis_ele_z.initialize(m_nXcc,m_nZce,0.0); 
   etaVis_ele_xz.initialize(m_nXce,m_nZce,0.0); 
   divPiex_z.initialize(nXcc,m_nZce,0.0);
   divPiez_x.initialize(m_nXce,nZcc,0.0);


}

void addMembersToDataFile( HDF5dataFile&  dataFile )
{
   dataFile.add(N, "N", 1);      // density 
   dataFile.add(Mx, "Mx", 1);      // momentum density 
   dataFile.add(My, "My", 1);      // momentum density 
   dataFile.add(Bx, "Bx", 1);      // magnetic field
   dataFile.add(By, "By", 1);      // magnetic field
   dataFile.add(Bz, "Bz", 1);      // magnetic field
   dataFile.add(Bxcc, "Bxcc", 1);      // magnetic field
   dataFile.add(Bycc, "Bycc", 1);      // magnetic field
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
   dataFile.add(qey,"qey", 1);  
   dataFile.add(qex0,"qex0", 1); 
   dataFile.add(qey0,"qey0", 1);  
   dataFile.add(chie_perp,"chie_perp", 1);  
   dataFile.add(kappae,"kappae", 1);  
   dataFile.add(divqe,"divqe", 1);  
   
   // ion heat flux
   //
   dataFile.add(qix,"qix", 1); 
   dataFile.add(qiy,"qiy", 1);  
   dataFile.add(qix0,"qix0", 1); 
   dataFile.add(qiy0,"qiy0", 1);  
   dataFile.add(kappai,"kappai", 1);  
   dataFile.add(divqi,"divqi", 1);  
   
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
   dataFile.add(Ne_xy,"Ne_xy", 1); 


   //dataFile.add(Vx, "Vx", 1);      // velocity
   //dataFile.add(Vy, "Vy", 1);      // velocity
   dataFile.add(Jx, "Jx", 1);     // x-current density
   dataFile.add(Jy, "Jy", 1);     // y-current density
   dataFile.add(Jz, "Jz", 1);     // z-current density
   dataFile.add(Jxcc, "Jxcc", 1); // x-current density at cell-center
   dataFile.add(Jycc, "Jycc", 1); // y-current density at cell-center
   dataFile.add(Jzcc, "Jzcc", 1); // z-current density at cell-center
   dataFile.add(Jx0, "Jx0", 1);   // curl of B
   dataFile.add(Jy0, "Jy0", 1);   // curl of B
   dataFile.add(Jz0, "Jz0", 1);   // curl of B
   dataFile.add(Ex, "Ex", 1);          // x-electric field
   dataFile.add(Ey, "Ey", 1);          // y-electric field
   dataFile.add(Ez, "Ez", 1);          // z-electric field
   dataFile.add(Excc, "Excc", 1);          // x-electric field
   dataFile.add(Eycc, "Eycc", 1);          // y-electric field
   dataFile.add(Ezcc, "Ezcc", 1);          // z-electric field
   dataFile.add(Eidealx, "Eidealx", 1);  // x-ideal electric field
   dataFile.add(Eidealy, "Eidealy", 1);  // y-ideal electric field
   dataFile.add(Eidealz, "Eidealz", 1);  // z-ideal electric field
   dataFile.add(Ehallz, "Ehallz", 1);  // z-Hall electric field
   dataFile.add(Ehallx, "Ehallx", 1);  // x-Hall electric field
   dataFile.add(Egradz, "Egradz", 1);  // z-grad electric field
   dataFile.add(Egradx, "Egradx", 1);  // x-grad electric field
   dataFile.add(Vhallx, "Vhallx", 1);  // Vhallx = - Li0/Ne*Jx
   dataFile.add(Vhally, "Vhally", 1);  // Vhally = - Li0/Ne*Jy
   //dataFile.add(Vhallz, "Vhallz", 1);  // Vhallz = - Li0/Ne*Jx
   dataFile.add(Vex, "Vex", 1);    // Vex = Vx - Li0/Ne*Jx
   dataFile.add(Vez, "Vez", 1);    // Vey = Vy - Li0/Ne*Jz
   //dataFile.add(eJxFlux_x, "eJxFlux_x", 1);
   //dataFile.add(eJxFlux_z, "eJxFlux_z", 1);
   //dataFile.add(eJzFlux_x, "eJzFlux_x", 1);
   //dataFile.add(eJzFlux_z, "eJzFlux_z", 1);
   //dataFile.add(divJxstress, "divJxstress", 1);  // x-comp of Ohm's law stress tensor
   //dataFile.add(divJzstress, "divJzstress", 1);  // z-comp of Ohm's law stress tensor
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
   dataFile.add(FluxMy_x, "FluxMy_x", 1);  
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
   dataFile.add(hy_x,"hy_x",0);
   dataFile.add(hy_y,"hy_y",0);
   dataFile.add(hy_xy,"hy_xy",0);

}

void parseInputFile( const domainGrid& Xgrid, const Json::Value& a_root )
{ 
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
      Json::Value taueMinVal = Phys.get("taueMin",defValue);
     
  
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
            for (auto i=0; i<m_nXce; i++) {
               for (auto j=0; j<m_nZce; j++) {
                  if(i<nXcc && j<nZcc) hy_cc(i,j) = Xgrid.Xcc.at(i);
                  if(j<nZcc) hy_x(i,j) = Xgrid.Xce2.at(i);
                  if(i<nXcc) hy_y(i,j) = Xgrid.Xcc.at(i);
                  hy_xy(i,j) = Xgrid.Xce2.at(i);
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
      
      taueMin = taueMinVal.asDouble();
      if(procID==0) cout << "min taue for visc = " << taueMin << endl;

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
   double Cmax_x, Cmax_y, nue_izn_max;
   double Vhall_max_x, Vhall_max_y;
   //Cchar = abs(Vx)+Cs;
   //Cmax = max(Cchar);
   Cmax_x = max(Cspeed_x);
   Cmax_y = max(Cspeed_y);
   Vhall_max_x = max(abs(Vhallx));
   Vhall_max_y = max(abs(Vhally));
   nue_izn_max = max(abs(nue_izn));
   //cout << "Cmax_x = " << Cmax_x << endl;
   //cout << "abs(Vx) = " << max(abs(Vx))<< endl;
   //cout << "Cs = " << max(Cs)<< endl;

   const double dX = Xgrid.dX;
   const double dZ = Xgrid.dZ;

   double dtCFL_sound_x = m_dX/Cmax_x;
   double dtCFL_sound_y = m_dY/Cmax_y;
   double dtCFL_sound = min(dtCFL_sound_x,dtCFL_sound_y);
   
   double dtCFL_hall_x = m_dX/Vhall_max_x;
   double dtCFL_hall_y = m_dY/Vhall_max_y;
   double dtCFL_hall = min(dtCFL_hall_x,dtCFL_hall_y);

   double dt_izn = 1.0/nue_izn_max;

   double dtCFL_light = m_dX*m_dY/sqrt(m_dX*m_dX + m_dY*m_dY)*sqrt(delta0);
   //double dtCFL_light = m_dX*sqrt(delta0);

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

