/***
 * 
 * physics module for 2D resistive MHD pinch
 *
 * dN/dt  + d(Mx)/dx = 0
 * dMx/dt + d(Mx*Ux + P )/dx = J*B
 * dS/dt  + d(S*Ux)/dx = eta*(gamma-1)/N^(gamma-1)*J^2  
 * dBy/dt + d(-Ex)/dx = 0
 * dEz/dt + d(-By)/dx = -J
 *
 * J = [E+U*B]/eta
 * S = P/N^(gamma-1)
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

string advScheme0;    // advection differencing scheme
double gamma0;        // adiabatic coefficient
double eta0=0.08;     // resistivity (may need to decrease dt more if eta0 really low)
double delta0=1.0e-4; // relaxation const (v/c)^2
double B0 = 0.0;      // boundary value of magnetic field
int Nsub;             // time-solver subcycle steps
vector<double> N, M, S, B, Ez;   // time-evolving variables
vector<double> eta, Cs, V, P, T, J, J0; // derived variables
vector<double> etace, Jcc, VBce;
vector<double> Nold, Mold, Sold, Bold, Ezold, Ezhalf;
vector<double> FluxRatio, FluxLim;
vector<double> FluxR, FluxL;  // flux at cell-edges   
vector<double> FluxN, FluxM, FluxS, FluxB, FluxEz;

//vector<vector<double>> matrixA(20, vector<double>(21));
vector<vector<double>> matW, matrixA, matrixB, matrixC;


// Set lower threshold values for N, T, and S
double Nthresh=1.0e-4, Tthresh=2.5e-2, Sthresh;

void computeFluxes(const domainGrid&, const int);
void setXminBoundary(vector<double>&, const double, const double);
void setXminBoundary(vector<vector<double>>&, const double, const double);
void setXminExtrap(vector<double>&);
void setXminExtrap(vector<vector<double>>&);
void setXmaxBoundary(vector<double>&, const double, const double);
void setXmaxBoundary(vector<vector<double>>&, const double, const double);

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

   const int nZcc = Xgrid.Zcc.size(); // number of z-cells
   const int nZce = Xgrid.Zce.size();

   //matrixA.assign(nXce,N);
   vector<double> vectorA;
   vectorA.assign(nZcc,1.0);
   //matrixA.assign(Zcc,Xgrid.Xcc);
   matW.assign(nXcc,vector<double>(nZcc)); // cell-center
   matrixA.assign(nXcc,vector<double>(nZcc)); // cell-center
   matrixB.assign(nXcc,vector<double>(nZce)); // stag-Z
   matrixC.assign(nXcc,vector<double>(nZce)); // stag-X
   //matrixA.assign(nXcc,Xgrid.Zcc); // cell-center
   for (auto i=0; i<nXcc; i++) {
      for (auto j=0; j<nZcc; j++) {
         matrixA[i][j] = Xgrid.Xcc[i]*Xgrid.Zcc[j];
      }
   }
   //Xgrid.DDX(matrixC,matrixA);
   //Xgrid.InterpToCellCenter(matrixC,matrixA);
   Xgrid.InterpToCellEdges(matrixC,matrixA,matrixA,"C2",1);
   Xgrid.communicate(matrixC);
   if(procID==0) setXminBoundary(matrixC, 0.0, 1.0);   
   //if(procID==0) setXminExtrap(matrixC);   
   if(procID==numProcs-1) setXmaxBoundary(matrixC, 0.0, 1.0);   
   
   //matrixA = 1.0/matrixA; 
   //matrixA = 3.0+matrixA;
   //matrixA = exp(matrixA);
   //double minA;
   //minA = min(matrixA);
   //cout << "minA = " << minA << endl;

   N.assign(nXcc,0.0);
   Nold.assign(nXcc,0.0);
   M.assign(nXcc,0.0);
   Mold.assign(nXcc,0.0);
   S.assign(nXcc,0.0);
   Sold.assign(nXcc,0.0);
   B.assign(nXcc,0.0);
   Bold.assign(nXcc,0.0);
   P.assign(nXcc,0.0);
   T.assign(nXcc,0.0);
   eta.assign(nXcc,0.0);
   Cs.assign(nXcc,0.0);
   V.assign(nXcc,0.0);
   J.assign(nXce,0.0); // J defined at cell edges
   Jcc.assign(nXcc,0.0);
   etace.assign(nXce,0.0);
   VBce.assign(nXce,0.0);
   J0.assign(nXce,0.0);
   // Ez is defined on cell edges
   Ez.assign(nXce,0.0);
   Ezold.assign(nXce,0.0);
   Ezhalf.assign(nXce,0.0);
   //
   FluxRatio.assign(nXce,0.0);
   FluxLim.assign(nXce,0.0);
   FluxN.assign(nXce,0.0);
   FluxM.assign(nXce,0.0);
   FluxS.assign(nXce,0.0);
   FluxB.assign(nXce,0.0);
   FluxEz.assign(nXcc,0.0); // Flux for Ez on cell-center
   FluxR.assign(nXce,0.0);
   FluxL.assign(nXce,0.0);
   //
   const Json::Value defValue; // used for default reference
   const Json::Value Phys = root.get("Physics",defValue);
   if(Phys.isObject()) {
      if(procID==0) printf("\nInitializing Physics ...\n");
      Json::Value advScheme = Phys.get("advScheme",defValue);
      if(advScheme == defValue) {
         cout << "ERROR: advScheme or gamma " << endl;
         cout << "or Nsub or delta0 is " << endl;
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

      Json::Value gammaVal  = Phys.get("gammaC",defValue);
      gamma0 = gammaVal.asDouble();
      if(gammaVal.isObject()) gamma0 = gammaVal.asDouble();
      else gamma0 = 5.0/3.0;
      if(procID==0) cout << "adiabatic coefficent = " << gamma0 << endl;
      if(gamma0 < 1.0) {
         printf("ERROR: adiabatic coefficient can't be < 1\n");
         exit (EXIT_FAILURE);
      }
      Sthresh = Nthresh*Tthresh/pow(Nthresh,gamma0-1);
      
      Json::Value NsubVal = Phys.get("Nsub",defValue);
      if(NsubVal.isObject()) Nsub = NsubVal.asInt();
      else Nsub = 2;
      if(procID==0) cout << "Nsub = " << Nsub << endl;
      if(Nsub < 1) {
         printf("ERROR: Nsub must be int >= 1\n");
         exit (EXIT_FAILURE);
      }
      
      Json::Value deltaVal  = Phys.get("delta0",defValue);
      if(deltaVal.isObject()) delta0 = deltaVal.asDouble();
      else delta0 = 1.0e-4;
      if(procID==0) cout << "relaxation constant = " << delta0 << endl;
      if(delta0 >= 1.0) {
         printf("ERROR: delta0>=1 ==> cvac<=V \n");
         exit (EXIT_FAILURE);
      }
      

   }
   else {
      cout << "value for key \"Physics\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
  

   //   get initial profile for matrix Variable
   //
   const Json::Value Wvar = Phys.get("W",defValue);
   if(Wvar.isObject()) { 
      Xgrid.setInitialProfile(matW,Wvar);
      if(procID==0) setXminBoundary(matW, 0.0, 1.0);   
      if(procID==numProcs-1) setXmaxBoundary(matW, 0.0, 1.0);   
      Xgrid.communicate(matW);
      //Wold  = matW;
   } else {
      cout << "value for Physics variable \"W\" is not object type !" << endl;
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
      S = P/pow(N,gamma0-1.0);
      Sold  = S;
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

   T = P/N;
   Cs = sqrt(gamma0*P/N);
   Xgrid.DDX(J,B); 
   Xgrid.DDX(Jcc,B); 
   Xgrid.communicate(J);
   Xgrid.communicate(Jcc);
   J0 = J;
   eta = eta0/T/sqrt(T);
   Xgrid.InterpToCellEdges(VBce,V*B,B,"C2");
   Xgrid.InterpToCellEdges(etace,eta,eta,"C2");
   Ez = etace*J-VBce; 
   Ezhalf = Ez;

   // need inital flux calculation before add to output
   // and before first call to Physics.advance   
   computeFluxes(Xgrid, 1); // 1 for first order calculation   
  

   // add stuff to output files
   //
   dataFile.add(matW, "matW", 0);  // matrixTest 
   dataFile.add(matrixA, "matrixA", 0);  // matrixTest 
   dataFile.add(matrixB, "matrixB", 0);  // matrixTest 
   dataFile.add(matrixC, "matrixC", 0);  // matrixTest 
   dataFile.add(N, "N", 1);  // density 
   dataFile.add(M, "M", 1);  // momentum density 
   dataFile.add(S, "S", 1);  // entropy density
   dataFile.add(B, "B", 1);  // magnetic field
   dataFile.add(P, "P", 1);  // pressure
   dataFile.add(T, "T", 1);  // temperature
   dataFile.add(eta, "eta", 1);  // resistivity
   dataFile.add(V, "V", 1);  // velocity
   dataFile.add(J, "J", 1);  // current density
   dataFile.add(J0, "J0", 1);  // curl of B
   dataFile.add(Ez, "Ez", 1);  // z-electric field
   dataFile.add(Cs,"Cs",1);  // sound speed
   dataFile.add(gamma0,"gamma0",0); 
   //
   dataFile.add(FluxN, "FluxN", 1);  
   dataFile.add(FluxM, "FluxM", 1);  
   dataFile.add(FluxS, "FluxS", 1);  
   dataFile.add(FluxB, "FluxB", 1);  
   dataFile.add(FluxEz, "FluxEz", 1);  
   //
   dataFile.add(FluxRatio, "FluxRatio", 1);  
   dataFile.add(FluxLim, "FluxLim", 1);  
   dataFile.add(FluxR, "FluxR", 1);
   dataFile.add(FluxL, "FluxL", 1);

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
   double thisdt, expFact;
   for (auto n=1; n<Nsub+1; n++) {
      thisdt = dt*n/Nsub;

      // Update N, M, S, and B terms from n to n+1
      for (auto i=nXg; i<nMax-nXg; i++) {
	 
	 N.at(i) = Nold.at(i) - thisdt*(FluxN.at(i)-FluxN.at(i-1))/Xgrid.dX;
         M.at(i) = Mold.at(i) - thisdt*(FluxM.at(i)-FluxM.at(i-1))/Xgrid.dX
		 - thisdt*Jcc.at(i)*Bold.at(i);
		// - thisdt*(Jcc.at(i)-J0.at(i))*Bold.at(i);
         S.at(i) = Sold.at(i) - thisdt*(FluxS.at(i)-FluxS.at(i-1))/Xgrid.dX
		 + thisdt*(gamma0-1.0)*eta.at(i)*pow(Nold.at(i),1.0-gamma0)*Jcc.at(i)*Jcc.at(i);
         
	 B.at(i) = Bold.at(i) - thisdt*(FluxB.at(i)-FluxB.at(i-1))/Xgrid.dX;
	 //B.at(i) = Bold.at(i) + thisdt*(Ezold.at(i)-Ezold.at(i-1))/Xgrid.dX;
	 
	 if(N.at(i)<Nthresh) N.at(i) = Nthresh;
	 if(M.at(i)<0.0) M.at(i) = 0.0;
	 if(S.at(i)<Sthresh) S.at(i) = Sthresh;

      }

      timeDomain* tmesh = timeDomain::tmesh;
      //cout << "thist = " << tmesh->tSim << endl;
      double thist = tmesh->tSim;
      B0 = thist*40.0;
      if(B0>4.0) B0 = 4.0;
      
      if(procID==0) {
         setXminBoundary(N, N.at(2), 0.0);   
         //setXminBoundary(M, 0.0, -1.0);   
         setXminBoundary(M, M.at(2), 0.0);   
         setXminBoundary(S, S.at(2), 0.0);
         //setXminExtrap(N);
         //setXminExtrap(M);
         //setXminExtrap(S);
         setXminBoundary(B, B0, 0.0);   
	 if(N.at(1)<Nthresh) {
            N.at(0) = Nthresh;
            N.at(1) = Nthresh;
            //cout << "Are we here?" << endl;
	 }
	 if(N.at(0)<Nthresh) {
            //cout << "Are we here?" << endl;
            N.at(0) = Nthresh;
	 }
	 if(S.at(1)<Sthresh) {
            S.at(0) = Sthresh/pow(Nthresh,gamma0-2);
            S.at(1) = Sthresh;
	 }
	 if(S.at(0)<Sthresh) {
            S.at(0) = Sthresh;
	 }
	 
	 if(M.at(1)<0.0) {
            M.at(0) = 0.0*M.at(2);
            M.at(1) = 0.0*M.at(2);
	 }
	 if(M.at(0)<0.0) {
            M.at(0) = 0.0*M.at(1);
	 }
	 
      }
     
      if(procID==numProcs-1) {
         setXmaxBoundary(N, 0.0, 1.0);   
         setXmaxBoundary(M, 0.0, -1.0);   
         setXmaxBoundary(S, 0.0, 1.0);   
         setXmaxBoundary(B, 0.0, 0.0);   
      }
      Xgrid.communicate(N);
      Xgrid.communicate(M);
      Xgrid.communicate(S);
      Xgrid.communicate(B);

      // compute fluxes using partialy updated fields at n+1/2 for 2nd order
      // or n+1 for 1st order
      computeFluxes(Xgrid, Nsub); // compute second order fluxes at n+1/2

   } // finish subcycle steps for N, M, S, and B


   // Now update electric field (stag in time wrt others)
   //
   //Xgrid.InterpToCellEdges(VBce,M/N*B,B,"C2");
   for (auto i=nXg-1; i<nMax-nXg; i++) {
      if(thisdt/delta0/etace.at(i)<=1.0e-3) {
         expFact = 1.0-thisdt/delta0/etace.at(i)+0.5*pow(thisdt/delta0/etace.at(i),2);
      }
      else {
         expFact = exp(-thisdt/delta0/etace.at(i));
      }
      Ez.at(i) = Ezold.at(i)*expFact 
    	       - ( VBce.at(i) - etace.at(i)*(B.at(i+1)-B.at(i))/Xgrid.dX 
	         )*(1.0-expFact); 
   }
   Xgrid.communicate(Ez);

   // update old fields
   Nold = N;
   Mold = M;
   Sold = S;
   Bold = B;
   Ezhalf = (Ez+Ezold)/2.0;
   Ezold = Ez;
   
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
   vector<double> Cspeed, FluxNcc, FluxMcc, FluxScc;
   vector<double> Cspeed2, Cvaceff, FluxB0cc, Eprime; 
   FluxNcc.assign(nCC,0.0);
   FluxMcc.assign(nCC,0.0);
   FluxScc.assign(nCC,0.0);
   FluxB0cc.assign(nCC,0.0);
   Eprime.assign(nCE,0.0);

   //  define derived variables
   //
   V  = M/N;
   P  = S*pow(N,gamma0-1);
   T  = P/N;
   Cs = sqrt(gamma0*P/N);
   eta = eta0/T/sqrt(T); //*(1.0+10*Nthresh/N);
   Xgrid.InterpToCellEdges(etace,eta,eta,"C2");
   Xgrid.communicate(etace);

   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeed  = abs(V + Cs); // adv flux jacobian
   Cspeed2 = abs(V - Cs); // adv flux jacobian
   const int nXcc = Cspeed.size();
   for (auto i=0; i<nXcc; i++) {
      if(Cspeed2.at(i)>Cspeed.at(i)) Cspeed.at(i) = Cspeed2.at(i);
   }
   Cvaceff.assign(nCC,1.0/sqrt(delta0));
   FluxNcc = M;
   FluxMcc = M*V + P + 0.0*B*B/2.0;
   FluxScc = V*S;
   FluxB0cc = V*B;
   FluxEz = -B;

   // compute advective flux using
   // specified scheme from input file
   //
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxN,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxNcc,Cspeed,N,Nsub);
      //Xgrid.computeFluxTVD(FluxM,FluxL,FluxR,FluxRatio,FluxLim,
      //                     FluxMcc,Cspeed,M,Nsub);
      Xgrid.computeFluxTVD(FluxS,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxScc,Cspeed,S,Nsub);
      Xgrid.computeFluxTVD(FluxM,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxMcc,Cspeed,M,Nsub);
      //Xgrid.computeFluxTVD(VBce,FluxL,FluxR,FluxRatio,FluxLim,
      //                     FluxB0cc,Cspeed,B,Nsub);
   }
   else {
      Xgrid.InterpToCellEdges(FluxN,FluxNcc,V,advScheme0);
      Xgrid.InterpToCellEdges(FluxM,FluxMcc,V,advScheme0);
      Xgrid.InterpToCellEdges(FluxS,FluxScc,V,advScheme0);
      //Xgrid.InterpToCellEdges(FluxB,FluxBcc,V,advScheme0);
      //Xgrid.InterpToCellEdges(FluxEz,FluxEzcc,Ez,advScheme0);
      //Xgrid.InterpToCellEdges(FluxEz,FluxEzcc,Ez,"C2");
   } 
   Xgrid.InterpToCellEdges(VBce,FluxB0cc,B,"C2");

   if(procID==0) {
      setXminBoundary(FluxN, 0.0, 0.0);   
      //setXminBoundary(FluxM, (P.at(2)+P.at(1))/2.0, 0.0);   
      //setXminBoundary(FluxM, (FluxMcc.at(2)+FluxMcc.at(1))/2.0, 0.0);   
      setXminBoundary(FluxS, 0.0, 0.0);
      //setXminBoundary(VBce, (FluxB0cc.at(2)+FluxB0cc.at(1))/2.0, 0.0);
   }   
   //FluxB = FluxB-Eprime;
   FluxB = -Ez;

   Xgrid.communicate(VBce);
   Eprime = Ez+VBce;
   J = (Ezhalf + VBce)/etace;
   Xgrid.communicate(J);
   Xgrid.InterpToCellCenter(Jcc,J);
   Xgrid.communicate(Jcc);
   Xgrid.DDX(J0,B); 
   Xgrid.communicate(J0);

   Xgrid.communicate(FluxL);   
   Xgrid.communicate(FluxR);   
   Xgrid.communicate(FluxRatio);   
   Xgrid.communicate(FluxLim);   

} // end computeFluxes


void setXminBoundary(vector<double>& var, const double C0, const double C1)
{
   
   domainGrid* mesh = domainGrid::mesh;
   //F0.front() = 2.0*C-F0.at(1); 
   //F0.front() = C; 
   for (auto i=0; i<mesh->nXg; i++) {
      //var.at(i) = C0;
      var.at(i) = C0 + C1*var.at(mesh->nXg+1-i);
   }
   /*
   if(C0==0) {
      var.at(1) = 2.0*var.at(2) - var.at(3);
      var.at(0) = var.at(1);
   }
   */
}

void setXminBoundary(vector<vector<double>>& var, 
		     const double C0, const double C1)
{
   const int thisnZ = var[0].size();

   domainGrid* mesh = domainGrid::mesh;
   //F0.front() = 2.0*C-F0.at(1); 
   //F0.front() = C; 
   for (auto i=0; i<mesh->nXg; i++) {
      for (auto j=0; j<thisnZ; j++) {
         //var.at(i) = C0;
         var[i][j] = C0 + C1*var[mesh->nXg+1-i][j];
      }
   }

}

void setXminExtrap(vector<double>& var)
{
   
   //var.at(1) = 2.0*var.at(2) - var.at(3);
   //var.at(0) = 2.0*var.at(1) - var.at(2);
   var.at(1) = 3.0*(var.at(2) - var.at(3)) + var.at(4);
   var.at(0) = 3.0*(var.at(1) - var.at(2)) + var.at(3);
   //var.at(0) = var.at(1);
}

void setXminExtrap(vector<vector<double>>& var)
{
   //const int thisnX = var.size();
   const int thisnZ = var[0].size();

   for (auto j=0; j<thisnZ; j++) {
      //var.at(1) = 2.0*var.at(2) - var.at(3);
      //var.at(0) = 2.0*var.at(1) - var.at(2);
      var[1][j] = 3.0*(var[2][j] - var[3][j]) + var[4][j];
      var[0][j] = 3.0*(var[1][j] - var[2][j]) + var[3][j];
      //var.at(0) = var.at(1);
   }

}


void setXmaxBoundary(vector<double>& var, const double C0, const double C1)
{
   
   const int thisnX = var.size();
   
   domainGrid* mesh = domainGrid::mesh;
   const int ishift=thisnX-mesh->nXg;
   
   //F0[nXsub+1] = 2.0*C-F0[nXsub]; 
   //F0[nXsub+1] = C;
   for (auto i=ishift; i<mesh->nXcc; i++) {
      var.at(i) = C0 + C1*var.at(2.0*ishift-i-1);
   }
      
}

void setXmaxBoundary(vector<vector<double>>& var, 
		     const double C0, const double C1)
{
   const int thisnX = var.size();
   const int thisnZ = var[0].size();

   domainGrid* mesh = domainGrid::mesh;
   const int ishift=thisnX-mesh->nXg;

   //F0[nXsub+1] = 2.0*C-F0[nXsub]; 
   //F0[nXsub+1] = C;
   for (auto i=ishift; i<thisnX; i++) {
      for (auto j=0; j<thisnZ; j++) {
         var[i][j] = C0 + C1*var[2.0*ishift-i-1][j];
      }
   }
      
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

