/***
 * 
 * physics module for 1D burgers equation
 *
 * Right now you have to manually include
 * this file in main.cpp...
 *
***/


#include "../domainGrid.h"
#include "../timeDomain.h"
#include "../variables.h"
#include "../vectorMath.h"

using namespace std;

string advScheme0;    // advection differencing scheme
double K;             // diffusion coefficient
vector<double> F0, F0old;    // function
vector<double> FluxRatio, FluxLim;
vector<double> Flux, FluxR, FluxL;  // flux at cell-edges   

string type0;    // initial function type
double a, b, c;  // initial function params
   
void computeFluxes(const domainGrid&);
void setXminBoundary(const domainGrid&, const double&);
void setXmaxBoundary(const domainGrid&, const double&);


void variables::initialize(const domainGrid& Xgrid, const Json::Value& root, 
                      HDF5dataFile& dataFile)
{
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   vector<double> Xshift;
   const int nXcc = Xgrid.Xcc.size();
   const int nXce = Xgrid.Xce.size();
   F0.assign(nXcc,0.0);
   F0old.assign(nXcc,0.0);
   FluxRatio.assign(nXce,0.0);
   FluxLim.assign(nXce,0.0);
   Flux.assign(nXce,0.0);
   FluxR.assign(nXce,0.0);
   FluxL.assign(nXce,0.0);
   const Json::Value defValue; // used for default reference
   const Json::Value Vars = root.get("Variables",defValue);
   if(Vars.isObject()) {
      if(procID==0) printf("\nInitializing Variables ...\n");
      Json::Value aVal = Vars.get("a",defValue);
      Json::Value bVal = Vars.get("b",defValue);
      Json::Value cVal = Vars.get("c",defValue);
      Json::Value type = Vars.get("type",defValue);
      Json::Value advScheme = Vars.get("advScheme",defValue);
      Json::Value KVal = Vars.get("diffC",defValue);
      if(aVal == defValue || bVal == defValue || 
         cVal == defValue || type == defValue || 
         advScheme == defValue || KVal == defValue) {
         cout << "ERROR: at least 1 Variables value is " << endl;
         cout << "not declared in input file" << endl;
         exit (EXIT_FAILURE);
      } 
      a = aVal.asDouble();
      b = bVal.asDouble();
      c = cVal.asDouble();
      
      type0 = type.asString();
      if(type0=="gaussian" || type0=="Gaussian") {
         
         Xshift = (Xgrid.Xcc-b)/c;
         F0 = a*exp(-Xshift*Xshift/2.0);
         
         if(c < 0.0) {
            cout << "ERROR: gaussian width c is not positive" << endl; 
            cout << "value in input file" << endl;
            exit (EXIT_FAILURE);
         }
         if(procID==0) {
            cout << "Initial F0 is Gaussian with amplitude = " << a << endl;
            cout << "center at x = " << b << endl;
            cout << "and width = " << c << endl;
         }

      }
      else if(type0=="tanh") {
         
         Xshift = (Xgrid.Xcc-b)/c;
         F0 = a*tanh(Xshift);
         
         if(c < 0.0) {
            cout << "ERROR: tanh width c is not positive" << endl; 
            cout << "value in input file" << endl;
            exit (EXIT_FAILURE);
         }
         if(procID==0) {
            cout << "Initial F0 is tanh with amplitude = " << a << endl;
            cout << "center at x = " << b << endl;
            cout << "and width = " << c << endl;
         }
      }
      else {
         cout << "Initial Variable type = " << type0 << " is not valid " << endl;
         exit (EXIT_FAILURE);
      }

      advScheme0 = advScheme.asString();
      if(advScheme0=="C2" || advScheme0=="U1" || 
         advScheme0=="QUICK" || advScheme0=="TVD") {
         if(procID==0) {
            cout << "advection diff/interp scheme is " << advScheme0 << endl;
            //cout << endl;
         }
      }
      else {
         cout << "advection scheme " << advScheme0 << " is not valid " << endl;
         cout << "valid types are C2, U1, QUICK, and TVD " << endl;
         exit (EXIT_FAILURE);
      }

      K = KVal.asDouble();
      if(procID==0) cout << "diffusion coefficent = " << K << endl;
      if(K < 0.0) {
         printf("ERROR: diffusion coefficient can't be < 0\n");
         exit (EXIT_FAILURE);
      }

   }
   else {
      cout << "value for key \"Variables\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
  
   cout << endl;  


   if(procID==0) setXminBoundary(Xgrid, 0.0);   
   if(procID==numProcs-1) setXmaxBoundary(Xgrid, 0.0);   
   Xgrid.communicate(F0);
   F0old  = F0;
   computeFluxes(Xgrid); // inital calculation before add to output   
   dataFile.add(F0, "F0", 1); // function   
   dataFile.add(FluxRatio, "FluxRatio", 1); // function   
   dataFile.add(FluxLim, "FluxLim", 1); // function   
   dataFile.add(Flux, "Flux", 1); // total flux function   
   dataFile.add(FluxR, "FluxR", 1); // right going flux
   dataFile.add(FluxL, "FluxL", 1); // left going flux

}


void computeFluxes(const domainGrid& Xgrid)
{
   // compute Flux = F0^2/2 - K*dF0/dX
   //              = FluxAdv + FluxDif
   //
   // FluxAdv = (FluxR+FluxL)/2 is computed using upwinding schemes
   // FluxDif is computed using standard centered scheme


   const int nCE = Flux.size();
   const int nCC = F0.size();
   vector<double> Cspeed, FluxAdvCC, FluxAdv, FluxDif;
   FluxAdvCC.assign(nCC,0.0);
   FluxAdv.assign(nCE,0.0);
   FluxDif.assign(nCE,0.0);
   

   // set flux freezing speed and 
   // compute advection flux at cell center
   //
   Cspeed = F0; // adv flux jacobian
   FluxAdvCC = F0*F0*0.5;
   /*
   for (auto i=0; i<nCC; i++) {
      FluxAdvCC.at(i) = F0.at(i)*F0.at(i)/2.0;
   }
   */


   // compute diffusive flux using 
   // standard centered scheme
   //
   Xgrid.DDX(FluxDif,F0);
   //FluxDif = DDX(F0,Xgrid.dX);
   transform(FluxDif.begin(), FluxDif.end(), FluxDif.begin(), 
             bind1st(multiplies<double>(),-K)); 


   // compute advective flux using
   // specified scheme from input file
   //
   if(advScheme0 == "TVD") {
      Xgrid.computeFluxTVD(FluxAdv,FluxL,FluxR,FluxRatio,FluxLim,
                           FluxAdvCC,Cspeed,F0);
   }      
   else {
      Xgrid.InterpToCellEdges(FluxAdv,FluxAdvCC,Cspeed,advScheme0);
   } 


   // add advective and diffusive flux together
   //
   Flux = FluxAdv + FluxDif;
   //transform(Flux.begin(), Flux.end(), FluxDif.begin(), 
   //               Flux.begin(), plus<double>());


}


void setXminBoundary(const domainGrid& Xgrid, const double& C)
{
   
   //F0.front() = 2.0*C-F0.at(1); 
   F0.front() = C; 
   
}


void setXmaxBoundary(const domainGrid& Xgrid, const double& C)
{
   
   //F0[nXsub+1] = 2.0*C-F0[nXsub]; 
   //F0[nXsub+1] = C; 
   F0.back() = C; 
      

}


void variables::advanceF0(const domainGrid& Xgrid, const double& dt)
{
   const int nMax = F0.size();
   const int nXg = Xgrid.nXg;
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   // Explicit forward advance from n to n+1/2 using Flux(n)
   // (calc of Flux(t=0) is done during initilization)
   //
   for (auto n=nXg; n<nMax-nXg; n++) {
      F0.at(n) = F0old.at(n) - dt/2.0*(Flux.at(n)-Flux.at(n-1))/Xgrid.dX;
   }

   // apply boundary conditions and communicate
   //
   if(procID==0) setXminBoundary(Xgrid, 0.0);   
   if(procID==numProcs-1) setXmaxBoundary(Xgrid, 0.0);   
   Xgrid.communicate(F0);

   // compute RHS using F0(n+1/2)
   //
   computeFluxes(Xgrid);

   // Explicit forward advance from n to n+1 using Flux(n+1/2)
   //
   for (auto n=nXg; n<nMax-nXg; n++) {
      F0.at(n) = F0old.at(n) - dt*(Flux.at(n)-Flux.at(n-1))/Xgrid.dX;
   }
   
   // apply boundary conditions and communicate
   //
   if(procID==0) setXminBoundary(Xgrid, 0.0);   
   if(procID==numProcs-1) setXmaxBoundary(Xgrid, 0.0);   
   Xgrid.communicate(F0);
   computeFluxes(Xgrid);
   
   //  update F0old
   //
   F0old = F0;

}

void variables::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   
   const double dX = Xgrid.dX;
   double dtmaxDif = 0.5*dX*dX/K;
   double dtmaxAdv = dX;
   double dtmax = min(dtmaxDif,dtmaxAdv);
   dtSim = min(dtmax/tDom.dtFrac,tDom.dtOut);
   if(procID==0) {
      cout << endl; 
      cout << "max stable time step is " << dtmax << endl;
      cout << endl; 
   }
}

