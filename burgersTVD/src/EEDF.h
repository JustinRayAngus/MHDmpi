/***
 * 
 * electron-energy-distribution function class
 *
***/

#ifndef EEDF_h
#define EEDF_h
//#ifndef __EEDF_H_INCLUDED__
//#define __EEDF_H_INCLUDED__


#include "domainGrid.h"
#include "timeDomain.h"

using namespace std;

class EEDF
{
public:
   string advScheme0;    // advection differencing scheme
   double K;             // diffusion coefficient
   vector<double> F0, F0old;    // function
   vector<double> FluxRatio, FluxLim;
   vector<double> Flux, FluxR, FluxL;  // flux at cell-edges
   
   void initialize(const domainGrid&, const Json::Value&);
   void computeFluxes(const domainGrid&);
   void advanceF0(const domainGrid&, const double&);
   void setXminBoundary(const domainGrid&, const double&);
   void setXmaxBoundary(const domainGrid&, const double&);
   void setdtSim(double&, const timeDomain&, const domainGrid&);

private:
   string type0;    // initial function type
   double a, b, c;  // initial function params
};


void EEDF::initialize(const domainGrid& Xgrid, const Json::Value& root)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);

   double Xshift;
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
   const Json::Value EEDF = root.get("EEDF",defValue);
   if(EEDF.isObject()) {
      if(procID==0) printf("\nInitializing EEDF ...\n");
      Json::Value aVal = EEDF.get("a",defValue);
      Json::Value bVal = EEDF.get("b",defValue);
      Json::Value cVal = EEDF.get("c",defValue);
      Json::Value type = EEDF.get("type",defValue);
      Json::Value advScheme = EEDF.get("advScheme",defValue);
      Json::Value KVal = EEDF.get("diffC",defValue);
      if(aVal == defValue || bVal == defValue || 
         cVal == defValue || type == defValue || 
         advScheme == defValue || KVal == defValue) {
         cout << "ERROR: at least 1 EEDF value is " << endl;
         cout << "not declared in input file" << endl;
         exit (EXIT_FAILURE);
      } 
      a = aVal.asDouble();
      b = bVal.asDouble();
      c = cVal.asDouble();
      if(c < 0.0) {
         printf("ERROR: gaussian width c is not a positive value in input file\n");
         exit (EXIT_FAILURE);
      }
      
      type0 = type.asString();
      if(type0=="gaussian" || type0=="Gaussian") {
         
         for (auto n=0; n<Xgrid.nXsub+2; n++) {
            //F0[n] = 1.0-0.25*Xgrid.Xcc[n]*Xgrid.Xcc[n];
            Xshift = (Xgrid.Xcc[n]-b)/c;
            F0.at(n) = a*exp(-Xshift*Xshift/2.0);
         }
         //transform(F0.begin(), F0.end(), F0.begin(), 
         //bind1st(multiplies<double>(),1.0/zeroMom)); 
         if(procID==0) {
            cout << "Initial F0 is Gaussian with amplitude = " << a << endl;
            cout << "center at x = " << b << endl;
            cout << "and width = " << c << endl;
         }

      }
      else {
         cout << "Initial EEDF type = " << type0 << " is not valid " << endl;
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
      cout << "value for key \"Xgrid\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
  
   cout << endl;  
}


void EEDF::computeFluxes(const domainGrid& Xgrid)
{
   const double dX = Xgrid.dX;
   const int nCE = Flux.size();
   const vector<double> U=F0;
   const int nCC = U.size();
 

   // compute Flux = U^2/2 - K*dFU/dX
   //              = FluxAdv + FluxDif
   //
   // FluxAdv = (FluxR+FluxL)/2 is computed using upwinding schemes
   // FluxDif is computed using standard centered scheme


   // Step 1: set flux freezing speed and 
   // compute advection flux at cell center
   
   vector<double> Cspeed, FluxAdvCC, FluxAdv, FluxDif;
   vector<double> FluxR1st, FluxL1st, DeltaFluxR, DeltaFluxL;
   FluxAdvCC.assign(nCC,0.0);
   FluxAdv.assign(nCE,0.0);
   DeltaFluxL.assign(nCE,0.0);
   DeltaFluxR.assign(nCE,0.0);
   FluxL1st.assign(nCE,0.0);
   FluxR1st.assign(nCE,0.0);
   Cspeed = U; // Local flux freezing speed 
   for (auto i=0; i<nCC; i++) {
      FluxAdvCC.at(i) = U.at(i)*U.at(i)/2.0;
   }


   // compute diffusive flux using standard centered scheme
   //
   FluxDif.assign(nCE,0.0);
   //domainGrid domGrid;
   Xgrid.DDX(FluxDif,U);
   transform(FluxDif.begin(), FluxDif.end(), FluxDif.begin(), 
             bind1st(multiplies<double>(),-K)); 


   // Step 2: compute first order upwind fluxes
   // for FluxL and FluxR at cell edges
   

   if(advScheme0 == "TVD") {
      
      // Step 2: compute first order upwind fluxes
      // for FluxL and FluxR at cell edges
      //
      for (auto i=0; i<nCE; i++) {
         FluxR1st.at(i) = FluxAdvCC.at(i)   + Cspeed.at(i)*U.at(i);
         FluxL1st.at(i) = FluxAdvCC.at(i+1) - Cspeed.at(i+1)*U.at(i+1);
         //FluxDif.at(i) = -K*(U.at(i+1)-U.at(i))/dX;
      }

      // Step 3: compute 2nd order corrections to FluxR and FluxL
      //
      for (auto i=1; i<nCE-1; i++) {
         DeltaFluxR.at(i) = 0.5*(FluxR1st.at(i+1)-FluxR1st.at(i)); 
         DeltaFluxL.at(i) = 0.5*(FluxL1st.at(i)-FluxL1st.at(i-1)); 
         //FluxRatio.at(i) = DeltaFluxL.at(i)/DeltaFluxR.at(i);
         FluxRatio.at(i) = (U.at(i+1) - U.at(i-1))/(U.at(i+2) - U.at(i)); // divide by zero?
         //FluxLim.at(i) = (abs(FluxRatio.at(i))+FluxRatio.at(i))/(abs(FluxRatio.at(i)) + 1.0); // van Leer
         FluxLim.at(i) = 2.0;
         if(FluxRatio.at(i)<=2.0) FluxLim.at(i) = FluxRatio.at(i);
         if(FluxRatio.at(i)<=1.0) FluxLim.at(i) = 1.0;
         if(FluxRatio.at(i)<=0.5) FluxLim.at(i) = 2.0*FluxRatio.at(i);
         if(FluxRatio.at(i)<=0.0) FluxLim.at(i) = 0.0;
         
         FluxR.at(i) = FluxR1st.at(i) + FluxLim.at(i)*DeltaFluxR.at(i);
         FluxL.at(i) = FluxL1st.at(i) + FluxLim.at(i)*DeltaFluxL.at(i);
         FluxAdv.at(i) = (FluxR.at(i) + FluxL.at(i))/2.0;
         Flux.at(i) = FluxAdv.at(i) + FluxDif.at(i);
      }
      //FluxAdv.at(0) = (FluxR.at(0) + FluxL.at(0))/2.0;
      //Flux.at(0) = FluxAdv.at(0) + FluxDif.at(0);
      //FluxAdv.at(nCE-1) = (FluxR.at(nCE-1) + FluxL.at(nCE-1))/2.0;
      //Flux.at(nCE-1) = FluxAdv.at(nCE-1) + FluxDif.at(nCE-1);
      

   }      
   else if(advScheme0 == "C2") {
      // Use second order central differencing/interpolation
      //
      for (auto i=0; i<nCE; i++) {
         Flux.at(i) = (U.at(i+1)+U.at(i))/2.0*(U.at(i+1)+U.at(i))/2.0/2.0 
                    - K*(U.at(i+1)-U.at(i))/dX;
      }
   } 
   else if(advScheme0 == "U1") {
      // Use first order upwinding
      //
      double Ui, ap, am;
      for (auto i=0; i<nCE; i++) {
      
         Ui = U.at(i)/2.0;
         ap = 1.0;
         am = 0.0;
         if(Ui<0.0) {
            ap = 0.0;
            am = 1.0; 
         }

         Flux.at(i) = ap*U.at(i)*U.at(i)/2.0 
                    + am*U.at(i+1)*U.at(i+1)/2.0 
                    - K*(U.at(i+1)-U.at(i))/dX;
   
      } // end for loop
   }
   else if(advScheme0 == "QUICK") {   
      // Use 2nd order QUICK upwinding
      //
      double Ui, ap, am;
      double a0 = 3.0/4.0, a1 = 3.0/8.0, a2 = 1.0/8.0;
      for (auto i=0; i<nCE; i++) {
      
         Ui = U.at(i)/2.0;
         ap = 1.0;
         am = 0.0;
         if(Ui<0.0) {
            ap = 0.0;
            am = 1.0; 
         }

         if(i==0 || i==nCE-1) {
            Flux.at(i) = ap*U.at(i)*U.at(i)/2.0 
                       + am*U.at(i+1)*U.at(i+1)/2.0 
                       - K*(U.at(i+1)-U.at(i))/dX;
   
         } else {
            Flux.at(i) = ap*( a0*U.at(i)*U.at(i) 
                           +  a1*U.at(i+1)*U.at(i+1) 
                           -  a2*U.at(i-1)*U.at(i-1) )/2.0 
                       + am*( a0*U.at(i+1)*U.at(i+1) 
                           +  a1*U.at(i)*U.at(i)
                           -  a2*U.at(i+2)*U.at(i+2) )/2.0 
                       - K*(U.at(i+1)-U.at(i))/dX;
         }

      } // end for loop
   } else {
      cout << "advection scheme not valid," << endl;
      cout << "should be caught on initilization!" << endl;
      exit (EXIT_FAILURE);
   }

}


void EEDF::setXminBoundary(const domainGrid& Xgrid, const double& C)
{
   
   //F0.front() = 2.0*C-F0.at(1); 
   F0.front() = C; 
   
}


void EEDF::setXmaxBoundary(const domainGrid& Xgrid, const double& C)
{
   
   //F0[nXsub+1] = 2.0*C-F0[nXsub]; 
   //F0[nXsub+1] = C; 
   F0.back() = C; 
      

}


void EEDF::advanceF0(const domainGrid& Xgrid, const double& dt)
{
   //const int nXsub = Xgrid.nXsub;
   const int nMax = F0.size();
   const int nXg = Xgrid.nXg;
   
   // Explicit forward advance
   //
   for (auto n=nXg; n<nMax-nXg; n++) {
      F0.at(n) = F0old.at(n) - dt*(Flux.at(n)-Flux.at(n-1))/Xgrid.dX;
   }

}

void EEDF::setdtSim(double& dtSim, const timeDomain& tDom, const domainGrid& Xgrid)
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


#endif
