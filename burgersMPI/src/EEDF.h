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


using namespace std;

class EEDF
{
public:
   double a, b, c;       // zero moment
   string type0;         // initial type of EEDF
   string advScheme0;    // advection differencing scheme
   vector<double> F0, F0old, F0half; // EEDF
   vector<double> Flux;  // Energy space adv, diff, and flux at cell-edge
   
   void initialize(const domainGrid&, const Json::Value&);
   void computeFluxes(const domainGrid&, const double&);
   void communicate(const domainGrid&);
   void advanceF0(const domainGrid&, const double&);
   void setXminBoundary(const domainGrid&, const double&);
   void setXmaxBoundary(const domainGrid&, const double&);

private:
   //double Te;
};


void EEDF::initialize(const domainGrid& Xgrid, const Json::Value& root)
{
   int procID;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);

   double Xshift;
   const int nXcc = Xgrid.Xcc.size();
   const int nXce = Xgrid.Xce.size();
   F0.assign(nXcc,0.0);
   Flux.assign(nXce,0.0);
   const Json::Value defValue; // used for default reference
   const Json::Value EEDF = root.get("EEDF",defValue);
   if(EEDF.isObject()) {
      if(procID==0) printf("\nInitializing EEDF ...\n");
      Json::Value aVal = EEDF.get("a",defValue);
      Json::Value bVal = EEDF.get("b",defValue);
      Json::Value cVal = EEDF.get("c",defValue);
      Json::Value type = EEDF.get("type",defValue);
      Json::Value advScheme = EEDF.get("advScheme",defValue);
      if(aVal == defValue || bVal == defValue || 
         cVal == defValue || type == defValue || 
         advScheme == defValue) {
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
      if(advScheme0=="C2" || advScheme0=="U1" || advScheme0=="QUICK") {
         if(procID==0) {
            cout << "advection diff/interp scheme is " << advScheme0 << endl;
            //cout << endl;
         }
      }
      else {
         cout << "advection scheme " << advScheme0 << " is not valid " << endl;
         cout << "valid types are C2, U1, and QUICK " << endl;
         exit (EXIT_FAILURE);
      }
   }
   else {
      cout << "value for key \"Xgrid\" is not object type !" << endl;
      exit (EXIT_FAILURE);
   }
   F0old = F0;
   F0half = F0;
  
   cout << endl;  
}


void EEDF::computeFluxes(const domainGrid& Xgrid, const double& Kappa)
{
   //const vector<double> Xcempi = Xgrid.Xcempi;
   const double dX = Xgrid.dX;
   const int iMax = Flux.size();


   // compute Flux = F0^2/2 - K*dF0/dX
   //

   if(advScheme0 == "C2") {
      // Use second order central differencing/interpolation
      //
      for (auto i=0; i<iMax; i++) {
         Flux.at(i) = (F0.at(i+1)+F0.at(i))/2.0*(F0.at(i+1)+F0.at(i))/2.0/2.0 
                    - Kappa*(F0.at(i+1)-F0.at(i))/dX;
      }
   } 
   else if(advScheme0 == "U1") {
      // Use first order upwinding
      //
      double Ui, ap, am;
      for (auto i=0; i<iMax; i++) {
      
         Ui = F0.at(i)/2.0;
         ap = 1.0;
         am = 0.0;
         if(Ui<0.0) {
            ap = 0.0;
            am = 1.0; 
         }

         Flux.at(i) = ap*F0.at(i)*F0.at(i)/2.0 
                    + am*F0.at(i+1)*F0.at(i+1)/2.0 
                    - Kappa*(F0.at(i+1)-F0.at(i))/dX;
   
      } // end for loop
   }
   else if(advScheme0 == "QUICK") {   
      // Use 2nd order QUICK upwinding
      //
      double Ui, ap, am;
      double a0 = 3.0/4.0, a1 = 3.0/8.0, a2 = 1.0/8.0;
      for (auto i=0; i<iMax; i++) {
      
         Ui = F0.at(i)/2.0;
         ap = 1.0;
         am = 0.0;
         if(Ui<0.0) {
            ap = 0.0;
            am = 1.0; 
         }

         if(i==0 || i==iMax-1) {
            Flux.at(i) = ap*F0.at(i)*F0.at(i)/2.0 
                       + am*F0.at(i+1)*F0.at(i+1)/2.0 
                       - Kappa*(F0.at(i+1)-F0.at(i))/dX;
   
         } else {
            Flux.at(i) = ap*( a0*F0.at(i)*F0.at(i) 
                           +  a1*F0.at(i+1)*F0.at(i+1) 
                           -  a2*F0.at(i-1)*F0.at(i-1) )/2.0 
                       + am*( a0*F0.at(i+1)*F0.at(i+1) 
                           +  a1*F0.at(i)*F0.at(i)
                           -  a2*F0.at(i+2)*F0.at(i+2) )/2.0 
                       - Kappa*(F0.at(i+1)-F0.at(i))/dX;
         }

      } // end for loop
   } else {
      cout << "advection scheme not valid," << endl;
      cout << "should be caught on initilization!" << endl;
      exit (EXIT_FAILURE);
   }

}


void EEDF::communicate(const domainGrid& Xgrid)
{
   const int nMax = F0.size(); // number of cell-center points
   const int nXg = Xgrid.nXg; // number of guard cells each end
   int procID, numProcs;
   double Fsend[nXg], Frecv[nXg];
   
   MPI_Status status;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   //  send F0[nXg+i], i<nXg for ID>0 to proc ID-1
   //
   if (procID>0) {
      int tag = 1;
      for (auto i=0; i<nXg; i++) {
         Fsend[i] = F0[nXg+i];
      }
      MPI_Send(Fsend, nXg, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD);
      //MPI_Send(&F0.at(nXg), 1, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD);
   }  
  
   // receive: F0[nXg+i,ID+1] => F0[nMax-nXg+i,ID]
   //
   if (procID<numProcs-1) {
      int tag = 1;
      MPI_Recv(Frecv, nXg, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD, &status);
      //MPI_Recv(&F0[nXsub+1], 1, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD, &status);
      for (auto i=0; i<nXg; i++) {
         F0.at(nMax-nXg+i) = Frecv[i];
      }
   }

   // send: F0[nMax-2*nXg+i,ID] => ID+1
   //
   if (procID<numProcs-1) {
      int tag = 2;
      for (auto i=0; i<nXg; i++) {
         Fsend[i] = F0.at(nMax-2*nXg+i);
      }
      MPI_Send(Fsend, nXg, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD);
      //MPI_Send(&F0[nXsub], 1, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD);
   }  

   // receive: F0[nMax-2*nXg+i,ID] => F0[i,ID+1]
   //
   if (procID>0) {
      int tag = 2;
      MPI_Recv(Frecv, nXg, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD, &status);
      //MPI_Recv(&F0[0], 1, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD, &status);
      for (auto i=0; i<nXg; i++) {
         F0.at(i) = Frecv[i];
      }

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
      F0.at(n) = F0.at(n) - dt*(Flux.at(n)-Flux.at(n-1))/Xgrid.dX;
   }

}


#endif
