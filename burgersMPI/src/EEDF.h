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
   double Xshift;
   F0.assign(Xgrid.nXsub+2,0.0);
   Flux.assign(Xgrid.nXsub+1,0.0);
   const Json::Value defValue; // used for default reference
   const Json::Value EEDF = root.get("EEDF",defValue);
   if(EEDF.isObject()) {
      printf("Initializing EEDF ...\n");
      Json::Value aVal = EEDF.get("a",defValue);
      Json::Value bVal = EEDF.get("b",defValue);
      Json::Value cVal = EEDF.get("c",defValue);
      Json::Value type = EEDF.get("type",defValue);
      if(aVal == defValue || bVal == defValue || 
         cVal == defValue || type == defValue) {
         cout << "ERROR: amplitude, center, width," << endl;
         cout << "or type not declared in input file" << endl;
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
         cout << "Initial F0 is Gaussian with amplitude = " << a << endl;
         cout << "center at x = " << b << endl;
         cout << "and width = " << c << endl;
         for (auto n=0; n<Xgrid.nXsub+2; n++) {
            //F0[n] = 1.0-0.25*Xgrid.Xcc[n]*Xgrid.Xcc[n];
            Xshift = (Xgrid.Xcc[n]-b)/c;
            F0.at(n) = a*exp(-Xshift*Xshift/2.0);
         }
         //transform(F0.begin(), F0.end(), F0.begin(), 
         //bind1st(multiplies<double>(),1.0/zeroMom)); 
      }
      else {
         cout << "Initial EEDF type = " << type0 << " is not valid " << endl;
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

   // Use central differencing
   //
   //for (auto i=0; i<iMax; i++) {
   //   Flux.at(i) = (F0.at(i+1)+F0.at(i))*(F0.at(i+1)+F0.at(i))/2.0 
   //              - Kappa*(F0.at(i+1)-F0.at(i))/dX;
   //}
   

   // Use first order upwinding
   //
   double Ui, ap, am;
   for (auto i=0; i<iMax; i++) {
      
      //if(i==0) { cout << "F0[0] =" << F0.at(0) << endl; }
      //if(i==0) { cout << "F0[-1] =" << F0.at(-1) << endl; }
      //if(i==0) { cout << "F0[iMax+2] =" << F0.at(iMax+2) << endl; }
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
   
   }

}


void EEDF::communicate(const domainGrid& Xgrid)
{
   const int nXsub = Xgrid.nXsub; // number of cell-center points
   int procID, numProcs;
   
   MPI_Status status;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   //  send F0[1] for ID>0 to proc ID-1
   //
   if (procID>0) {
      int tag = 1;
      MPI_Send(&F0[1], 1, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD);
   }  
  
   // receive: F0[1,ID+1] => F0[nXsub+1,ID]
   //
   if (procID<numProcs-1) {
      int tag = 1;
      MPI_Recv(&F0[nXsub+1], 1, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD, &status);
   }

   // send: F0[nXsub,ID] => ID+1
   //
   if (procID<numProcs-1) {
      int tag = 2;
      MPI_Send(&F0[nXsub], 1, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD);
   }  

   // receive: F0[nXsub,ID] => F0[0,ID+1]
   //
   if (procID>0) {
      int tag = 2;
      MPI_Recv(&F0[0], 1, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD, &status);
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
   
   // Explicit forward advance
   //
   for (auto n=1; n<nMax-1; n++) {
      F0.at(n) = F0.at(n) - dt*(Flux.at(n)-Flux.at(n-1))/Xgrid.dX;
   }

}


#endif
