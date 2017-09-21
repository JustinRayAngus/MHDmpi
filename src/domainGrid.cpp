/***
 *
 * domainGrid class source file
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

using namespace std;


void domainGrid::initialize(const Json::Value& root)
{
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   const Json::Value defValue; // used for default reference
   const Json::Value Xgrid = root.get("Xgrid",defValue);
   if(Xgrid.isObject()) {
      if(procID==0) printf("\nInitializing domain grid ...\n");
      Json::Value XminVal = Xgrid.get("Xmin",defValue);
      Json::Value XmaxVal = Xgrid.get("Xmax",defValue);
      Json::Value nXVal   = Xgrid.get("nX",defValue);
      Json::Value nXgVal  = Xgrid.get("nXg",defValue);
      if(XminVal == defValue || XmaxVal == defValue || 
           nXVal == defValue || nXgVal == defValue) {
         printf("ERROR: Xmin, Xmax, nX, or nXg not declared in input file\n");
         exit (EXIT_FAILURE);
      } 
      nX = nXVal.asInt();
      if(nX != nXVal.asDouble() || nX < 1) {
         printf("ERROR: nX is not set as a positive integer in input file\n");
         exit (EXIT_FAILURE);
      }
      nXg = nXgVal.asInt();
      if(nXg != nXgVal.asDouble() || nXg < 1 || nXg>4) {
         cout << "ERROR: number of guard cells nXg" << endl;
         cout << "is not set correctly in input file" << endl;
         exit (EXIT_FAILURE);
      }
      Xmin = XminVal.asDouble();
      Xmax = XmaxVal.asDouble();
      if(Xmin >= Xmax) {
         cout << "ERROR: Xmin > Xmax in input file" << endl;
         exit (EXIT_FAILURE);
      }
      dX = (Xmax-Xmin)/(double)nX;
      nXsub = nX/numProcs;
      if(procID==0) {
         double nXsubTest = nX/(double)numProcs;
         if(floor(nXsubTest)==ceil(nXsubTest)) {
            //cout << "nXsub = " << nXsub << endl;
         }
         else {
            printf("ERROR: nXsub=nX/numProcs is not an integer!!!!\n");
            exit (EXIT_FAILURE); 
         }
      }
 
      if(procID==0) {
         cout << "Xmin = " << Xmin << endl;
         cout << "Xmax = " << Xmax << endl;
         cout << "nX = " << nX << endl;
         cout << "nXsub = " << nXsub << endl;
         cout << "dX = " << dX << endl;
         cout << "nXg = " << nXg << endl;
         cout << endl;
      }

   }
   else {
      cout << "value for key \"Xgrid\" is not object type !" << endl;
   }
  
   Xcc.assign(nXsub+2*nXg,0.0);
   nXcc = Xcc.size();
   Xce.assign(nXsub+2*nXg-1,0.0);
   nXce = Xce.size();
   double offset = procID*(Xmax-Xmin)/numProcs-(0.5+nXg-1.0)*dX;
   const int nMax = Xce.size();
   for (auto n=0; n<nMax; n++) {
      Xcc.at(n) = Xmin + offset + n*dX;
      Xce.at(n) = Xmin + offset + n*dX + 0.5*dX;
   }
   Xcc.at(nMax) = Xcc.at(nMax-1)+dX;
   //Xcc[nX/numProcs+1] = Xmin + offset + (nXsub+1)*dX;

}

void domainGrid::setInitialProfile(vector<double> &var, 
                   const Json::Value &varRoot) const
{
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   double a, b, c, d;
   vector<double> Xshift;
   string type0;

   const int Nvar = var.size();
   assert(Nvar == nXcc);

   const Json::Value defValue; // used for default reference
   Json::Value aVal = varRoot.get("a",defValue);
   Json::Value bVal = varRoot.get("b",defValue);
   Json::Value cVal = varRoot.get("c",defValue);
   Json::Value dVal = varRoot.get("d",defValue);
   Json::Value type = varRoot.get("type",defValue);
   if(aVal == defValue || bVal == defValue || 
      cVal == defValue || dVal == defValue || 
      type == defValue) {
      cout << "ERROR: at least 1 initial profile value " << endl;
      cout << "not declared in input file" << endl;
      exit (EXIT_FAILURE);
   } 
   a = aVal.asDouble();
   b = bVal.asDouble();
   c = cVal.asDouble();
   d = dVal.asDouble();
      
   type0 = type.asString();
   if(type0=="gaussian" || type0=="Gaussian") {
      
      Xshift = (Xcc-b)/c;
      var = a*exp(-Xshift*Xshift/2.0) + d;
      
      if(c == 0.0) {
         cout << "ERROR: gaussian width c cannot be zero" << endl; 
         exit (EXIT_FAILURE);
      }
      if(procID==0) {
         cout << "Initial var is Gaussian with amplitude = " << a << endl;
         cout << "center at x = " << b << endl;
         cout << "width = " << c << endl;
         cout << "and y-shift = " << d << endl;
      }

   }
   else if(type0=="tanh") {
      
      Xshift = (Xcc-b)/c;
      var = a*tanh(Xshift) + d;
      
      if(c == 0.0) {
         cout << "ERROR: tanh width c cannot be zero" << endl; 
         exit (EXIT_FAILURE);
      }
      if(procID==0) {
         cout << "Initial F0 is tanh with amplitude = " << a << endl;
      }

   } else {
     var = 1.0 + 0.0*Xcc;
   }

} // end setInitialProfile

void domainGrid::communicate(vector<double> &F0) const {

   const int nMax = F0.size(); // number of cell-center points
   const int nXce = Xce.size(); // number of cell-center points
   int thisnXg = nXg;
   if(nMax == nXce) thisnXg = nXg-1; // receiving for fluxes 

   int procID, numProcs;
   double Fsend[nXg], Frecv[nXg];

   //cout << "nMax = " << nMax << endl;

   MPI_Status status;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   //  send F0[nXg+i], i<nXg for ID>0 to proc ID-1
   //
   if (procID>0) {
      int tag = 1;
      for (auto i=0; i<thisnXg; i++) {
         Fsend[i] = F0.at(nXg+i);
      }
      MPI_Send(Fsend, thisnXg, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD);
      //MPI_Send(&F0.at(nXg), 1, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD);
   }  
  
   // receive: F0[nXg+i,ID+1] => F0[nMax-nXg+i,ID]
   //
   if (procID<numProcs-1) {
      int tag = 1;
      MPI_Recv(Frecv, thisnXg, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD, &status);
      //MPI_Recv(&F0[nXsub+1], 1, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD, &status);
      for (auto i=0; i<thisnXg; i++) {
         F0.at(nMax-thisnXg+i) = Frecv[i];
      }
   }

   // send: F0[nMax-2*nXg+i,ID] => ID+1
   //
   if (procID<numProcs-1) {
      int tag = 2;
      for (auto i=0; i<thisnXg; i++) {
         Fsend[i] = F0.at(nMax-2*nXg+i);
         if(nMax==nXce) Fsend[i] = F0.at(nMax+1-2*nXg+i);
      }
      MPI_Send(Fsend, thisnXg, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD);
      //MPI_Send(&F0[nXsub], 1, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD);
   }  

   // receive: F0[nMax-2*nXg+i,ID] => F0[i,ID+1]
   //
   if (procID>0) {
      int tag = 2;
      MPI_Recv(Frecv, thisnXg, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD, &status);
      //MPI_Recv(&F0[0], 1, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD, &status);
      for (auto i=0; i<thisnXg; i++) {
         F0.at(i) = Frecv[i];
      }

   }

}

void domainGrid::DDX(vector<double> &Fout, const vector<double> &Fin) const {

   // check that Fin (at cell center) and Fout (cell edges) are proper size
   
   const int Nout = Fout.size();
   const int Nin  = Fin.size();
   assert(Nout == nXce);
   assert(Nin == nXcc);

   for (auto i=0; i<Nout; i++) {
      Fout.at(i) = (Fin.at(i+1)-Fin.at(i))/dX;
   }

}

void domainGrid::InterpToCellEdges(vector<double> &Fout, 
                                   const vector<double> &Fin,
                                   const vector<double> &upC,
                                   const string& METHOD) const {

   // this function interpolates cell-center Fin to cell
   // to Fout, defined at cell edges
   // upC is the local maximum characteristic speed
   // METHOD referes to interpolation scheme


   // check that vectors in call are proper size
   //
   const int Nout = Fout.size();
   const int Nin  = Fin.size();
   const int NupC = upC.size();
   assert(Nout == nXce);
   assert(Nin  == nXcc);
   assert(NupC == nXcc);


   //  interpolate using specifed method
   //
   if(METHOD == "C2") { // 2nd order central

      for (auto i=0; i<Nout; i++) {
         Fout.at(i) = (Fin.at(i+1)+Fin.at(i))/2.0;
      }

   } // end METHOD=C2

   else if(METHOD == "U1") { // first order upwind
      
      for (auto i=0; i<Nout; i++) {
      
         if(upC.at(i)<0.0) {
            Fout.at(i) = Fin.at(i+1);
         } else {
            Fout.at(i) = Fin.at(i);
         }
   
      }
   
   } // end METHOD=U1

   else if(METHOD == "QUICK") {   
      // Use 2nd order QUICK upwinding
      //
      double Ui, ap, am;
      double a0 = 3.0/4.0, a1 = 3.0/8.0, a2 = 1.0/8.0;
      for (auto i=0; i<Nout; i++) {
      
         Ui = upC.at(i);
         ap = 1.0;
         am = 0.0;
         if(Ui<0.0) {
            ap = 0.0;
            am = 1.0; 
         }

         if(i==0 || i==Nout-1) {
            Fout.at(i) = ap*Fin.at(i) + am*Fin.at(i+1);
         } else {
            Fout.at(i) = ap*( a0*Fin.at(i) 
                           +  a1*Fin.at(i+1) 
                           -  a2*Fin.at(i-1) ) 
                       + am*( a0*Fin.at(i+1) 
                           +  a1*Fin.at(i)
                           -  a2*Fin.at(i+2) ); 
         }
      }

   } // end METHOD=QUICK 
   
   else {
      cout << "advection scheme not valid," << endl;
      cout << "should be caught on initilization!" << endl;
      exit (EXIT_FAILURE);
   }


} // end function InterpToCellEdges


void domainGrid::computeFluxTVD(vector<double> &Flout, 
                                vector<double> &FloutL,  vector<double> &FloutR, 
                                vector<double> &Flratio, vector<double> &FlLim, 
                                const vector<double> &Flin,
                                const vector<double> &upC,
                                const vector<double> &fin) const {

   // this function interpolates cell-center flux Flin
   // to Fout, defined at cell edges, using TVD scheme 
   // FloutL and FloutR are left and right going fluxes
   // upC is the local maximum characteristic speed
   // and fin is the function being fluxed
   
   
   const int Nout = Flout.size();
   assert(Nout == nXce);
   assert(FloutL.size()  == Xce.size());
   assert(FloutR.size()  == Xce.size());
   assert(Flratio.size() == Xce.size());
   assert(FlLim.size()   == Xce.size());
   assert(Flin.size() == Xcc.size());
   assert(upC.size()  == Xcc.size());
   assert(fin.size()  == Xcc.size());

 
   vector<double> FluxL1st, FluxR1st;
   vector<double> DeltaFluxL, DeltaFluxR;
   FluxL1st.assign(Nout,0.0);   
   FluxR1st.assign(Nout,0.0); 
   DeltaFluxL.assign(Nout,0.0);
   DeltaFluxR.assign(Nout,0.0);  

   // compute first order left and right fluxes
   //
   for (auto i=0; i<Nout; i++) {
      FluxR1st.at(i) = Flin.at(i)   + upC.at(i)*fin.at(i);
      FluxL1st.at(i) = Flin.at(i+1) - upC.at(i+1)*fin.at(i+1);
   }

   // compute 2nd order left and right fluxes 
   // using flux limiter
   //
   for (auto i=1; i<Nout-1; i++) {
      DeltaFluxR.at(i) = 0.5*(FluxR1st.at(i+1)-FluxR1st.at(i)); 
      DeltaFluxL.at(i) = 0.5*(FluxL1st.at(i)-FluxL1st.at(i-1)); 
      //Flratio.at(i) = DeltaFluxL.at(i)/DeltaFluxR.at(i);
      Flratio.at(i) = (fin.at(i+1) - fin.at(i-1))/(fin.at(i+2) - fin.at(i)); // divide by zero?
      //FlLim.at(i) = (abs(Flratio.at(i))+Flratio.at(i))/(abs(Flratio.at(i)) + 1.0); // van Leer
      FlLim.at(i) = 2.0;
      if(Flratio.at(i)<=2.0) FlLim.at(i) = Flratio.at(i);
      if(Flratio.at(i)<=1.0) FlLim.at(i) = 1.0;
      if(Flratio.at(i)<=0.5) FlLim.at(i) = 2.0*Flratio.at(i);
      if(Flratio.at(i)<=0.0) FlLim.at(i) = 0.0;
      
      FloutR.at(i) = FluxR1st.at(i) + FlLim.at(i)*DeltaFluxR.at(i);
      FloutL.at(i) = FluxL1st.at(i) + FlLim.at(i)*DeltaFluxL.at(i);
      Flout.at(i) = (FloutR.at(i) + FloutL.at(i))/2.0;
   }
   //Flout.at(0) = (FloutR.at(0) + FloutL.at(0))/2.0;
   //Flout.at(nCE-1) = (FluxR.at(nCE-1) + FluxL.at(nCE-1))/2.0;
 

} // end function computeFluxTVD

/*
vector<double> DDX(const vector<double> &Fin, const double &dx) {

   const int Nout  = Fin.size()-1;
   vector<double> Fout;
   Fout.assign(Nout,0.0);
   for (auto i=0; i<Nout; i++) {
      Fout.at(i) = (Fin.at(i+1)-Fin.at(i))/dx;
   }

   return Fout;

}
*/

