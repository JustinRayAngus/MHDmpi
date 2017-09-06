/***
 * 
 * domain grid class
 *
***/

#ifndef domainGrid_h
#define domainGrid_h

//#include "EEDF.h"

using namespace std;

class domainGrid
{

public:
  int nX, nXsub, nXg, nXcc, nXce;
  int numProcs, procID;
  double Xmax, Xmin, dX;
  vector<double> Xcc, Xce; // spatial grid at cell-center and at cell-edge 

  void initialize(const Json::Value&);
  void communicate(vector<double>&) const;
  void DDX(vector<double>&, const vector<double>&) const;

};

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
      cout << "Xmin = " << Xmin << endl;
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

/*
void domainGrid::computeFluxes(EEDF& eedf)
{
   const double dX = Xgrid.dX;
   const int nCE = Flux.size();
   const int nCC = F0.size();

   // compute Flux = F0^2/2 - K*dF0/dX
   //              = FluxAdv + FluxDif
   //
   // FluxAdv = (FluxR+FluxL)/2 is computed using upwinding schemes
   // FluxDif is computed using standard centered scheme


   // Step 1: set flux freezing speed and 
   // compute advection flux at cell center
   
   vector<double> Cspeed, FluxAdvCC, FluxAdv, FluxDif;
   vector<double> DeltaFluxR, DeltaFluxL, FluxRatio, FluxLim;
   FluxAdvCC.assign(nCC,0.0);
   FluxAdv.assign(nCE,0.0);
   DeltaFluxL.assign(nCE,0.0);
   DeltaFluxR.assign(nCE,0.0);
   FluxRatio.assign(nCE,0.0);
   FluxLim.assign(nCE,0.0);
   FluxDif.assign(nCE,0.0);
   //Cspeed.swap(F0); // Local flux freezing speed
   //copy(F0.begin(), F0.end(), Cspeed.begin()); 
   Cspeed = F0; // Local flux freezing speed 
   for (auto i=0; i<nCC; i++) {
      FluxAdvCC.at(i) = F0.at(i)*F0.at(i)/2.0;
   }


   // Step 2: compute first order upwind fluxes
   // for FluxL and FluxR at cell edges
   

   if(advScheme0 == "TVD") {
      
      // Step 2: compute first order upwind fluxes
      // for FluxL and FluxR at cell edges
      //
      for (auto i=0; i<nCE; i++) {
         FluxR.at(i) = FluxAdvCC.at(i)   + Cspeed.at(i)*F0.at(i);
         FluxL.at(i) = FluxAdvCC.at(i+1) - Cspeed.at(i+1)*F0.at(i+1);
         FluxDif.at(i) = -K*(F0.at(i+1)-F0.at(i))/dX;
      }

      // Step 3: compute 2nd order corrections to FluxR and FluxL
      //
      for (auto i=1; i<nCE-1; i++) {
         DeltaFluxR.at(i) = 0.5*(FluxR.at(i+1)-FluxR.at(i)); 
         DeltaFluxL.at(i) = 0.5*(FluxL.at(i)-FluxL.at(i-1)); 
         //FluxRatio.at(i) = DeltaFluxL.at(i)/DeltaFluxR.at(i);
         FluxRatio.at(i) = (F0.at(i) - F0.at(i-1))/(F0.at(i+1) - F0.at(i)); // divide by zero?
         //FluxLim.at(i) = (abs(FluxRatio.at(i))+FluxRatio.at(i))/(abs(FluxRatio.at(i)) + 1.0); // van Leer
         FluxLim.at(i) = 2.0;
         if(FluxRatio.at(i)<=2.0) FluxLim.at(i) = FluxRatio.at(i);
         if(FluxRatio.at(i)<=1.0) FluxLim.at(i) = 1.0;
         if(FluxRatio.at(i)<=0.5) FluxLim.at(i) = 2.0*FluxRatio.at(i);
         
         FluxR.at(i) += FluxLim.at(i)*DeltaFluxR.at(i);
         FluxL.at(i) += FluxLim.at(i)*DeltaFluxL.at(i);
         FluxAdv.at(i) = (FluxR.at(i) + FluxL.at(i))/2.0;
         Flux.at(i) = FluxAdv.at(i) + FluxDif.at(i);
      }
      FluxAdv.at(0) = (FluxR.at(0) + FluxL.at(0))/2.0;
      Flux.at(0) = FluxAdv.at(0) + FluxDif.at(0);
      FluxAdv.at(nCE-1) = (FluxR.at(nCE-1) + FluxL.at(nCE-1))/2.0;
      Flux.at(nCE-1) = FluxAdv.at(nCE-1) + FluxDif.at(nCE-1);
      

   }      
   else if(advScheme0 == "C2") {
      // Use second order central differencing/interpolation
      //
      for (auto i=0; i<nCE; i++) {
         Flux.at(i) = (F0.at(i+1)+F0.at(i))/2.0*(F0.at(i+1)+F0.at(i))/2.0/2.0 
                    - K*(F0.at(i+1)-F0.at(i))/dX;
      }
   } 
   else if(advScheme0 == "U1") {
      // Use first order upwinding
      //
      double Ui, ap, am;
      for (auto i=0; i<nCE; i++) {
      
         Ui = F0.at(i)/2.0;
         ap = 1.0;
         am = 0.0;
         if(Ui<0.0) {
            ap = 0.0;
            am = 1.0; 
         }

         Flux.at(i) = ap*F0.at(i)*F0.at(i)/2.0 
                    + am*F0.at(i+1)*F0.at(i+1)/2.0 
                    - K*(F0.at(i+1)-F0.at(i))/dX;
   
      } // end for loop
   }
   else if(advScheme0 == "QUICK") {   
      // Use 2nd order QUICK upwinding
      //
      double Ui, ap, am;
      double a0 = 3.0/4.0, a1 = 3.0/8.0, a2 = 1.0/8.0;
      for (auto i=0; i<nCE; i++) {
      
         Ui = F0.at(i)/2.0;
         ap = 1.0;
         am = 0.0;
         if(Ui<0.0) {
            ap = 0.0;
            am = 1.0; 
         }

         if(i==0 || i==nCE-1) {
            Flux.at(i) = ap*F0.at(i)*F0.at(i)/2.0 
                       + am*F0.at(i+1)*F0.at(i+1)/2.0 
                       - K*(F0.at(i+1)-F0.at(i))/dX;
   
         } else {
            Flux.at(i) = ap*( a0*F0.at(i)*F0.at(i) 
                           +  a1*F0.at(i+1)*F0.at(i+1) 
                           -  a2*F0.at(i-1)*F0.at(i-1) )/2.0 
                       + am*( a0*F0.at(i+1)*F0.at(i+1) 
                           +  a1*F0.at(i)*F0.at(i)
                           -  a2*F0.at(i+2)*F0.at(i+2) )/2.0 
                       - K*(F0.at(i+1)-F0.at(i))/dX;
         }

      } // end for loop
   } else {
      cout << "advection scheme not valid," << endl;
      cout << "should be caught on initilization!" << endl;
      exit (EXIT_FAILURE);
   }

}
*/

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

   // check that Fin (at cell center)  and Fout (cell edges) are proper size
   
   const int Nout = Fout.size();
   const int Nin  = Fin.size();
   if(Nout != nXce) {  
         cout << "ERROR: output vector in call to domainGrid::DDX " << endl;
         cout << "is not proper size" << endl;
         cout << "Nout = " << Nout << endl;
         cout << "nXce = " << nXce << endl;
         exit (EXIT_FAILURE);
   }
   if(Nin != nXcc) {
         cout << "ERROR: input vector in call to domainGrid::DDX " << endl;
         cout << "is not proper size" << endl;
         exit (EXIT_FAILURE);

   } 


   for (auto i=0; i<Nout; i++) {
      Fout.at(i) = (Fin.at(i+1)-Fin.at(i))/dX;
   }


}

#endif
