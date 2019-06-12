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
#include "matrix2D.h"
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
   for (auto n=0; n<nXce; n++) {
      Xcc.at(n) = Xmin + offset + n*dX;
      Xce.at(n) = Xmin + offset + n*dX + 0.5*dX;
   }
   Xcc.at(nXce) = Xcc.at(nXce-1)+dX;
   //Xcc[nX/numProcs+1] = Xmin + offset + (nXsub+1)*dX;


   ///////////////////////////////////////////////////////////
   //   
   //   get Z-grid information from input file
   //
   ///////////////////////////////////////////////////////////
   
   const Json::Value Zgrid = root.get("Zgrid",defValue);
   if(Zgrid.isObject()) {
      Json::Value ZminVal = Zgrid.get("Zmin",defValue);
      Json::Value ZmaxVal = Zgrid.get("Zmax",defValue);
      Json::Value nZVal   = Zgrid.get("nZ",defValue);
      Json::Value nZgVal  = Zgrid.get("nZg",defValue);
      if(ZminVal == defValue || ZmaxVal == defValue || 
           nZVal == defValue || nZgVal == defValue) {
         printf("ERROR: Zmin, Zmax, nZ, or nZg not declared in input file\n");
         exit (EXIT_FAILURE);
      } 
      nZ = nZVal.asInt();
      if(nZ != nZVal.asDouble() || nZ < 1) {
         printf("ERROR: nZ is not set as a positive integer in input file\n");
         exit (EXIT_FAILURE);
      }
      nZg = nZgVal.asInt();
      if(nZg != nZgVal.asDouble() || nZg < 1 || nZg>4) {
         cout << "ERROR: number of guard cells nZg" << endl;
         cout << "is not set correctly in input file" << endl;
         exit (EXIT_FAILURE);
      }
      Zmin = ZminVal.asDouble();
      Zmax = ZmaxVal.asDouble();
      if(Zmin >= Zmax) {
         cout << "ERROR: Zmin > Zmax in input file" << endl;
         exit (EXIT_FAILURE);
      }
      dZ = (Zmax-Zmin)/(double)nZ;
      nZsub = nZ;
      
      /*
      nZsub = nZ/numProcsZ;
      if(procID==0) {
         double nZsubTest = nZ/(double)numProcsZ;
         if(floor(nZsubTest)==ceil(nZsubTest)) {
            //cout << "nXsub = " << nXsub << endl;
         }
         else {
            printf("ERROR: nZsub=nZ/numProcs is not an integer!!!!\n");
            exit (EXIT_FAILURE); 
         }
      }
      */

      if(procID==0) {
         cout << "Zmin = " << Zmin << endl;
         cout << "Zmax = " << Zmax << endl;
         cout << "nZ = " << nZ << endl;
         cout << "nZsub = " << nZsub << endl;
         cout << "dZ = " << dZ << endl;
         cout << "nZg = " << nZg << endl;
         cout << endl;
      }
      
   }
   else { 
      if(procID==0) { 	   
         cout << "value for key \"Zgrid\" is not object type !" << endl;
         cout << "no Z-grid, problem is 1D" << endl << endl;
      }
      Zmin=0, 
      Zmax=1, 
      nZ=1,
      nZsub=1,
      nZg=1,
      dZ = (Zmax-Zmin)/nZ;
   }
   
   Zcc.assign(nZsub+2*nZg,0.0);
   nZcc = Zcc.size();
   Zce.assign(nZsub+2*nZg-1,0.0);
   nZce = Zce.size();
   //double offsetZ = procID*(Zmax-Zmin)/numProcsZ-(0.5+nZg-1.0)*dZ;
   double offsetZ = -(0.5+nZg-1.0)*dZ;
   for (auto n=0; n<nZce; n++) {
      Zcc.at(n) = Zmin + offsetZ + n*dZ;
      Zce.at(n) = Zmin + offsetZ + n*dZ + 0.5*dZ;
   }
   Zcc.at(nZce) = Zcc.at(nZce-1)+dZ;
   //Zcc[nZ/numProcsZ+1] = Zmin + offset + (nZsub+1)*dZ;
   
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

   } 
   else if(type0=="step") {
      
      int Nmax = Xcc.size();
      for (auto i=0; i<Nmax; i++) {
         (Xcc.at(i)<b) ? var.at(i)=a : var.at(i)=c;
      }   
   
   } 
   else { // quadratic
      Xshift = (Xcc-b);
      var = a*Xshift*Xshift + c*Xshift + d;
   }

} // end setInitialProfile


void domainGrid::setInitialProfile(vector<vector<double>> &var, 
                   const Json::Value &varRoot) const
{
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   double Xa, Xb, Xc, Xd;
   double Za, Zb, Zc, Zd;
   vector<double> Xprofile, Zprofile, Xvec, Zvec;
   string Xtype0, Ztype0;
   Xvec = Xcc;
   Zvec = Zcc;

   const int Nvar0 = var.size();
   const int Nvar1 = var[0].size();
   Xprofile.assign(Nvar0,0.0); 
   Zprofile.assign(Nvar1,0.0); 

   if(Nvar0==nXce) { 
      Xvec.clear();
      Xvec = Xce;
      //Xvec.swap(Xce);
   }   
   if(Nvar1==nZce) {
      Zvec.clear();
      Zvec = Zce;
   }

   //  get X-direction stuff
   //
   const Json::Value defValue; // used for default reference
   Json::Value aVal = varRoot.get("Xa",defValue);
   Json::Value bVal = varRoot.get("Xb",defValue);
   Json::Value cVal = varRoot.get("Xc",defValue);
   Json::Value dVal = varRoot.get("Xd",defValue);
   Json::Value type = varRoot.get("Xtype",defValue);
   if(aVal == defValue || bVal == defValue || 
      cVal == defValue || dVal == defValue || 
      type == defValue) {
      cout << "ERROR: at least 1 initial profile value " << endl;
      cout << "not declared in input file" << endl;
      exit (EXIT_FAILURE);
   } 
   Xa = aVal.asDouble();
   Xb = bVal.asDouble();
   Xc = cVal.asDouble();
   Xd = dVal.asDouble();   
   Xtype0 = type.asString();

   //  set X-profile
   //
   setInitialProfileArbDir(Xprofile,Xvec,Xmin,Xmax,Xa,Xb,Xc,Xd,Xtype0);


   //  get Z-direction stuff
   //
   aVal = varRoot.get("Za",defValue);
   bVal = varRoot.get("Zb",defValue);
   cVal = varRoot.get("Zc",defValue);
   dVal = varRoot.get("Zd",defValue);
   type = varRoot.get("Ztype",defValue);
   if(aVal == defValue || bVal == defValue || 
      cVal == defValue || dVal == defValue || 
      type == defValue) {
      cout << "ERROR: at least 1 initial profile value " << endl;
      cout << "not declared in input file" << endl;
      exit (EXIT_FAILURE);
   } 
   Za = aVal.asDouble();
   Zb = bVal.asDouble();
   Zc = cVal.asDouble();
   Zd = dVal.asDouble();   
   Ztype0 = type.asString();
   
   //  set Z-profile
   //
   setInitialProfileArbDir(Zprofile,Zvec,Zmin,Zmax,Za,Zb,Zc,Zd,Ztype0);

   
   //  multiply X and Z profiles together
   //
   for (auto i=0; i<Nvar0; i++) {
      for (auto j=0; j<Nvar1; j++) {
         var[i][j] = Xprofile[i]*Zprofile[j];
      }
   }

}  // end setInitialProfile

void domainGrid::setInitialProfile(matrix2D<double> &var, 
                   const Json::Value &varRoot) const
{
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   double Xa, Xb, Xc, Xd;
   double Za, Zb, Zc, Zd;
   vector<double> Xprofile, Zprofile, Xvec, Zvec;
   string Xtype0, Ztype0;
   Xvec = Xcc;
   Zvec = Zcc;

   const int Nvar0 = var.size0();
   const int Nvar1 = var.size1();
   Xprofile.assign(Nvar0,0.0); 
   Zprofile.assign(Nvar1,0.0); 

   if(Nvar0==nXce) { 
      Xvec.clear();
      Xvec = Xce;
      //Xvec.swap(Xce);
   }   
   if(Nvar1==nZce) {
      Zvec.clear();
      Zvec = Zce;
   }

   //  get X-direction stuff
   //
   const Json::Value defValue; // used for default reference
   Json::Value aVal = varRoot.get("Xa",defValue);
   Json::Value bVal = varRoot.get("Xb",defValue);
   Json::Value cVal = varRoot.get("Xc",defValue);
   Json::Value dVal = varRoot.get("Xd",defValue);
   Json::Value type = varRoot.get("Xtype",defValue);
   if(aVal == defValue || bVal == defValue || 
      cVal == defValue || dVal == defValue || 
      type == defValue) {
      cout << "ERROR: at least 1 initial profile value " << endl;
      cout << "not declared in input file" << endl;
      exit (EXIT_FAILURE);
   } 
   Xa = aVal.asDouble();
   Xb = bVal.asDouble();
   Xc = cVal.asDouble();
   Xd = dVal.asDouble();   
   Xtype0 = type.asString();

   //  set X-profile
   //
   setInitialProfileArbDir(Xprofile,Xvec,Xmin,Xmax,Xa,Xb,Xc,Xd,Xtype0);


   //  get Z-direction stuff
   //
   aVal = varRoot.get("Za",defValue);
   bVal = varRoot.get("Zb",defValue);
   cVal = varRoot.get("Zc",defValue);
   dVal = varRoot.get("Zd",defValue);
   type = varRoot.get("Ztype",defValue);
   if(aVal == defValue || bVal == defValue || 
      cVal == defValue || dVal == defValue || 
      type == defValue) {
      cout << "ERROR: at least 1 initial profile value " << endl;
      cout << "not declared in input file" << endl;
      exit (EXIT_FAILURE);
   } 
   Za = aVal.asDouble();
   Zb = bVal.asDouble();
   Zc = cVal.asDouble();
   Zd = dVal.asDouble();   
   Ztype0 = type.asString();
   
   //  set Z-profile
   //
   setInitialProfileArbDir(Zprofile,Zvec,Zmin,Zmax,Za,Zb,Zc,Zd,Ztype0);


   //  multiply X and Z profiles together
   //
   for (auto i=0; i<Nvar0; i++) {
      for (auto j=0; j<Nvar1; j++) {
         var(i,j) = Xprofile[i]*Zprofile[j];
      }
   }

}





void domainGrid::setInitialProfileArbDir(vector<double> &var,
	           const vector<double>& Xvec,	
                   const double Xmin, const double Xmax,
                   const double a, const double b,
		   const double c, const double d,
		   const string& type0) const
{
   int procID, numProcs;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   const int thisNvar = var.size();
   const int thisNvec = Xvec.size();
   assert(thisNvar==thisNvec);

   double Lvec = Xmax-Xmin;
   vector<double> Xshift;

   if(type0=="gaussian" || type0=="Gaussian") {
      
      Xshift = (Xvec-b)/c;
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
      
      Xshift = (Xvec-b)/c;
      var = a*tanh(Xshift) + d;
      
      if(c == 0.0) {
         cout << "ERROR: tanh width c cannot be zero" << endl; 
         exit (EXIT_FAILURE);
      }
      if(procID==0) {
         cout << "Initial F0 is tanh with amplitude = " << a << endl;
      }

   } 
   else if(type0=="cos") {
      
      var = a*cos(2.0*3.141592653589793*c*Xvec/Lvec+b) + d;

      if(procID==0) {
         cout << "Initial F0 is cos with mode number  = " << c << endl;
      }

   } 
   else if(type0=="step") {
      
      int Nmax = Xvec.size();
      for (auto i=0; i<Nmax; i++) {
         (Xvec.at(i)<b) ? var.at(i)=a : var.at(i)=c;
      }   
   
   }
   else if(type0=="bennettP") {
      
      Xshift = (Xvec-b)/c;
      int Nmax = Xvec.size();
      for (auto i=0; i<Nmax; i++) {
         var.at(i) = a/pow(1.0+Xshift.at(i)*Xshift.at(i),2);
      }

   } 
   else if(type0=="bennettB") {
      
      Xshift = (Xvec-b)/c;
      int Nmax = Xvec.size();
      for (auto i=0; i<Nmax; i++) {
         var.at(i) = a*sqrt(2)*Xshift.at(i)/(1.0+Xshift.at(i)*Xshift.at(i));
      }

   } 
   else if(type0=="contStratP") {
      
      Xshift = (Xvec-b)/c;
      double thisExponent, s;
      s = d;  // shape factor
      int Nmax = Xvec.size();
      for (auto i=0; i<Nmax; i++) {
         thisExponent = pow(Xshift.at(i),2*s);
         var.at(i) = a*(s+1.0)/(2.0*s)*exp(-thisExponent);
      }

   } 
   else if(type0=="contStratB") {
      
      Xshift = (Xvec-b)/c;
      double thisExponent, s, r, uigamma;
      s = d;   // shape factor
      int Nmax = Xvec.size();
      for (auto i=0; i<Nmax; i++) {
         r = Xshift.at(i);
         thisExponent = pow(r,2*s);
         if(s==1.0) uigamma = exp(-thisExponent);
         if(s==2.0) uigamma = tgamma(1.0/s) - sqrt(3.1415926536)*erf(r*r);
         if(s!=1.0 && s!=2.0) {
            cout << "ERROR: contStratB() requires shape factor (d) = 1 or 2" << endl;
            cout << "have not generalized further yet since don't have incomplete" << endl;
            cout << "gamma function call defined" << endl;
            exit (EXIT_FAILURE);
         }
         var.at(i) = a*(s+1.0)/s*((tgamma(1.0/s)-uigamma)/(s*r*r)-exp(-thisExponent));
         var.at(i) = sqrt(var.at(i));
      }

   } 
   else {
      Xshift = (Xvec-b);
      var = a*Xshift*Xshift + c*Xshift + d;
   }

}


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

void domainGrid::communicate(vector<vector<double>> &F0) const {

   const int nMax = F0.size(); // number of cell-center points
   const int thisnZ = F0[0].size(); // number of Z points
   const int nXce = Xce.size(); // number of cell-center points
   int thisnXg = nXg;
   if(nMax == nXce) thisnXg = nXg-1; // receiving for fluxes 

   int procID, numProcs;
   double Fsend[nXg], Frecv[nXg];

   MPI_Status status;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   for (auto j=0; j<thisnZ; j++) { // loop over z-direction

   //  send F0[nXg+i][j], i<nXg for ID>0 to proc ID-1
   //
   if (procID>0) {
      int tag = 1;
      for (auto i=0; i<thisnXg; i++) {
         Fsend[i] = F0[nXg+i][j];
      }
      MPI_Send(Fsend, thisnXg, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD);
      //MPI_Send(&F0.at(nXg), 1, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD);
   }  
  
   // receive: F0[nXg+i,ID+1][j] => F0[nMax-nXg+i,ID]
   //
   if (procID<numProcs-1) {
      int tag = 1;
      MPI_Recv(Frecv, thisnXg, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD, &status);
      //MPI_Recv(&F0[nXsub+1], 1, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD, &status);
      for (auto i=0; i<thisnXg; i++) {
         F0[nMax-thisnXg+i][j] = Frecv[i];
      }
   }

   // send: F0[nMax-2*nXg+i,ID][j] => ID+1
   //
   if (procID<numProcs-1) {
      int tag = 2;
      for (auto i=0; i<thisnXg; i++) {
         Fsend[i] = F0[nMax-2*nXg+i][j];
         if(nMax==nXce) Fsend[i] = F0[nMax+1-2*nXg+i][j];
      }
      MPI_Send(Fsend, thisnXg, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD);
      //MPI_Send(&F0[nXsub], 1, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD);
   }  

   // receive: F0[nMax-2*nXg+i,ID][j] => F0[i,ID+1][j]
   //
   if (procID>0) {
      int tag = 2;
      MPI_Recv(Frecv, thisnXg, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD, &status);
      //MPI_Recv(&F0[0], 1, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD, &status);
      for (auto i=0; i<thisnXg; i++) {
         F0[i][j] = Frecv[i];
      }
   }

   }

}

void domainGrid::communicate(matrix2D<double> &F0) const {

   const int nMax = F0.size0();     // number of X points
   const int thisnZ = F0.size1();  // number of Z points
   const int nXce = Xce.size(); // number of cell-center points
   int thisnXg = nXg;
   if(nMax == nXce) thisnXg = nXg-1; // receiving for fluxes 

   int procID, numProcs;
   double Fsend[thisnZ], Frecv[thisnZ];
   //double Fsend[nXg], Frecv[nXg];

   MPI_Status status;
   MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   for (auto i=0; i<thisnXg; i++) { // loop over z-direction
   //for (auto j=0; j<thisnZ; j++) { // loop over z-direction

   //  send F0[nXg+i][j], i<nXg for ID>0 to proc ID-1
   //
   if (procID>0) {
      int tag = 1;
      for (auto j=0; j<thisnZ; j++) { // loop over z-direction
         Fsend[j] = F0(nXg+i,j);
      }
      MPI_Send(Fsend, thisnZ, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD);
      //MPI_Send(&F0.at(nXg), 1, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD);
   }  
  
   // receive: F0[nXg+i,ID+1][j] => F0[nMax-nXg+i,ID]
   //
   if (procID<numProcs-1) {
      int tag = 1;
      MPI_Recv(Frecv, thisnZ, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD, &status);
      //MPI_Recv(&F0[nXsub+1], 1, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD, &status);
      for (auto j=0; j<thisnZ; j++) {
         //F0[nMax-thisnXg+i][j] = Frecv[i];
         F0(nMax-thisnXg+i,j) = Frecv[j];
      }
   }

   // send: F0[nMax-2*nXg+i,ID][j] => ID+1
   //
   if (procID<numProcs-1) {
      int tag = 2;
      for (auto j=0; j<thisnZ; j++) {
         Fsend[j] = F0(nMax-2*nXg+i,j);
         //if(nMax==nXce) Fsend[i] = F0[nMax+1-2*nXg+i][j];
         if(nMax==nXce) Fsend[j] = F0(nMax+1-2*nXg+i,j);
      }
      MPI_Send(Fsend, thisnZ, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD);
      //MPI_Send(&F0[nXsub], 1, MPI_DOUBLE, procID+1, tag, MPI_COMM_WORLD);
   }  

   // receive: F0[nMax-2*nXg+i,ID][j] => F0[i,ID+1][j]
   //
   if (procID>0) {
      int tag = 2;
      MPI_Recv(Frecv, thisnZ, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD, &status);
      //MPI_Recv(&F0[0], 1, MPI_DOUBLE, procID-1, tag, MPI_COMM_WORLD, &status);
      for (auto j=0; j<thisnZ; j++) {
         //F0[i][j] = Frecv[i];
         F0(i,j) = Frecv[j];
      }
   }

   }

}


void domainGrid::DDX(vector<double>& Fout, const vector<double>& Fin) const {

   // check that Fin (at cell center) and 
   // use Fout size to determine how to take DDX
   
   const int Nout = Fout.size();
   const int Nin  = Fin.size();
   //assert(Nout == nXce);
   //assert(Nin == nXcc);

   if(Nout == nXce && Nin == nXcc) { // Fout is at cell edges
      for (auto i=0; i<Nout; i++) {
         Fout.at(i) = (Fin.at(i+1)-Fin.at(i))/dX;
      }   
   } 
   else if (Nout == nXcc && Nin == nXcc) {
      for (auto i=1; i<Nout-1; i++) {
         Fout.at(i) = (Fin.at(i+1)-Fin.at(i-1))/2.0/dX;
      }   
   }
   else if (Nout == nXcc && Nin == nXce) {
      for (auto i=1; i<Nout-1; i++) {
         Fout.at(i) = (Fin.at(i)-Fin.at(i-1))/dX;
      }   
   }
   else {
      printf("ERROR: Nout and Nin combo failed for DDX");
      exit (EXIT_FAILURE);
   }

}

void domainGrid::DDX(vector<vector<double>>& Fout, 
		const vector<vector<double>>& Fin) const {

   // check that Fin (at cell center) and 
   // use Fout size to determine how to take DDX
   
   const int Nout0 = Fout.size();
   const int Nout1 = Fout[0].size();
   const int Nin0  = Fin.size();
   const int Nin1  = Fin[0].size();
   assert(Nin0 == nXcc);
   assert(Nin1 == Nout1);

   if(Nout0 == nXce) { // Fout is at cell edges in X-direction
      for (auto i=0; i<Nout0; i++) {
         for (auto j=0; j<Nout1; j++) {
            Fout[i][j] = (Fin[i+1][j]-Fin[i][j])/dX;
	 }
      }   
   } 
   else {  // Fout is at cell center in X-direction
      for (auto i=1; i<Nout0-1; i++) {
         for (auto j=0; j<Nout1; j++) {
            Fout[i][j] = (Fin[i+1][j]-Fin[i-1][j])/2.0/dX;
         }
      }	 
   }

}

void domainGrid::DDX(matrix2D<double>& Fout, 
		const matrix2D<double>& Fin) const {

   // check that Fin (at cell center) and 
   // use Fout size to determine how to take DDX
   
   const int Nout0 = Fout.size0();
   const int Nout1 = Fout.size1();
   const int Nin0  = Fin.size0();
   const int Nin1  = Fin.size1();
   //assert(Nin0 == nXcc);
   assert(Nin1 == Nout1);

   if(Nout0 == nXce && Nin0 == nXcc) { // Fout is at cell edges in X-direction
      for (auto i=0; i<Nout0; i++) {
         for (auto j=0; j<Nout1; j++) {
            Fout(i,j) = (Fin(i+1,j)-Fin(i,j))/dX;
	 }
      }   
   } 
   else if (Nout0 == nXcc && Nin0 == nXcc) {  // Fout is at cell center in X-direction
      for (auto i=1; i<Nout0-1; i++) {
         for (auto j=0; j<Nout1; j++) {
            Fout(i,j) = (Fin(i+1,j)-Fin(i-1,j))/2.0/dX;
         }
      }	 
   }
   else if (Nout0 == nXcc && Nin0 == nXce) {
      for (auto i=1; i<Nout0-1; i++) {
         for (auto j=0; j<Nout1; j++) {
            Fout(i,j) = (Fin(i,j)-Fin(i-1,j))/dX;
         }   
      }
   }
   else {
      printf("ERROR: Nout and Nin combo failed for DDX");
      exit (EXIT_FAILURE);
   }


}

void domainGrid::DDZ(vector<vector<double>>& Fout, 
		const vector<vector<double>>& Fin) const {

   // check that Fin (at cell center) and 
   // use Fout size to determine how to take DDZ
   
   const int Nout0 = Fout.size();
   const int Nout1 = Fout[0].size();
   const int Nin0  = Fin.size();
   const int Nin1  = Fin[0].size();
   assert(Nout0 == Nin0);
   assert(Nin1 == nZcc);

   if(Nout1 == nZce) { // Fout is at cell edges in Z-direction
      for (auto i=0; i<Nout0; i++) {
         for (auto j=0; j<Nout1; j++) {
            Fout[i][j] = (Fin[i][j+1]-Fin[i][j])/dZ;
	 }
      }   
   } 
   else {  // Fout is at cell center in Z-direction
      for (auto i=0; i<Nout0; i++) {
         for (auto j=1; j<Nout1-1; j++) {
            Fout[i][j] = (Fin[i][j+1]-Fin[i][j-1])/2.0/dZ;
         }
      }	 
   }

}

void domainGrid::DDZ(matrix2D<double>& Fout, 
		const matrix2D<double>& Fin) const {

   // check that Fin (at cell center) and 
   // use Fout size to determine how to take DDZ
   
   const int Nout0 = Fout.size0();
   const int Nout1 = Fout.size1();
   const int Nin0  = Fin.size0();
   const int Nin1  = Fin.size1();
   assert(Nout0 == Nin0);
   //assert(Nin1 == nZcc);

   if(Nout1 == nZce && Nin1 == nZcc) { // Fout is at cell edges in Z-direction
      for (auto i=0; i<Nout0; i++) {
         for (auto j=0; j<Nout1; j++) {
            Fout(i,j) = (Fin(i,j+1)-Fin(i,j))/dZ;
	 }
      }   
   } 
   else if(Nout1 == nZcc && Nin1 == nZcc) {  // Fout is at cell center in Z-direction
      for (auto i=0; i<Nout0; i++) {
         for (auto j=1; j<Nout1-1; j++) {
            Fout(i,j) = (Fin(i,j+1)-Fin(i,j-1))/2.0/dZ;
         }
      }	 
   }
   else if (Nout1 == nZcc && Nin1 == nZce) {
      for (auto i=0; i<Nout0; i++) {
         for (auto j=1; j<Nout1-1; j++) {
            Fout(i,j) = (Fin(i,j)-Fin(i,j-1))/dZ;
         }   
      }
   }
   else {
      printf("ERROR: Nout and Nin combo failed for DDZ");
      exit (EXIT_FAILURE);
   }


}

void domainGrid::D2DZ2(matrix2D<double>& Fout, 
    	 	 const matrix2D<double>& Fin) const {

   // check that Fin and Fout at cell center 
   
   const int Nout0 = Fout.size0();
   const int Nout1 = Fout.size1();
   const int Nin0  = Fin.size0();
   const int Nin1  = Fin.size1();
   assert(Nout0 == Nin0);
   assert(Nin1 == nZcc);
   assert(Nout1 == nZcc);

   for (auto i=0; i<Nout0; i++) {
      for (auto j=1; j<Nout1-1; j++) {
         Fout(i,j) = (Fin(i,j+1) - 2.0*Fin(i,j) + Fin(i,j-1))/dZ/dZ;
      }
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
      
      assert(nXg >= 2);
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
   
   else if(METHOD == "WENO5") {   
      
      assert(nXg >= 3);
      // Use 5th order WENO upwinding (see Ren 2003)
      //
      double ep, d0, d1, d2, f0, f1, f2, b0, b1, b2;
      double df, alpha, w0, w1, w2, w0b, w1b, w2b, sumwb;
      ep = 1.0e-6;

      for (auto i=nXg-1; i<Nout-nXg+1; i++) {
      
         //df = Fin.at(i+1)/upC.at(i+1)-Fin.at(i)/upC.at(i);
         //if(df==df && df != 0.0) {
         //   alpha = (Fin.at(i+1)-Fin.at(i))/df;
         //}
         //else {
         //   alpha = upC.at(i+1)+upC.at(i);
         //}
         alpha = upC.at(i+1)+upC.at(i);
        
	 if(alpha>=0.0) {
            d0 = 0.3;
            d1 = 0.6;
            d2 = 0.1;
	    //
            f0 = (2.0*Fin.at(i)+5.0*Fin.at(i+1)-1.0*Fin.at(i+2))/6.0;
            f1 = (-1.0*Fin.at(i-1)+5.0*Fin.at(i)+2.0*Fin.at(i+1))/6.0;
            f2 = (2.0*Fin.at(i-2)-7.0*Fin.at(i-1)+11.0*Fin.at(i))/6.0;
            //
            b0 = 13.0/12.0*pow(Fin.at(i)-2.0*Fin.at(i+1)+Fin.at(i+2),2)
               + 1.0/4.0*pow(3.0*Fin.at(i)-4.0*Fin.at(i+1)+Fin.at(i+2),2);
            b1 = 13.0/12.0*pow(Fin.at(i-1)-2.0*Fin.at(i)+Fin.at(i+1),2)
               + 1.0/4.0*pow(Fin.at(i-1)-Fin.at(i+1),2);
            b2 = 13.0/12.0*pow(Fin.at(i-2)-2.0*Fin.at(i-1)+Fin.at(i),2)
               + 1.0/4.0*pow(Fin.at(i-2)-4.0*Fin.at(i-1)+3.0*Fin.at(i),2);
	 }
	 else {
            d0 = 0.1;
            d1 = 0.6;
            d2 = 0.3;
	    //
            f0 = (11.0*Fin.at(i+1)-7.0*Fin.at(i+2)+2.0*Fin.at(i+3))/6.0;
            f1 = (2.0*Fin.at(i)+5.0*Fin.at(i+1)-1.0*Fin.at(i+2))/6.0;
            f2 = (-1.0*Fin.at(i-1)+5.0*Fin.at(i)+2.0*Fin.at(i+1))/6.0;
            //
            b0 = 13.0/12.0*pow(Fin.at(i+1)-2.0*Fin.at(i+2)+Fin.at(i+3),2)
               + 1.0/4.0*pow(3.0*Fin.at(i+1)-4.0*Fin.at(i+2)+Fin.at(i+3),2);
            b1 = 13.0/12.0*pow(Fin.at(i)-2.0*Fin.at(i+1)+Fin.at(i+2),2)
               + 1.0/4.0*pow(Fin.at(i)-Fin.at(i+2),2);
            b2 = 13.0/12.0*pow(Fin.at(i-1)-2.0*Fin.at(i)+Fin.at(i+1),2)
               + 1.0/4.0*pow(Fin.at(i-1)-4.0*Fin.at(i)+3.0*Fin.at(i+1),2);
	 }

	 w0b = d0/pow(ep+b0,2);
	 w1b = d1/pow(ep+b1,2);
	 w2b = d2/pow(ep+b2,2);

         sumwb = w0b+w1b+w2b;
         w0 = w0b/sumwb;
         w1 = w1b/sumwb;
         w2 = w2b/sumwb;

         Fout.at(i) = w0*f0 + w1*f1 + w2*f2; 
      
      }

   } // end METHOD=WENO5
   
   else {
      cout << "advection scheme not valid," << endl;
      cout << "should be caught on initilization!" << endl;
      exit (EXIT_FAILURE);
   }


} // end function InterpToCellEdges


void domainGrid::InterpToCellEdges(vector<vector<double>> &Fout, 
                             const vector<vector<double>> &Fin,
                             const vector<vector<double>> &upC,
                             const string& METHOD,
			     const int dir) const {

   // this function interpolates cell-center Fin to cell
   // to Fout, defined at cell edges
   // upC is the local maximum characteristic speed
   // METHOD referes to interpolation scheme


   // check that vectors in call are proper size
   //
   const int Nout0 = Fout.size();
   const int Nout1 = Fout[0].size();
   const int Nin0  = Fin.size();
   const int Nin1  = Fin[0].size();
   if(dir==0) {
      assert(Nin0==nXcc);  // input cell-center in X 
      assert(Nout0==nXce); // output cell-edge in X
      assert(Nin1==Nout1); // input and output same size in Z
   }   
   else if(dir==1) {
      assert(Nin1==nZcc);  // input cell-center in Z 
      assert(Nout1==nZce); // output cell-edge in Z
      assert(Nin0==Nout0); // input and output same size in X
   }   
   else {
      cout << "dir in call to InterpToCellEdges most be 0 or 1" << endl;
      exit (EXIT_FAILURE);
   }

   //  interpolate using specified method
   //
   if(METHOD == "C2") { // 2nd order central
   
      if(dir==0) { // interp to X-faces
         for (auto i=0; i<Nout0; i++) {
            for (auto j=0; j<Nout1; j++) {
               Fout[i][j] = (Fin[i+1][j] + Fin[i][j])/2.0;
            }
         }
      }
      else {       // interp to Z-faces
         for (auto i=0; i<Nout0; i++) {
            for (auto j=0; j<Nout1; j++) {
               Fout[i][j] = (Fin[i][j+1] + Fin[i][j])/2.0;
            }
         }
      }

   } // end METHOD=C2

   else if(METHOD == "U1") { // first order upwind
  
      if(dir==0) { // interp to X-faces
         for (auto i=0; i<Nout0; i++) {
            for (auto j=0; j<Nout1; j++) {
      
               if(upC[i][j]<0.0) {
                  Fout[i][j] = Fin[i+1][j];
               } else {
                  Fout[i][j] = Fin[i][j];
               }

	    }
         }
      }
      else { // interp to Z-faces
         for (auto i=0; i<Nout0; i++) {
            for (auto j=0; j<Nout1; j++) {
      
               if(upC[i][j]<0.0) {
                  Fout[i][j] = Fin[i][j+1];
               } else {
                  Fout[i][j] = Fin[i][j];
               }

	    }
         }
      }
   
   } // end METHOD=U1
   
   else {
      cout << "advection scheme not valid," << endl;
      cout << "should be caught on initilization!" << endl;
      exit (EXIT_FAILURE);
   }

}

void domainGrid::InterpToCellEdges(matrix2D<double> &Fout, 
                             const matrix2D<double> &Fin,
                             const matrix2D<double> &upC,
                             const string& METHOD,
			     const int dir) const {

   // this function interpolates cell-center Fin to cell
   // to Fout, defined at cell edges
   // upC is the local maximum characteristic speed
   // METHOD referes to interpolation scheme


   // check that vectors in call are proper size
   //
   const int Nout0 = Fout.size0();
   const int Nout1 = Fout.size1();
   const int Nin0  = Fin.size0();
   const int Nin1  = Fin.size1();
   if(dir==0) {
      assert(Nin0==nXcc);  // input cell-center in X 
      assert(Nout0==nXce); // output cell-edge in X
      assert(Nin1==Nout1); // input and output same size in Z
   }   
   else if(dir==1) {
      assert(Nin1==nZcc);  // input cell-center in Z 
      assert(Nout1==nZce); // output cell-edge in Z
      assert(Nin0==Nout0); // input and output same size in X
   }   
   else {
      cout << "dir in call to InterpToCellEdges most be 0 or 1" << endl;
      exit (EXIT_FAILURE);
   }

   //  interpolate using specified method
   //
   if(METHOD == "C2") { // 2nd order central
   
      if(dir==0) { // interp to X-faces
         for (auto i=0; i<Nout0; i++) {
            for (auto j=0; j<Nout1; j++) {
               Fout(i,j) = (Fin(i+1,j) + Fin(i,j))/2.0;
            }
         }
      }
      else {       // interp to Z-faces
         for (auto i=0; i<Nout0; i++) {
            for (auto j=0; j<Nout1; j++) {
               Fout(i,j) = (Fin(i,j+1) + Fin(i,j))/2.0;
            }
         }
      }

   } // end METHOD=C2

   else if(METHOD == "U1") { // first order upwind
  
      if(dir==0) { // interp to X-faces
         for (auto i=0; i<Nout0; i++) {
            for (auto j=0; j<Nout1; j++) {
      
               if(upC(i,j)<0.0) {
                  Fout(i,j) = Fin(i+1,j);
               } else {
                  Fout(i,j) = Fin(i,j);
               }

	    }
         }
      }
      else { // interp to Z-faces
         for (auto i=0; i<Nout0; i++) {
            for (auto j=0; j<Nout1; j++) {
      
               if(upC(i,j)<0.0) {
                  Fout(i,j) = Fin(i,j+1);
               } else {
                  Fout(i,j) = Fin(i,j);
               }

	    }
         }
      }
   
   } // end METHOD=U1
   
   else if(METHOD == "QUICK") { // 2nd? order upwind
      
      double Ui, ap, am;
      double a0 = 3.0/4.0, a1 = 3.0/8.0, a2 = 1.0/8.0;   

      if(dir==0) { // interp to X-faces
         assert(nXg >= 2);
         for (auto i=0; i<Nout0; i++) {
            for (auto j=0; j<Nout1; j++) {
               Ui = upC(i,j);
               ap = 1.0;
               am = 0.0;
               if(Ui<0.0) {
                   ap = 0.0;
                   am = 1.0;
               }

               if(i==0.0 || i==Nout0-1) {
                  Fout(i,j) = ap*Fin(i,j) + am*Fin(i+1,j);
               } else {
                  Fout(i,j) = ap*( a0*Fin(i,j)
                                +  a1*Fin(i+1,j)
                                -  a2*Fin(i-1,j) )
                            + am*( a0*Fin(i+1,j)
                                +  a1*Fin(i,j)
                                -  a2*Fin(i+2,j) );
               }

	    }
         }
      }
      else { // interp to Z-faces
         assert(nZg >= 2);
         for (auto i=0; i<Nout0; i++) {
            for (auto j=0; j<Nout1; j++) {
               Ui = upC(i,j);
               ap = 1.0;
               am = 0.0;
               if(Ui<0.0) {
                   ap = 0.0;
                   am = 1.0;
               }
      
               if(j==0.0 || j==Nout1-1) {
                  Fout(i,j) = ap*Fin(i,j) + am*Fin(i,j+1);
               } else {
                  Fout(i,j) = ap*( a0*Fin(i,j)
                                +  a1*Fin(i,j+1)
                                -  a2*Fin(i,j-1) )
                            + am*( a0*Fin(i,j+1)
                                +  a1*Fin(i,j)
                                -  a2*Fin(i,j+2) );
               }

	    }
         }
      }
   
   } // end METHOD=QUICK
   
   else if(METHOD == "WENO5") {  
      // Use 5th order WENO upwinding (see Ren 2003)
      //
      double ep, d0, d1, d2, f0, f1, f2, b0, b1, b2;
      double df, alpha, w0, w1, w2, w0b, w1b, w2b, sumwb;
      ep = 1.0e-6;

      if(dir==0) { // interp to X-faces

         assert(nXg >= 3);

         for (auto i=nXg-1; i<Nout0-nXg+1; i++) {
            for (auto j=0; j<Nout1; j++) {
               
               alpha = upC(i+1,j)+upC(i,j);  
               if(alpha>=0.0) {
                  d0 = 0.3;
                  d1 = 0.6;
                  d2 = 0.1;
                  //
                  f0 = (2.0*Fin(i,j) + 5.0*Fin(i+1,j) - 1.0*Fin(i+2,j))/6.0;
                  f1 = (-1.0*Fin(i-1,j) + 5.0*Fin(i,j) + 2.0*Fin(i+1,j))/6.0;
                  f2 = (2.0*Fin(i-2,j) - 7.0*Fin(i-1,j) + 11.0*Fin(i,j))/6.0;
                  //
                  b0 = 13.0/12.0*pow(Fin(i,j)-2.0*Fin(i+1,j)+Fin(i+2,j),2)
                     + 1.0/4.0*pow(3.0*Fin(i,j)-4.0*Fin(i+1,j)+Fin(i+2,j),2);
                  b1 = 13.0/12.0*pow(Fin(i-1,j)-2.0*Fin(i,j)+Fin(i+1,j),2)
                     + 1.0/4.0*pow(Fin(i-1,j)-Fin(i+1,j),2);
                  b2 = 13.0/12.0*pow(Fin(i-2,j)-2.0*Fin(i-1,j)+Fin(i,j),2)
                     + 1.0/4.0*pow(Fin(i-2,j)-4.0*Fin(i-1,j)+3.0*Fin(i,j),2);
               } else {
                  d0 = 0.1;
                  d1 = 0.6;
                  d2 = 0.3;
                  //
                  f0 = (11.0*Fin(i+1,j) - 7.0*Fin(i+2,j) + 2.0*Fin(i+3,j))/6.0;
                  f1 = (2.0*Fin(i,j) + 5.0*Fin(i+1,j) - 1.0*Fin(i+2,j))/6.0;
                  f2 = (-1.0*Fin(i-1,j) + 5.0*Fin(i,j) + 2.0*Fin(i+1,j))/6.0;
                  //
                  b0 = 13.0/12.0*pow(Fin(i+1,j)-2.0*Fin(i+2,j)+Fin(i+3,j),2)
                     + 1.0/4.0*pow(3.0*Fin(i+1,j)-4.0*Fin(i+2,j)+Fin(i+3,j),2);
                  b1 = 13.0/12.0*pow(Fin(i,j)-2.0*Fin(i+1,j)+Fin(i+2,j),2)
                     + 1.0/4.0*pow(Fin(i,j)-Fin(i+2,j),2);
                  b2 = 13.0/12.0*pow(Fin(i-1,j)-2.0*Fin(i,j)+Fin(i+1,j),2)
                     + 1.0/4.0*pow(Fin(i-1,j)-4.0*Fin(i,j)+3.0*Fin(i+1,j),2);
               }

               w0b = d0/pow(ep+b0,2);
               w1b = d1/pow(ep+b1,2);
               w2b = d2/pow(ep+b2,2);

               sumwb = w0b+w1b+w2b;
               w0 = w0b/sumwb;
               w1 = w1b/sumwb;
               w2 = w2b/sumwb;

               Fout(i,j) = w0*f0 + w1*f1 + w2*f2;
         
	    }
         }

      }
      else { // interp to Z-faces

         assert(nZg >= 3);

         for (auto j=nZg-1; j<Nout1-nZg+1; j++) {
            for (auto i=0; i<Nout0; i++) {

               alpha = upC(i,j+1)+upC(i,j);  
               if(alpha>=0.0) {
                  d0 = 0.3;
                  d1 = 0.6;
                  d2 = 0.1;
                  //
                  f0 = (2.0*Fin(i,j) + 5.0*Fin(i,j+1) - 1.0*Fin(i,j+2))/6.0;
                  f1 = (-1.0*Fin(i,j-1) + 5.0*Fin(i,j) + 2.0*Fin(i,j+1))/6.0;
                  f2 = (2.0*Fin(i,j-2) - 7.0*Fin(i,j-1) + 11.0*Fin(i,j))/6.0;
                  //
                  b0 = 13.0/12.0*pow(Fin(i,j)-2.0*Fin(i,j+1)+Fin(i,j+2),2)
                     + 1.0/4.0*pow(3.0*Fin(i,j)-4.0*Fin(i,j+1)+Fin(i,j+2),2);
                  b1 = 13.0/12.0*pow(Fin(i,j-1)-2.0*Fin(i,j)+Fin(i,j+1),2)
                     + 1.0/4.0*pow(Fin(i,j-1)-Fin(i,j+1),2);
                  b2 = 13.0/12.0*pow(Fin(i,j-2)-2.0*Fin(i,j-1)+Fin(i,j),2)
                     + 1.0/4.0*pow(Fin(i,j-2)-4.0*Fin(i,j-1)+3.0*Fin(i,j),2);
               } else {
                  d0 = 0.1;
                  d1 = 0.6;
                  d2 = 0.3;
                  //
                  f0 = (11.0*Fin(i,j+1) - 7.0*Fin(i,j+2) + 2.0*Fin(i,j+3))/6.0;
                  f1 = (2.0*Fin(i,j) + 5.0*Fin(i,j+1) - 1.0*Fin(i,j+2))/6.0;
                  f2 = (-1.0*Fin(i,j-1) + 5.0*Fin(i,j) + 2.0*Fin(i,j+1))/6.0;
                  //
                  b0 = 13.0/12.0*pow(Fin(i,j+1)-2.0*Fin(i,j+2)+Fin(i,j+3),2)
                     + 1.0/4.0*pow(3.0*Fin(i,j+1)-4.0*Fin(i,j+2)+Fin(i,j+3),2);
                  b1 = 13.0/12.0*pow(Fin(i,j)-2.0*Fin(i,j+1)+Fin(i,j+2),2)
                     + 1.0/4.0*pow(Fin(i,j)-Fin(i,j+2),2);
                  b2 = 13.0/12.0*pow(Fin(i,j-1)-2.0*Fin(i,j)+Fin(i,j+1),2)
                     + 1.0/4.0*pow(Fin(i,j-1)-4.0*Fin(i,j)+3.0*Fin(i,j+1),2);
               }

               w0b = d0/pow(ep+b0,2);
               w1b = d1/pow(ep+b1,2);
               w2b = d2/pow(ep+b2,2);

               sumwb = w0b+w1b+w2b;
               w0 = w0b/sumwb;
               w1 = w1b/sumwb;
               w2 = w2b/sumwb;

               Fout(i,j) = w0*f0 + w1*f1 + w2*f2;

	    }
         }

      }
   
   } // end METHOD=WENO5
   
   else {
      cout << "advection scheme not valid," << endl;
      cout << "should be caught on initilization!" << endl;
      exit (EXIT_FAILURE);
   }

}



void domainGrid::InterpToCellCenter(vector<double> &Fout, 
                                   const vector<double> &Fin) const {

   // this function interpolates cell-edge Fin to cell
   // to Fout, defined at cell center

   // check that vectors in call are proper size
   //
   const int Nout = Fout.size();
   const int Nin  = Fin.size();
   assert(Nout == nXcc);
   assert(Nin  == nXce);


   for (auto i=1; i<Nout-1; i++) {
      Fout.at(i) = (Fin.at(i)+Fin.at(i-1))/2.0;
   }


} // end InterpToCellCenter


void domainGrid::InterpToCellCenter(vector<vector<double>> &Fout, 
                                    const vector<vector<double>> &Fin) const {

   // this function interpolates cell-edge Fin
   // to Fout, defined at cell center

   // check that vectors in call are proper size
   //
   const int Nout0 = Fout.size();
   const int Nout1 = Fout[0].size();
   const int Nin0  = Fin.size();
   const int Nin1  = Fin[0].size();
   assert(Nout0 == nXcc);
   assert(Nout1 == nZcc);
   assert(Nin0==nXce || Nin1==nZce); // at least one direction stag

   if(Nin0==nXce) { // input stag in X-direction
   
      for (auto i=1; i<Nout0-1; i++) {
         for (auto j=0; j<Nout1; j++) {
            Fout[i][j] = (Fin[i][j] + Fin[i-1][j])/2.0;
         }
      }

   }
   else {  // input stag in Z-direction
      
      for (auto i=0; i<Nout0; i++) {
         for (auto j=1; j<Nout1-1; j++) {
            Fout[i][j] = (Fin[i][j] + Fin[i][j-1])/2.0;
         }
      }

   }


} // end InterpToCellCenter

void domainGrid::InterpToCellCenter(matrix2D<double> &Fout, 
                                    const matrix2D<double> &Fin) const {

   // this function interpolates cell-edge Fin
   // to Fout, defined at cell center

   // check that vectors in call are proper size
   //
   const int Nout0 = Fout.size0();
   const int Nout1 = Fout.size1();
   const int Nin0  = Fin.size0();
   const int Nin1  = Fin.size1();
   //cout << " Nout0 = " << Nout0 << endl;
   //cout << " Nout1 = " << Nout1 << endl;
   //cout << " nXcc = " << nXcc << endl;
   //cout << " nZcc = " << nZcc << endl;
   assert(Nout0 == nXcc);
   assert(Nout1 == nZcc);
   assert(Nin0==nXce || Nin1==nZce); // at least one direction stag

   if(Nin0==nXce && Nin1==nZcc) { // input stag in X-direction
   
      for (auto i=1; i<Nout0-1; i++) {
         for (auto j=0; j<Nout1; j++) {
            Fout(i,j) = (Fin(i,j) + Fin(i-1,j))/2.0;
         }
      }

   }
   if(Nin0==nXcc && Nin1==nZce) { // input stag in Z-direction
   //else {  // input stag in Z-direction
      
      for (auto i=0; i<Nout0; i++) {
         for (auto j=1; j<Nout1-1; j++) {
            Fout(i,j) = (Fin(i,j) + Fin(i,j-1))/2.0;
         }
      }

   }
   if(Nin0==nXce && Nin1==nZce) { // input stag in both directions
   
      for (auto i=1; i<Nout0-1; i++) {
         for (auto j=1; j<Nout1-1; j++) {
            Fout(i,j) = (Fin(i,j) + Fin(i-1,j) + Fin(i,j-1) + Fin(i-1,j-1))/4.0;
         }
      }

   }


} // end InterpToCellCenter


void domainGrid::computeFluxTVD(vector<double> &Flout, 
                                vector<double> &FloutL,  vector<double> &FloutR, 
                                vector<double> &Flratio, vector<double> &FlLimR, 
                                const vector<double> &Flin,
                                const vector<double> &upC,
                                const vector<double> &fin,
                                const string &LIMITER,
				const int order) const {

   // this function interpolates cell-center flux Flin
   // to Flout, defined at cell edges, using TVD scheme 
   // FloutL and FloutR are left and right going fluxes
   // upC is the local maximum characteristic speed
   // and fin is the function being fluxed
   
   
   const int Nout = Flout.size();
   assert(Nout == nXce);
   assert(FloutL.size()  == Xce.size());
   assert(FloutR.size()  == Xce.size());
   assert(Flratio.size() == Xce.size());
   assert(FlLimR.size()   == Xce.size());
   assert(Flin.size() == Xcc.size());
   assert(upC.size()  == Xcc.size());
   assert(fin.size()  == Xcc.size());

 
   vector<double> FluxL1st, FluxR1st;
   double DeltaFluxRL, DeltaFluxRR;
   double DeltaFluxLL, DeltaFluxLR;
   FluxL1st.assign(Nout,0.0);   
   FluxR1st.assign(Nout,0.0); 
   
   vector<double> FlLimL, FlinR, FlinL;
   FlLimL.assign(Nout,0.0);   
   FlinR.assign(nXcc,0.0);   
   FlinL.assign(nXcc,0.0);   

   FlinR = 0.5*(Flin + upC*fin);
   FlinL = 0.5*(Flin - upC*fin);


   // compute first order left and right fluxes
   //
   for (auto i=0; i<Nout; i++) {
      FluxR1st.at(i) = FlinR.at(i);
      FluxL1st.at(i) = FlinL.at(i+1);
   }

   // compute 2nd order left and right flux corrections 
   // using flux limiter
   //
   for (auto i=1; i<Nout-1; i++) {
    
      assert(nXg >= 2);

      // calculate flux corrections for right going wave
      //
      DeltaFluxRL = 0.5*(FlinR.at(i)   - FlinR.at(i-1)); // B2 scheme 
      DeltaFluxRR = 0.5*(FlinR.at(i+1) - FlinR.at(i));   // C2 scheme

      // calculate flux corrections for left going wave
      //
      DeltaFluxLL = -0.5*(FlinL.at(i+1) - FlinL.at(i));   // C2 scheme 
      DeltaFluxLR = -0.5*(FlinL.at(i+2) - FlinL.at(i+1)); // F2 scheme
      
      // vanleer(a,b)  = 2*a*b/(a+b) (a and b same sign) 
      //               = 0 otherwise
      // minmod(a,b)   = 1/2(sign(a) + sign(b))*min(|a|,|b|)
      // superbee(a,b) = minmod(a,2b) if |a|>=|b|
      //               = minmod(2a,b) if otherwise
      //
      if(LIMITER=="minmod") { 
          FlLimR.at(i) = minmod(DeltaFluxRL,DeltaFluxRR);
          FlLimL.at(i) = minmod(DeltaFluxLL,DeltaFluxLR);
      } else if(LIMITER=="superbee") {       
          FlLimR.at(i) = superbee(DeltaFluxRL,DeltaFluxRR);
	  FlLimL.at(i) = superbee(DeltaFluxLL,DeltaFluxLR);
      } else {
	  FlLimR.at(i) = vanleer(DeltaFluxRL,DeltaFluxRR);
	  FlLimL.at(i) = vanleer(DeltaFluxLL,DeltaFluxLR);
      }

      if(order==1) {
         FloutR.at(i) = FluxR1st.at(i);
         FloutL.at(i) = FluxL1st.at(i);
      }
      else {
         FloutR.at(i) = FluxR1st.at(i) + FlLimR.at(i);
         FloutL.at(i) = FluxL1st.at(i) + FlLimL.at(i);
      }
      Flout.at(i) = FloutR.at(i) + FloutL.at(i);
   }
   //Flout.at(0) = FloutR.at(0) + FloutL.at(0);
   //Flout.at(nCE-1) = FluxR.at(nCE-1) + FluxL.at(nCE-1);
   

} // end function computeFluxTVD

void domainGrid::computeFluxTVDsimple(vector<double> &Flout, 
                                vector<double> &FloutL,  vector<double> &FloutR, 
                                const vector<double> &Flin,
                                const vector<double> &upC,
                                const vector<double> &fin,
                                const string& METHOD) const {

   // this function interpolates cell-center flux Flin
   // to Flout, defined at cell edges, using TVD scheme 
   // FloutL and FloutR are left and right going fluxes
   // upC is the local maximum characteristic speed
   // and fin is the function being fluxed
   
   
   const int Nout = Flout.size();
   assert(Nout == nXce);
   assert(FloutL.size()  == Xce.size());
   assert(FloutR.size()  == Xce.size());
   assert(Flin.size() == Xcc.size());
   assert(upC.size()  == Xcc.size());
   assert(fin.size()  == Xcc.size());

   vector<double> FlinR, FlinL;
   FlinR.assign(nXcc,0.0);   
   FlinL.assign(nXcc,0.0);   

   FlinR = 0.5*(Flin + upC*fin);
   FlinL = 0.5*(Flin - upC*fin);

   InterpToCellEdges(FloutR,FlinR,upC,METHOD);
   InterpToCellEdges(FloutL,FlinL,-upC,METHOD);

   Flout = FloutR + FloutL;

}

void domainGrid::computeFluxTVDnew(matrix2D<double> &Flout, 
                                matrix2D<double> &FloutL,
				matrix2D<double> &FloutR, 
                                matrix2D<double> &FlLimL, 
                                matrix2D<double> &FlLimR, 
                                const matrix2D<double> &Flin,
                                const matrix2D<double> &upC,
                                const matrix2D<double> &fin,
				const int dir, 
				const int order) const {
   
   // same as computeFluxTVD, but set up to use different arguements
   // for flux limiter calculation

   // this function interpolates cell-center flux Flin
   // to Flout, defined at cell edges, using TVD scheme 
   // FloutL and FloutR are left and right going fluxes
   // upC is the local maximum characteristic speed
   // and fin is the function being fluxed
   // dir is the directoin (0 for X, 1 for Z)
   // order is scheme order (1 or 2)
   
   
   const int Nout0 = Flout.size0();
   const int Nout1 = Flout.size1();
   const int Nin0 = Flin.size0();
   const int Nin1 = Flin.size1();

   if(dir==0) {
      assert(Nout0 == nXce && Nin0 == nXcc);
   }
   else if(dir==1) {
      assert(Nout1 == nZce && Nin1 == nZcc);
   }
   else {
      cout << "dir in call to computeFluxTVD most be 0 or 1" << endl;
      exit (EXIT_FAILURE);
   }
   assert(upC.size0() == Flin.size0());
   assert(fin.size0() == Flin.size0());
   assert(upC.size1() == Flin.size1());
   assert(fin.size1() == Flin.size1());

   //  define some additional matrix values for calculation
   //
   matrix2D<double> FluxL1st(Nout0,Nout1,0.0); 
   matrix2D<double> FluxR1st(Nout0,Nout1,0.0); 
   matrix2D<double> FluxLC2(Nout0,Nout1,0.0); 
   matrix2D<double> FluxRC2(Nout0,Nout1,0.0); 
   double ratioR, ratioL;
   matrix2D<double> FlinR(Nin0,Nin1,0.0);
   matrix2D<double> FlinL(Nin0,Nin1,0.0);

   FlinR = 0.5*(Flin + upC*fin);
   FlinL = 0.5*(Flin - upC*fin);
   int iup=1, jup=0;
   if(dir==1) iup=0, jup=1;

   //  compute U1 and C2 left and right fluxes
   //
   for (auto i=0; i<Nout0; i++) {
      for (auto j=0; j<Nout1; j++) {
         FluxR1st(i,j) = FlinR(i,j);
         FluxL1st(i,j) = FlinL(i+iup,j+jup);
         FluxRC2(i,j) = (FlinR(i,j) + FlinR(i+iup,j+jup))/2.0;
         FluxLC2(i,j) = (FlinL(i,j) + FlinL(i+iup,j+jup))/2.0;
      }
   }
   
   //  compute second order left and right flux corrections
   //  using flux limiter
   //
   for (auto i=1; i<Nout0-1; i++) {
      for (auto j=1; j<Nout1-1; j++) {
         if(dir==0) { // X-Flux

            assert(nXg >= 2);

            // calculate flux correction for right going wave
	    //
	    //ratioR = (fin(i,j) - fin(i-1,j))/(fin(i+1,j) - fin(i,j));
	    ratioR = (FlinR(i,j) - FlinR(i-1,j))/(FlinR(i+1,j) - FlinR(i,j));
            FlLimR(i,j) = (ratioR+abs(ratioR))/(1.0+abs(ratioR));
            
	    // calculate flux correction for left going wave
	    //
	    //ratioL = (fin(i+2,j) - fin(i+1,j))/(fin(i+1,j) - fin(i,j));
	    ratioL = (FlinL(i+2,j) - FlinL(i+1,j))/(FlinL(i+1,j) - FlinL(i,j));
            FlLimL(i,j) = (ratioL+abs(ratioL))/(1.0+abs(ratioL));

	 }
	 if(dir==1) { // Z-flux
         
            assert(nZg >= 2);

	    // calculate flux correction for right going wave
	    //
	    //ratioR = (fin(i,j) - fin(i,j-1))/(fin(i,j+1) - fin(i,j));
	    ratioR = (FlinR(i,j) - FlinR(i,j-1))/(FlinR(i,j+1) - FlinR(i,j));
            FlLimR(i,j) = (ratioR+abs(ratioR))/(1.0+abs(ratioR));
            
	    // calculate flux correction for left going wave
	    //
	    //ratioL = (fin(i,j+2) - fin(i,j+1))/(fin(i,j+1) - fin(i,j));
	    ratioL = (FlinL(i,j+2) - FlinL(i,j+1))/(FlinL(i,j+1) - FlinL(i,j));
            FlLimL(i,j) = (ratioL+abs(ratioL))/(1.0+abs(ratioL));

	 }
	   
      }
   }
   
   //  compute total flux
   //
   if(order==1) {
      FloutR = FluxR1st;
      FloutL = FluxL1st;
   }
   else {
      FloutR = FluxR1st + FlLimR*(FluxRC2 - FluxR1st);
      FloutL = FluxL1st + FlLimL*(FluxLC2 - FluxL1st);
   }

   Flout = FloutR+FloutL;

} // end TVD flux new calculation


void domainGrid::computeFluxTVDsimple(matrix2D<double> &Flout, 
                                matrix2D<double> &FloutL,
				matrix2D<double> &FloutR, 
                                const matrix2D<double> &Flin,
                                const matrix2D<double> &upC,
                                const matrix2D<double> &fin,
                                const string& METHOD,
				const int dir) const {

   // similar to computeFluxTVD, but uses METHOD to compute
   // left and right going fluxes rather than MUSCL
   // Note that it may not be TVD scheme now

   // this function interpolates cell-center flux Flin
   // to Flout, defined at cell edges, using TVD scheme 
   // FloutL and FloutR are left and right going fluxes
   // upC is the local maximum characteristic speed
   // and fin is the function being fluxed
   // dir is the directoin (0 for X, 1 for Z)
   // order is scheme order (1 or 2)
   
   
   const int Nout0 = Flout.size0();
   const int Nout1 = Flout.size1();
   const int Nin0 = Flin.size0();
   const int Nin1 = Flin.size1();

   if(dir==0) {
      assert(Nout0 == nXce && Nin0 == nXcc);
   }
   else if(dir==1) {
      assert(Nout1 == nZce && Nin1 == nZcc);
   }
   else {
      cout << "dir in call to computeFluxTVD most be 0 or 1" << endl;
      exit (EXIT_FAILURE);
   }
   assert(upC.size0() == Flin.size0());
   assert(fin.size0() == Flin.size0());
   assert(upC.size1() == Flin.size1());
   assert(fin.size1() == Flin.size1());

   //  define some additional matrix values for calculation
   //
   matrix2D<double> FlinR(Nin0,Nin1,0.0);
   matrix2D<double> FlinL(Nin0,Nin1,0.0);

   FlinR = 0.5*(Flin + upC*fin);
   FlinL = 0.5*(Flin - upC*fin);
  
   InterpToCellEdges(FloutR,FlinR,upC,METHOD,dir);
   InterpToCellEdges(FloutL,FlinL,-upC,METHOD,dir);
  
   Flout = FloutR+FloutL;

}  // end TVD flux simple calculation


void domainGrid::computeFluxTVD(matrix2D<double> &Flout, 
                                matrix2D<double> &FloutL,
				matrix2D<double> &FloutR, 
                                matrix2D<double> &FlLimL, 
                                matrix2D<double> &FlLimR, 
                                const matrix2D<double> &Flin,
                                const matrix2D<double> &upC,
                                const matrix2D<double> &fin,
			        const string &LIMITER,
			 	const int dir, 
				const int order) const {

   // this function interpolates cell-center flux Flin
   // to Flout, defined at cell edges, using TVD scheme 
   // FloutL and FloutR are left and right going fluxes
   // upC is the local maximum characteristic speed
   // and fin is the function being fluxed
   // dir is the directoin (0 for X, 1 for Z)
   // order is scheme order (1 or 2)
   
   
   const int Nout0 = Flout.size0();
   const int Nout1 = Flout.size1();
   const int Nin0 = Flin.size0();
   const int Nin1 = Flin.size1();

   if(dir==0) {
      assert(Nout0 == nXce && Nin0 == nXcc);
   }
   else if(dir==1) {
      assert(Nout1 == nZce && Nin1 == nZcc);
   }
   else {
      cout << "dir in call to computeFluxTVD most be 0 or 1" << endl;
      exit (EXIT_FAILURE);
   }
   assert(upC.size0() == Flin.size0());
   assert(fin.size0() == Flin.size0());
   assert(upC.size1() == Flin.size1());
   assert(fin.size1() == Flin.size1());

  
   //  define some additional matrix values for calculation
   //
   matrix2D<double> FluxL1st(Nout0,Nout1,0.0); 
   matrix2D<double> FluxR1st(Nout0,Nout1,0.0); 
   //matrix2D<double> FluxLC2(Nout0,Nout1,0.0); 
   //matrix2D<double> FluxRC2(Nout0,Nout1,0.0); 
   double DeltaFluxRL, DeltaFluxRR;
   double DeltaFluxLL, DeltaFluxLR;
   matrix2D<double> FlinR(Nin0,Nin1,0.0);
   matrix2D<double> FlinL(Nin0,Nin1,0.0);

   FlinR = 0.5*(Flin + upC*fin);
   FlinL = 0.5*(Flin - upC*fin);
   //FlinR = Flin + fin;
   //FlinL = Flin - fin;
   int iup=1, jup=0;
   if(dir==1) iup=0, jup=1;

   //  compute U1 and C2 left and right fluxes
   //
   for (auto i=0; i<Nout0; i++) {
      for (auto j=0; j<Nout1; j++) {
         FluxR1st(i,j) = FlinR(i,j);
         FluxL1st(i,j) = FlinL(i+iup,j+jup);
         //if(dir==0) FluxL1st(i,j) = FlinL(i+1,j);
         //if(dir==1) FluxL1st(i,j) = FlinL(i,j+1);
      }
   }

   //  compute second order left and right flux corrections
   //  using flux limiter
   //
   for (auto i=1; i<Nout0-1; i++) {
      for (auto j=1; j<Nout1-1; j++) {
         if(dir==0) { // X-Flux

            assert(nXg >= 2);

            // calculate flux correction for right going wave
	    //
            DeltaFluxRL = 0.5*(FlinR(i,j)   - FlinR(i-1,j)); // B2 scheme
            DeltaFluxRR = 0.5*(FlinR(i+1,j) - FlinR(i,j));   // C2 scheme

	    // calculate flux correction for left going wave
	    //
	    DeltaFluxLL = -0.5*(FlinL(i+1,j) - FlinL(i,j));    // C2 scheme
            DeltaFluxLR = -0.5*(FlinL(i+2,j) - FlinL(i+1,j));  // F2 scheme

            //  calculate 2nd order flux correction using limiter
	    //
	    if(LIMITER=="minmod") {
               FlLimR(i,j) = minmod(DeltaFluxRL,DeltaFluxRR);
               FlLimL(i,j) = minmod(DeltaFluxLL,DeltaFluxLR);
            } else if (LIMITER=="superbee") {
	       FlLimR(i,j) = superbee(DeltaFluxRL,DeltaFluxRR);
	       FlLimL(i,j) = superbee(DeltaFluxLL,DeltaFluxLR);
            } else {
	       FlLimR(i,j) = vanleer(DeltaFluxRL,DeltaFluxRR);
	       FlLimL(i,j) = vanleer(DeltaFluxLL,DeltaFluxLR);
            }

	 }
	 if(dir==1) { // Z-flux
         
            assert(nZg >= 2);

	    // calculate flux correction for right going wave
	    //
	    DeltaFluxRL = 0.5*(FlinR(i,j)   - FlinR(i,j-1)); // B2 scheme 
            DeltaFluxRR = 0.5*(FlinR(i,j+1) - FlinR(i,j));   // C2 scheme
            
	    // calculate flux correction for left going wave
	    //
	    DeltaFluxLL = -0.5*(FlinL(i,j+1) - FlinL(i,j));    // C2 scheme
            DeltaFluxLR = -0.5*(FlinL(i,j+2) - FlinL(i,j+1));  // F2 scheme

            //  calculate 2nd order flux correction using limiter
	    //
	    if(LIMITER=="minmod") {
               FlLimR(i,j) = minmod(DeltaFluxRL,DeltaFluxRR);
               FlLimL(i,j) = minmod(DeltaFluxLL,DeltaFluxLR);
            } else if (LIMITER=="superbee") {
	       FlLimR(i,j) = superbee(DeltaFluxRL,DeltaFluxRR);
	       FlLimL(i,j) = superbee(DeltaFluxLL,DeltaFluxLR);
            } else {
	       FlLimR(i,j) = vanleer(DeltaFluxRL,DeltaFluxRR);
	       FlLimL(i,j) = vanleer(DeltaFluxLL,DeltaFluxLR);
            }


	 }
	   
      }
   }

   //  compute total flux
   //
   if(order==1) {
      FloutR = FluxR1st;
      FloutL = FluxL1st;
   }
   else {
      FloutR = FluxR1st + FlLimR;
      FloutL = FluxL1st + FlLimL;
   }
   Flout = FloutR+FloutL;


} // end TVD flux calculation

