/***
 * 
 * domain grid class
 *
***/

#ifndef domainGrid_h
#define domainGrid_h

using namespace std;

class domainGrid
{

public:
  int nX, nXsub, nXg;
  int numProcs, procID;
  double Xmax, Xmin, dX;
  vector<double> Xcc, Xce; // spatial grid at cell-center and at cell-edge 

  void initialize(const Json::Value&);

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
   Xce.assign(nXsub+2*nXg-1,0.0);
   double offset = procID*(Xmax-Xmin)/numProcs-(0.5+nXg-1.0)*dX;
   const int nMax = Xce.size();
   for (auto n=0; n<nMax; n++) {
      Xcc.at(n) = Xmin + offset + n*dX;
      Xce.at(n) = Xmin + offset + n*dX + 0.5*dX;
   }
   Xcc.at(nMax) = Xcc.at(nMax-1)+dX;
   //Xcc[nX/numProcs+1] = Xmin + offset + (nXsub+1)*dX;

}


#endif
