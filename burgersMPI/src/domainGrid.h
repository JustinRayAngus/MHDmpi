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
      printf("\nInitializing domain grid ...\n");
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
      if(nXg != nXgVal.asDouble() || nXg < 1 || nXg>5) {
         cout << "ERROR: number of guard cells nXg" << endl;
         cout << "is not set correctly in input file" << endl;
         exit (EXIT_FAILURE);
      }
      Xmin = XminVal.asDouble();
      Xmax = XmaxVal.asDouble();
      dX = (Xmax-Xmin)/nX;
      nXsub = nX/numProcs;
      if(procID==0) {
         double nXsubTest = nX/(double)numProcs;
         if(floor(nXsubTest)==ceil(nXsubTest)) {
            cout << "nXsub = " << nXsub << endl;
         }
         else {
            printf("ERROR: nXsub=nX/numProcs is not an integer!!!!\n");
            exit (EXIT_FAILURE); 
         }
      }
 
      //
      cout << "nX = " << nX << endl;
      cout << "Xmax = " << Xmax << endl;
      cout << "dX = " << dX << endl;
      //for (auto n = 0; n < Xgrid.nX; n++) {
      //   cout << "Xcc[" << n << "] = " << Xgrid.Xcc[n] << endl;
      //   cout << "Xce[" << n << "] = " << Xgrid.Xce[n] << endl;
      //}
   }
   else {
      cout << "value for key \"Xgrid\" is not object type !" << endl;
   }
  
   Xcc.assign(nXsub+2,0.0);
   Xce.assign(nXsub+1,0.0);
   double offset = procID*Xmax/numProcs-0.5*dX;
   for (auto n=0; n<nXsub+1; n++) {
      Xcc[n] = offset + n*dX;
      Xce[n] = offset + n*dX + 0.5*dX;
   }
   Xcc[nX/numProcs+1] = offset + (nXsub+1)*dX;

}


#endif
