/***
 *
 * domainGrid class header file
 *
***/

#ifndef domainGrid_h
#define domainGrid_h


#include "json/json.h"

using namespace std;

//vector<double> DDX(const vector<double> &Fin, const double &dx);
// above funtion doesn't use anything from grid and doesn't
// belong defined here


class domainGrid
{

public:
  int nX, nXsub, nXg, nXcc, nXce;
  int numProcs, procID;
  double Xmax, Xmin, dX;
  vector<double> Xcc, Xce; // spatial grid at cell-center and at cell-edge 

  static domainGrid* mesh; // pointer to a domainGrid instance

  void initialize(const Json::Value&);
  void setInitialProfile(vector<double>&, const Json::Value&) const;
  void communicate(vector<double>&) const;
  void InterpToCellEdges(vector<double>&, const vector<double>&,
                         const vector<double>&, const string&) const;
  void InterpToCellCenter(vector<double>&, const vector<double>&) const;
  void computeFluxTVD(vector<double>&, vector<double>&, vector<double>&,
                      vector<double>&, vector<double>&,
                      const vector<double>&, const vector<double>&,
                      const vector<double>&,
		      const int) const;
  void DDX(vector<double>&, const vector<double>&) const;

//domainGrid() {};

};

//domainGrid* domainGrid::mesh = NULL;


#endif
