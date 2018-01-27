/***
 *
 * domainGrid class header file
 *
***/

#ifndef domainGrid_h
#define domainGrid_h

#include "matrix2D.h"
#include "json/json.h"

using namespace std;

//vector<double> DDX(const vector<double> &Fin, const double &dx);
// above funtion doesn't use anything from grid and doesn't
// belong defined here


class domainGrid
{

public:
  int nX, nXsub, nXg, nXcc, nXce;
  int nZ=1, nZsub=1, nZg=1, nZcc=1, nZce=2;
  int numProcs, procID;
  double Xmax, Xmin, dX;
  double Zmax=1, Zmin=0, dZ=1;
  vector<double> Xcc, Xce; // X-grid at cell-center and at cell-edge 
  vector<double> Zcc, Zce; // Z-grid at cell-center and at cell-edge 

  static domainGrid* mesh; // pointer to a domainGrid instance

  void initialize(const Json::Value&);
  void setInitialProfile(vector<double>&, const Json::Value&) const;
  void setInitialProfile(vector<vector<double>>&, const Json::Value&) const;
  void setInitialProfile(matrix2D<double>&, const Json::Value&) const;
  void setInitialProfileArbDir(vector<double>&, const vector<double>&,
		         const double, const double,
		         const double, const double,	
		         const string&) const;
  
  void communicate(vector<double>&) const;
  void communicate(vector<vector<double>>&) const;
  void communicate(matrix2D<double>&) const;

  void InterpToCellEdges(vector<double>&, const vector<double>&,
                         const vector<double>&, const string&) const;
  void InterpToCellEdges(vector<vector<double>>&, 
		   const vector<vector<double>>&,
                   const vector<vector<double>>&, 
		   const string&, const int) const;
  void InterpToCellEdges(matrix2D<double>&, 
		   const matrix2D<double>&,
                   const matrix2D<double>&, 
		   const string&, const int) const;
  
  
  void InterpToCellCenter(vector<double>&, const vector<double>&) const;
  void InterpToCellCenter(vector<vector<double>>&, 
		          const vector<vector<double>>&) const;
  void InterpToCellCenter(matrix2D<double>&, 
		          const matrix2D<double>&) const;
  
  void computeFluxTVD(vector<double>&, vector<double>&, vector<double>&,
                      vector<double>&, vector<double>&,
                      const vector<double>&, const vector<double>&,
                      const vector<double>&,
		      const int) const;
  void computeFluxTVD(vector<vector<double>>&, vector<vector<double>>&, 
		      vector<vector<double>>&, vector<vector<double>>&,
                      vector<vector<double>>&,
                      const vector<vector<double>>&, 
		      const vector<vector<double>>&,
		      const vector<vector<double>>&,
		      const int,
		      const int) const;
  void computeFluxTVD(matrix2D<double>&, matrix2D<double>&, 
                      matrix2D<double>&, matrix2D<double>&,
                      matrix2D<double>&,
                      const matrix2D<double>&, 
		      const matrix2D<double>&,
		      const matrix2D<double>&,
		      const int,
		      const int) const;
  
  
  void DDX(vector<double>&, const vector<double>&) const;
  void DDX(vector<vector<double>>&, const vector<vector<double>>&) const;
  void DDZ(vector<vector<double>>&, const vector<vector<double>>&) const;
  void DDX(matrix2D<double>&, const matrix2D<double>&) const;
  void DDZ(matrix2D<double>&, const matrix2D<double>&) const;

//domainGrid() {};

};

//domainGrid* domainGrid::mesh = NULL;


#endif
