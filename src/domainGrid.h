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
  int nX, nXsub, nXg, nXcc, nXce, nXce2;
  int nZ=1, nZsub=1, nZg=1, nZcc=1, nZce=2, nZce2=2;
  int numProcs, procID;
  double Xmax, Xmin, dX;
  double Zmax=1, Zmin=0, dZ=1;
  vector<double> Xcc, Xce, Xce2; // X-grid at cell-center and at cell-edge 
  vector<double> Zcc, Zce, Zce2; // Z-grid at cell-center and at cell-edge 

  static domainGrid* mesh; // pointer to a domainGrid instance

  void initialize(const Json::Value&);
  void setInitialProfile(vector<double>&, const Json::Value&) const;
  void setInitialProfile(matrix2D<double>&, const Json::Value&) const;
  void setInitialProfileArbDir(vector<double>&, const vector<double>&,
		         const double, const double,
		         const double, const double,
		         const double, const double,	
		         const string&) const;
  
  void communicate(vector<double>&) const;
  void communicate(matrix2D<double>&) const;

  void InterpToCellEdges(vector<double>&, const vector<double>&,
                         const vector<double>&, const string&) const;
  void InterpToCellEdges(matrix2D<double>&, 
		   const matrix2D<double>&,
                   const matrix2D<double>&, 
		   const string&, const int) const;
  
  void InterpCellToEdges( matrix2D<double>&, 
		   const  matrix2D<double>&,
                   const  matrix2D<double>&, 
		   const  string&, 
                   const  int ) const;

  void InterpEdgesToEdges( matrix2D<double>&, 
                     const matrix2D<double>&  ) const;
 
  void InterpNodesToEdges( matrix2D<double>&, 
                     const matrix2D<double>&,
                     const int                ) const; 
  
  
  void InterpToCellCenter(vector<double>&, const vector<double>&) const;
  void InterpToCellCenter(matrix2D<double>&, 
		          const matrix2D<double>&) const;
  
  void computeFluxTVD(vector<double>&, vector<double>&, vector<double>&,
                      vector<double>&, vector<double>&,
                      const vector<double>&, const vector<double>&,
                      const vector<double>&,
	              const string&,
		      const int) const;
  void computeFluxTVD(matrix2D<double>&, matrix2D<double>&, 
                      matrix2D<double>&, matrix2D<double>&,
                      matrix2D<double>&,
                      const matrix2D<double>&, 
		      const matrix2D<double>&,
		      const matrix2D<double>&,
	              const string&,
		      const int,
		      const int) const;
  void computeFluxTVDnew(matrix2D<double>&, matrix2D<double>&, 
                      matrix2D<double>&, matrix2D<double>&,
                      matrix2D<double>&,
                      const matrix2D<double>&, 
		      const matrix2D<double>&,
		      const matrix2D<double>&,
		      const int,
		      const int) const;
  void computeFluxTVDsimple(matrix2D<double>&, 
                      matrix2D<double> &,
		      matrix2D<double> &, 
                      const matrix2D<double>&, 
		      const matrix2D<double>&,
		      const matrix2D<double>&,
                      const string&,
		      const int) const;
  void computeFluxTVDsimple(vector<double>&, 
                      vector<double> &,
		      vector<double> &, 
                      const vector<double>&, 
		      const vector<double>&,
		      const vector<double>&,
		      const string&) const;
  
  void DDX(vector<double>&, const vector<double>&) const;
  void DDX(matrix2D<double>&, const matrix2D<double>&) const;
  void DDZ(matrix2D<double>&, const matrix2D<double>&) const;
  void D2DZ2(matrix2D<double>&, const matrix2D<double>&) const;

  void setXminFluxBC( vector<double>&, const double, const double ) const;
  void setXmaxFluxBC( vector<double>&, const double, const double ) const;
  void setXminBoundary( vector<double>&, const double, const double ) const;
  void setXmaxBoundary( vector<double>&, const double, const double ) const;
  void setXminBoundary_J( vector<double>&, const double, const double ) const;
  void setXmaxBoundary_J( vector<double>&, const double, const double ) const;

  void setXminFluxBC( matrix2D<double>&, const double, const double ) const;
  void setXmaxFluxBC( matrix2D<double>&, const double, const double ) const;
  void setXminFluxBC( matrix2D<double>&, const vector<double>& ) const;
  void setXmaxFluxBC( matrix2D<double>&, const vector<double>& ) const;
  void setXminBoundary( matrix2D<double>&, const vector<double>& ) const;
  void setXmaxBoundary( matrix2D<double>&, const vector<double>& ) const;
  void setXminBoundary( matrix2D<double>&, const double, const double ) const;
  void setXmaxBoundary( matrix2D<double>&, const double, const double ) const;
  void setXminBoundaryExtrap( matrix2D<double>& ) const;
  void setXmaxBoundaryExtrap( matrix2D<double>& ) const;
  void setXminBoundary_J( matrix2D<double>&, const double, const double ) const;
  void setXmaxBoundary_J( matrix2D<double>&, const double, const double ) const;
  void setXminBoundary_J( matrix2D<double>&, const vector<double>& ) const;
 
  void setZboundaryPeriodic( matrix2D<double>& ) const;
  void setZboundaryInlet( matrix2D<double>&, 
                    const int,
                    const double, 
                    const double ) const;
  void setZminBoundary( matrix2D<double>&, const double, const double ) const;
  void setZmaxBoundary( matrix2D<double>&, const double, const double ) const;
  void setZminFluxBC( matrix2D<double>&, const double ) const;
  void setZmaxFluxBC( matrix2D<double>&, const double ) const;
  void setZminFluxBC( matrix2D<double>&, const matrix2D<double>& ) const;
  void setZmaxFluxBC( matrix2D<double>&, const matrix2D<double>& ) const;

//domainGrid() {};

};

//domainGrid* domainGrid::mesh = NULL;


#endif
