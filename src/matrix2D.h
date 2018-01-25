/***
 *
 * matrix2D class header file
 *
 * see https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File
 *
***/

#ifndef matrix2D_h
#define matrix2D_h

#include <vector>

using namespace std;


template <typename T> 
class matrix2D 
{
private: 
   vector<T> VecXZ; // vector representation of matrix values
   unsigned nX;
   unsigned nZ;


public:
   matrix2D();
   matrix2D(unsigned thisnX, unsigned thisnZ, const T& C0);
   matrix2D(const matrix2D<T>& rhs);
   virtual ~matrix2D();
 
   // operator overload for standard math matrix operations
   //
   matrix2D<T>& operator=(const matrix2D<T>& rhs);


   // math operators for matrix
   //
   matrix2D<T> operator+(const matrix2D<T>& rhs);
   //matrix2D<T>& operator+=(const matrix2D<T>& rhs);
   matrix2D<T> operator-(const matrix2D<T>& rhs);
   matrix2D<T> operator*(const matrix2D<T>& rhs);
   matrix2D<T> operator/(const matrix2D<T>& rhs);


   // initialize matrix2D
   //
   void initialize(const int, const int, const T& C0);


   // get the number of nX and nZ points
   //
   unsigned size0() const;
   unsigned size1() const;


   // get the individual elements
   //
   T& operator()(const unsigned& iX, const unsigned& jZ);
   const T& operator()(const unsigned& iX, const unsigned& jZ) const;

};

#include "matrix2D.cpp"

#endif
