/***
 *
 * matrix2D class header file
 *
 * see https://www.quantstart.com/articles/Matrix-Classes-in-C-The-Header-File
 *
***/

#ifndef matrix2D_h
#define matrix2D_h

//#include "vectorMath.h" 
#include <vector>
#include <mpi.h>
#include <math.h>
#include <cmath>

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

   // math operators for matrix with matrix
   //
   matrix2D<T> operator+(const matrix2D<T>& rhs);
   matrix2D<T>& operator+=(const matrix2D<T>& rhs);
   matrix2D<T> operator-(const matrix2D<T>& rhs);
   matrix2D<T>& operator-=(const matrix2D<T>& rhs);
   matrix2D<T> operator*(const matrix2D<T>& rhs);
   matrix2D<T>& operator*=(const matrix2D<T>& rhs);
   matrix2D<T> operator/(const matrix2D<T>& rhs);
   matrix2D<T>& operator/=(const matrix2D<T>& rhs);

   // math operators for matrix with scalar
   //
   matrix2D<T> operator+(const T& rhs);
   matrix2D<T>& operator+=(const T& rhs);
   matrix2D<T> operator-(const T& rhs);
   matrix2D<T>& operator-=(const T& rhs);
   matrix2D<T> operator*(const T& rhs);
   matrix2D<T>& operator*=(const T& rhs);
   matrix2D<T> operator/(const T& rhs);
   matrix2D<T>& operator/=(const T& rhs);
   //matrix2D<T> operator*(const T& thisConst, const matrix2D<T>& thisMat);


   // initialize matrix2D
   //
   void initialize(const int, const int, const T& C0);


   // get the number of nX and nZ points
   //
   unsigned size0() const;
   unsigned size1() const;


   // get the individual elements
   //
   T& operator()(const unsigned iX, const unsigned jZ);
   const T operator()(const unsigned iX, const unsigned jZ) const;

};


// have to define non-member functions for lhs operations
//

// add scalar to matrix on the lhs
//
template <typename T>
matrix2D<T> operator+(const T& thisConst, const matrix2D<T>& thisMat)
{
  matrix2D<T> resultMat(thisMat);
  resultMat = resultMat+thisConst;

  return resultMat;
}

template <typename T>
matrix2D<T> operator+(const matrix2D<T>& matA, const matrix2D<T>& matB)
{
  matrix2D<T> resultMat(matA);
  resultMat += matB;

  return resultMat;
}

template <typename T>
matrix2D<T> operator-(const matrix2D<T>& matA, const matrix2D<T>& matB)
{
  matrix2D<T> resultMat(matA);
  resultMat -= matB;

  return resultMat;
}

template <typename T>
matrix2D<T> operator-(const matrix2D<T>& matA)
{
  matrix2D<T> resultMat(matA);
  resultMat *= -1.0;

  return resultMat;
}

template <typename T>
matrix2D<T> operator*(const matrix2D<T>& matA, const matrix2D<T>& matB)
{
  matrix2D<T> resultMat(matA);
  resultMat *= matB;

  return resultMat;
}

// subtract matrix from scalar on the lhs
//
template <typename T>
matrix2D<T> operator-(const T& thisConst, const matrix2D<T>& thisMat)
{
  matrix2D<T> resultMat(thisMat);
  resultMat = resultMat-thisConst;
  resultMat *= -1.0;

  return resultMat;
}

// multiply matrix by a scalar on the lhs
//
template <typename T>
matrix2D<T> operator*(const T& thisConst, const matrix2D<T>& thisMat)
{
  matrix2D<T> resultMat(thisMat);
  resultMat = resultMat*thisConst;

  return resultMat;
}

// divide scalar by matrix
//
template <typename T>
matrix2D<T> operator/(const T& thisConst, const matrix2D<T>& thisMat)
{
  matrix2D<T> resultMat(thisMat);
  matrix2D<T> thisConstMat(thisMat.size0(),thisMat.size1(),thisConst);
  resultMat = thisConstMat/resultMat;

  return resultMat;
}

// basic math operators
//
template <typename T>
matrix2D<T> exp(const matrix2D<T> &thisMat)
{
   matrix2D<T> resultMat(thisMat);
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();
   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
         resultMat(i,j) = exp(resultMat(i,j));
      }
   }
   return resultMat;
}
template <typename T>
matrix2D<T> sqrt(const matrix2D<T> &thisMat)
{
   matrix2D<T> resultMat(thisMat);
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();
   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
         resultMat(i,j) = sqrt(resultMat(i,j));
      }
   }
   return resultMat;
}
template <typename T>
matrix2D<T> tanh(const matrix2D<T> &thisMat)
{
   matrix2D<T> resultMat(thisMat);
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();
   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
         resultMat(i,j) = tanh(resultMat(i,j));
      }
   }
   return resultMat;
}
template <typename T>
matrix2D<T> log(const matrix2D<T> &thisMat)
{
   matrix2D<T> resultMat(thisMat);
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();
   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
         resultMat(i,j) = log(resultMat(i,j));
      }
   }
   return resultMat;
}
template <typename T>
matrix2D<T> cos(const matrix2D<T> &thisMat)
{
   matrix2D<T> resultMat(thisMat);
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();
   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
         resultMat(i,j) = cos(resultMat(i,j));
      }
   }
   return resultMat;
}
template <typename T>
matrix2D<T> sin(const matrix2D<T> &thisMat)
{
   matrix2D<T> resultMat(thisMat);
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();
   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
         resultMat(i,j) = sin(resultMat(i,j));
      }
   }
   return resultMat;
}
template <typename T>
matrix2D<T> pow(const matrix2D<T> &thisMat, const T exponent)
{
   matrix2D<T> resultMat(thisMat);
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();
   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
         resultMat(i,j) = pow(resultMat(i,j),exponent);
      }
   }
   return resultMat;
}
template <typename T>
matrix2D<T> abs(const matrix2D<T> &thisMat)
{
   matrix2D<T> resultMat(thisMat);
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();
   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
         resultMat(i,j) = abs(resultMat(i,j));
      }
   }
   return resultMat;
}
template <typename T>
T min(const matrix2D<T> &thisMat)
{
   T result, thisval, localresult;
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();

   localresult = thisMat(0,0);
   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
	 thisval = thisMat(i,j);
	 if(thisval<localresult) localresult = thisval;
      }
   }

   MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  
   return result;
}
template <typename T>
T max(const matrix2D<T> &thisMat)
{
   T result, thisval, localresult;
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();

   localresult = thisMat(0,0);
   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
	 thisval = thisMat(i,j);
	 if(thisval>localresult) localresult = thisval;
      }
   }

   MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  
   return result;
}
template <typename T>
matrix2D<T> max(const matrix2D<T> &thisMat, const T a0)
{
   matrix2D<T> result(thisMat);
   T  thisval;
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();

   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
	 result(i,j) = max(thisMat(i,j),a0);
	 //thisval = thisMat(i,j);
	 //if(thisval<a0) result(i,j) = a0;
      }
   }

   return result;
}
template <typename T>
matrix2D<T> min(const matrix2D<T> &thisMat, const T a0)
{
   matrix2D<T> result(thisMat);
   T  thisval;
   const int thisNx = thisMat.size0();
   const int thisNz = thisMat.size1();

   for (auto i=0; i<thisNx; i++) {
      for (auto j=0; j<thisNz; j++) {
	 result(i,j) = min(thisMat(i,j),a0);
	 //thisval = thisMat(i,j);
	 //if(thisval<a0) result(i,j) = a0;
      }
   }

   return result;
}



#include "matrix2D.cpp"


#endif
