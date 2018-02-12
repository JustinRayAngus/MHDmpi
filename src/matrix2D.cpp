/***
 *
 * matrix2D class source file
 *
 *https://stackoverflow.com/questions/4421706/what-are-the-basic-rules-and-idioms-for-operator-overloading/4421719#4421719
 *
***/

#ifndef matrix2D_cpp
#define matrix2D_cpp


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

#include "vectorMath.h"
//#include "domainGrid.h"
#include "matrix2D.h"

using namespace std;


template<typename T>
matrix2D<T>::matrix2D() 
{
   // do nothing for default constructor
}

template<typename T>
matrix2D<T>::matrix2D(unsigned thisnX, unsigned thisnZ, const T& C0) 
{

   nX = thisnX;
   nZ = thisnZ;
   const unsigned nVec = nX*nZ;
   
   VecXZ.assign(nVec,C0);   
   //VecXZ.reserve(nVec);   

}

template<typename T>
matrix2D<T>::matrix2D(const matrix2D<T>& rhs)
{
   VecXZ = rhs.VecXZ;
   nX = rhs.nX;
   nZ = rhs.nZ;
}

template<typename T>
matrix2D<T>::~matrix2D() {}


template<typename T>
void matrix2D<T>::initialize(const int thisnX, const int thisnZ, const T& C0)
{
   //MPI_Comm_rank (MPI_COMM_WORLD, &procID);
   //MPI_Comm_size (MPI_COMM_WORLD, &numProcs);

   nX = thisnX;
   nZ = thisnZ;
   const int nVec = nX*nZ;
   //cout << "VecXZ.size() = " << VecXZ.size() << endl;
   VecXZ.assign(nVec,C0);
   //cout << "VecXZ.size() = " << VecXZ.size() << endl;

}

template<typename T>
T& matrix2D<T>::operator()(const unsigned iX, const unsigned jZ)
{
   return VecXZ.at(iX*nZ + jZ);
}

template<typename T>
const T matrix2D<T>::operator()(const unsigned iX, const unsigned jZ) const
{
   return VecXZ.at(iX*nZ + jZ);
}


template<typename T>
matrix2D<T>& matrix2D<T>::operator=(const matrix2D& rhs)
{
   unsigned rhs_nX = rhs.size0();
   unsigned rhs_nZ = rhs.size1();
   if(rhs_nX != nX || rhs_nZ != nZ ) {
      cout << "tring to set matrix2D of different sizes equal to each other" << endl;
      exit(EXIT_FAILURE);
   }
   else {	   
      //unsigned rhs_nVec = rhs_nX*rhs_nZ;
      //VecXZ.assign(rhs_nVec,rhs.VecXZ);
      VecXZ = rhs.VecXZ;
      return *this;
      //return &rhs;
   }
   //return *this;
}

// update matrix by adding another to it using +=
//
template<typename T>
matrix2D<T>& matrix2D<T>::operator+=(const matrix2D& rhs)
{
   matrix2D<T> result = (*this) + rhs;
   (*this) = result;   
   return *this;
   
   //transform(this->VecXZ.begin(), this->VecXZ.end(), rhs.VecXZ.begin(),
   //	     back_inserter(resultVec),plus<T>());
}

// add two matricies together
//
template<typename T>
matrix2D<T> matrix2D<T>::operator+(const matrix2D& rhs)
{
   vector<T> resultVec;
   matrix2D<T> resultMat;

   unsigned rhs_nX = rhs.size0();
   unsigned rhs_nZ = rhs.size1();
   if(rhs_nX != nX || rhs_nZ != nZ ) {
      cout << "tring to add two matrix2D of different sizes" << endl;
      exit(EXIT_FAILURE);
   }
   else {	  
      resultMat.nX = rhs_nX;
      resultMat.nZ = rhs_nZ;
      resultVec = VecXZ+rhs.VecXZ;
      //resultMat.VecXZ = VecXZ+rhs.VecXZ;
      //transform(VecXZ.begin(), VecXZ.end(), rhs.VecXZ.begin(),
       // 	back_inserter(resultVec),plus<T>());
      resultMat.VecXZ = resultVec;
      return resultMat;
   }
   //return *this;
}

// update matrix by subtracting another from it using -=
//
template<typename T>
matrix2D<T>& matrix2D<T>::operator-=(const matrix2D& rhs)
{
   matrix2D<T> result = (*this) - rhs;
   (*this) = result;   
   return *this;
}


// subtract two matricies
//
template<typename T>
matrix2D<T> matrix2D<T>::operator-(const matrix2D& rhs)
{
   vector<T> resultVec;
   matrix2D<T> resultMat;

   unsigned rhs_nX = rhs.size0();
   unsigned rhs_nZ = rhs.size1();
   if(rhs_nX != nX || rhs_nZ != nZ ) {
      cout << "tring to subtract two matrix2D of different sizes" << endl;
      exit(EXIT_FAILURE);
   }
   else {	  
      resultMat.nX = rhs_nX;
      resultMat.nZ = rhs_nZ;
      resultVec = VecXZ-rhs.VecXZ;
      //resultMat.VecXZ = VecXZ-rhs.VecXZ;
      //transform(VecXZ.begin(), VecXZ.end(), rhs.VecXZ.begin(),
      //  	back_inserter(resultVec),minus<T>());
      resultMat.VecXZ = resultVec;
      return resultMat;
   }
   //return *this;
}


// multiply two matricies
//
template<typename T>
matrix2D<T> matrix2D<T>::operator*(const matrix2D& rhs)
{
   vector<T> resultVec;
   matrix2D<T> resultMat;

   unsigned rhs_nX = rhs.size0();
   unsigned rhs_nZ = rhs.size1();
   if(rhs_nX != nX || rhs_nZ != nZ ) {
      cout << "tring to multiply two matrix2D of different sizes" << endl;
      exit(EXIT_FAILURE);
   }
   else {	  
      resultMat.nX = rhs_nX;
      resultMat.nZ = rhs_nZ;
      resultVec = VecXZ*rhs.VecXZ;
      //resultMat.VecXZ = VecXZ*rhs.VecXZ;
      //transform(VecXZ.begin(), VecXZ.end(), rhs.VecXZ.begin(),
      //          back_inserter(resultVec),multiplies<T>());
      resultMat.VecXZ = resultVec;
      return resultMat;
   }
   //return *this;
}


// multiply matrix by another matrix using *=
//
template<typename T>
matrix2D<T>& matrix2D<T>::operator*=(const matrix2D& rhs)
{
   matrix2D<T> result = (*this) * rhs;
   (*this) = result;   
   return *this;
}

// divide two matricies
//
template<typename T>
matrix2D<T> matrix2D<T>::operator/(const matrix2D& rhs)
{
   vector<T> resultVec;
   matrix2D<T> resultMat;

   unsigned rhs_nX = rhs.size0();
   unsigned rhs_nZ = rhs.size1();
   if(rhs_nX != nX || rhs_nZ != nZ ) {
      cout << "tring to divide two matrix2D of different sizes" << endl;
      exit(EXIT_FAILURE);
   }
   else {	  
      resultMat.nX = rhs_nX;
      resultMat.nZ = rhs_nZ;
      resultVec = VecXZ/rhs.VecXZ;
      //resultMat.VecXZ = VecXZ/rhs.VecXZ;
      resultMat.VecXZ = resultVec;
      //transform(VecXZ.begin(), VecXZ.end(), rhs.VecXZ.begin(),
      //	  back_inserter(resultMat.VecXZ),divides<T>());
      return resultMat;
   }
   //return *this;
}


// update matrix by dividing it by another using /=
//
template<typename T>
matrix2D<T>& matrix2D<T>::operator/=(const matrix2D& rhs)
{
   matrix2D<T> result = (*this) / rhs;
   (*this) = result;   
   return *this;
}

// add a scalar value to a matrix on rhs
//
template<typename T>
matrix2D<T> matrix2D<T>::operator+(const T& rhs)
{
   vector<T> resultVec;
   matrix2D<T> resultMat;

   resultMat.nX = this->size0();
   resultMat.nZ = this->size1();
   resultVec = rhs+VecXZ;
   resultMat.VecXZ = resultVec;
   //resultMat.VecXZ = rhs+VecXZ;
   return resultMat;
}

// add scalar to matrix on rhs
//
template<typename T>
matrix2D<T>& matrix2D<T>::operator+=(const T& rhs)
{
   matrix2D<T> result = (*this) + rhs;
   (*this) = result;
   return *this;
}


template<typename T>
matrix2D<T> matrix2D<T>::operator-(const T& rhs)
{
   vector<T> resultVec;
   matrix2D<T> resultMat;

   resultMat.nX = this->size0();
   resultMat.nZ = this->size1();
   resultVec = VecXZ-rhs;
   resultMat.VecXZ = resultVec;
   return resultMat;
}

// subtract scalar from matrix on rhs
//
template<typename T>
matrix2D<T>& matrix2D<T>::operator-=(const T& rhs)
{
   matrix2D<T> result = (*this) - rhs;
   (*this) = result;
   return *this;
}

// multiply matrix by a scalar on rhs
//
template<typename T>
matrix2D<T> matrix2D<T>::operator*(const T& rhs)
{
   vector<T> resultVec;
   matrix2D<T> resultMat;

   resultMat.nX = this->size0();
   resultMat.nZ = this->size1();
   resultVec = rhs*VecXZ;
   resultMat.VecXZ = resultVec;
   //resultMat.VecXZ = rhs*VecXZ;
   return resultMat;
}

// multiply matrix by a scalar on rhs
//
template<typename T>
matrix2D<T>& matrix2D<T>::operator*=(const T& rhs)
{
   matrix2D<T> result = (*this) * rhs;
   (*this) = result;
   return *this;
}

// divide matrix by scalar on rhs
//
template<typename T>
matrix2D<T> matrix2D<T>::operator/(const T& rhs)
{
   vector<T> resultVec;
   matrix2D<T> resultMat;

   resultMat.nX = this->size0();
   resultMat.nZ = this->size1();
   resultVec = VecXZ/rhs;
   resultMat.VecXZ = resultVec;
   return resultMat;
}

// multiply matrix by a scalar on rhs
//
template<typename T>
matrix2D<T>& matrix2D<T>::operator/=(const T& rhs)
{
   matrix2D<T> result = (*this) / rhs;
   (*this) = result;
   return *this;
}


////


template<typename T>
unsigned matrix2D<T>::size0() const
{
   return this->nX;
   //return nX;
}

template<typename T>
unsigned matrix2D<T>::size1() const
{
   return this->nZ;
   //return nZ;
}

/*
template <typename T>
matrix2D<T> abs(const matrix2D<T> &fin)
{
   matrix2D<T> resultMat(fin);
   vector<T> = resultMat.VecXZ;
   return resultMat;
}
*/

#endif
