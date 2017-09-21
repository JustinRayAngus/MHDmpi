/***
 *
 * template for doing vector math oporations
 *  such as +, -, *, and /
 *
***/

#ifndef vectorMath_h
#define vectorMath_h

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

using namespace std;


// vector a + vector b
//
template <typename T>
vector<T> operator+(const vector<T>& a, const vector<T> &b)
{
   assert(a.size() == b.size());

   vector<T> result;
   result.reserve(a.size());
   transform(a.begin(), a.end(), b.begin(),
             back_inserter(result), plus<T>());

   return result;
}

// vector a + double b
//
template <typename T>
vector<T> operator+(const vector<T>& a, const T &b)
{

   vector<T> result;
   vector<T> bvec;
   result.reserve(a.size());
   bvec.assign(a.size(),b);
   transform(a.begin(), a.end(), bvec.begin(),
             back_inserter(result), plus<T>());

   return result;
} 

//double b + vector a
//
template <typename T>
vector<T> operator+(const T &b, const vector<T>& a)
{

   vector<T> result;
   result.reserve(a.size());
   result = a + b;

   return result;
}

// vector a - vector b
//
template <typename T>
vector<T> operator-(const vector<T> &a, const vector<T> &b)
{
   assert(a.size() == b.size());

   vector<T> result;
   result.reserve(a.size());
   transform(a.begin(), a.end(), b.begin(),
             back_inserter(result), minus<T>());

   return result;
}

// vector a - double b
//
template <typename T>
vector<T> operator-(const vector<T> &a, const T &b)
{

   vector<T> result;
   vector<T> bvec;
   result.reserve(a.size());
   bvec.assign(a.size(),b);
   transform(a.begin(), a.end(), bvec.begin(),
             back_inserter(result), minus<T>());

   return result;
}

// double b - vector a
//
template <typename T>
vector<T> operator-(const T &b, const vector<T> &a)
{

   vector<T> result;
   vector<T> bvec;
   result.reserve(a.size());
   bvec.assign(a.size(),b);
   result = bvec - a;

   return result;
}

// - vector a
//
template <typename T>
vector<T> operator-(const vector<T> &a)
{

   vector<T> result;
   result.reserve(a.size());
   result = -1.0*a;

   return result;
}

// vector a * vector b (dot product)
//
template <typename T>
vector<T> operator*(const vector<T>& a, const vector<T> &b)
{
   assert(a.size() == b.size());

   vector<T> result;
   result.reserve(a.size());
   transform(a.begin(), a.end(), b.begin(),
             back_inserter(result), multiplies<T>());

   return result;
}

// vector a * double b
//
template <typename T>
vector<T> operator*(const vector<T>& a, const T &b)
{
   //assert(a.size() == b.size());

   vector<T> result;
   result = a;
   transform(result.begin(), result.end(), result.begin(),
             bind1st(multiplies<T>(),b));

   return result;
}

// double b * vector a
template <typename T>
vector<T> operator*(const T &b, const vector<T>& a)
{
   vector<T> result;
   result.reserve(a.size());
   result = a*b;

   //transform(result.begin(), result.end(), result.begin(),
   //          bind1st(multiplies<T>(),b));

   return result;
}

// vector a / vector b
//
template <typename T>
vector<T> operator/(const vector<T>& a, const vector<T> &b)
{
   assert(a.size() == b.size());

   vector<T> result;
   result.reserve(a.size());
   transform(a.begin(), a.end(), b.begin(),
             back_inserter(result), divides<T>());

   return result;
}

// vector a / double b
//
template <typename T>
vector<T> operator/(const vector<T>& a, const T &b)
{

   vector<T> result;
   result.reserve(a.size());
   result = a;
   transform(result.begin(), result.end(), result.begin(),
             bind1st(multiplies<T>(),1.0/b));

   return result;
}

// double b / vector a
//
template <typename T>
vector<T> operator/(const T &b, const vector<T> &a)
{

   const int imax = a.size();
   vector<T> result;
   result.resize(a.size());
   for (auto i=1; i<imax; i++) {
      result.at(i) = b/a.at(i);
   }

   return result;
}


////////////////////////////////////////////////////////////////
//
// functions for other common math operators
// exp(), tanh, log(), cos(), sin() ...
//

vector<double> exp(const vector<double> &fin);

vector<double> tanh(const vector<double> &fin);

vector<double> log(const vector<double> &fin);

vector<double> cos(const vector<double> &fin);

vector<double> sin(const vector<double> &fin);

vector<double> abs(const vector<double> &fin);

double min(const vector<double> &fin);

double max(const vector<double> &fin);


#endif

