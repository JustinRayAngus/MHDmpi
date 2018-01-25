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


//////////////////////////////////////////
//
//    vector<vector> operators
//


// matrix a + matrix b
//
template <typename T>
vector<vector<T>> operator+(const vector<vector<T>>& a, const vector<vector<T>>& b)
{
   const int imax = a.size();
   const int jmax = a[0].size();
   assert(a.size() == b.size());
   assert(a[0].size() == b[0].size());
   
   /*
   vector<vector<T>> result;
   vector<T> rvec;
   
   for (auto i=0; i<imax; i++) {
       rvec.reserve(jmax);
       transform(a.at(i).begin(), a.at(i).end(), b.at(i).begin(),
                 back_inserter(rvec), plus<T>());
       if(i==0) {
	   result.resize(1,rvec);
       }
       else { 
           result.push_back(rvec);
       }
       rvec.clear();
   }
   */
   vector<vector<T>> result(imax,vector<T>(jmax));
   
   for (auto i=0; i<imax; i++) {
      for (auto j=0; j<jmax; j++) {
         result[i][j] = a[i][j]+b[i][j];
      }
   }
   
   return result;
}

// matrix a + double b
//
template <typename T>
vector<vector<T>> operator+(const vector<vector<T>>& a, const T& b)
{
   const int imax = a.size();
   const int jmax = a[0].size();

   vector<vector<T>> result(imax,vector<T>(jmax));
   const vector<T>   bvector(jmax,b);
   const vector<vector<T>> bmatrix(imax,bvector);
   result = a + bmatrix;

   return result;
}

// double b + matrix a
//
template <typename T>
vector<vector<T>> operator+(const T& b, const vector<vector<T>>& a)
{
   const int imax = a.size();
   const int jmax = a[0].size();

   vector<vector<T>> result(imax,vector<T>(jmax));
   result = a + b;

   return result;
}

// matrix a - matrix b
//
template <typename T>
vector<vector<T>> operator-(const vector<vector<T>>& a, const vector<vector<T>>& b)
{
   const int imax = a.size();
   const int jmax = a[0].size();
   assert(a.size() == b.size());
   assert(a[0].size() == b[0].size());

   /*
   vector< vector<T> > result;
   vector<T> rvec;
   
   for (auto i=0; i<imax; i++) {
       rvec.reserve(jmax);
       transform(a.at(i).begin(), a.at(i).end(), b.at(i).begin(),
                 back_inserter(rvec), minus<T>());
       if(i==0) {
	   result.resize(1,rvec);
       }
       else { 
           result.push_back(rvec);
       }
       rvec.clear();
   }
   */
   vector<vector<T>> result(imax,vector<T>(jmax));
   
   for (auto i=0; i<imax; i++) {
      for (auto j=0; j<jmax; j++) {
         result[i][j] = a[i][j]-b[i][j];
      }
   }

   return result;
}

// matrix a - double b
//
template <typename T>
vector<vector<T>> operator-(const vector<vector<T>>& a, const T& b)
{
   const int imax = a.size();
   const int jmax = a[0].size();

   vector<vector<T>> result(imax,vector<T>(jmax));
   const vector<T>   bvector(jmax,b);
   const vector<vector<T>> bmatrix(imax,bvector);
   result = a - bmatrix;

   return result;
}

// double b - matrix a
//
template <typename T>
vector<vector<T>> operator-(const T& b, const vector<vector<T>>& a)
{
   const int imax = a.size();
   const int jmax = a[0].size();

   vector<vector<T>> result(imax,vector<T>(jmax));
   const vector<T>   bvector(jmax,b);
   const vector<vector<T>> bmatrix(imax,bvector);
   result = bmatrix - a;

   return result;
}

// - matrix a
//
template <typename T>
vector<vector<T>> operator-(const vector<vector<T>> &a)
{
   const int imax = a.size();
   const int jmax = a[0].size();
   
   vector<vector<T>> result(imax,vector<T>(jmax));
   result = -1.0*a;

   return result;
}

// matrix a * matrix b
//
template <typename T>
vector<vector<T>> operator*(const vector<vector<T>>& a, const vector<vector<T>>& b)
{
   const int imax = a.size();
   const int jmax = a[0].size();
   assert(a.size() == b.size());
   assert(a[0].size() == b[0].size());

   /*
   vector<vector<T>> result;
   vector<T> rvec;
   
   for (auto i=0; i<imax; i++) {
       rvec.reserve(jmax);
       transform(a.at(i).begin(), a.at(i).end(), b.at(i).begin(),
                 back_inserter(rvec), multiplies<T>());
       if(i==0) {
	   result.resize(1,rvec);
       }
       else { 
           result.push_back(rvec);
       }
       rvec.clear();
   }
   */

   vector<vector<T>> result(imax,vector<T>(jmax));
   
   for (auto i=0; i<imax; i++) {
      for (auto j=0; j<jmax; j++) {
         result[i][j] = a[i][j]*b[i][j];
      }
   }

   return result;
}

// double b * matrix a 
// 
template <typename T>
vector<vector<T>> operator*(const T b, const vector<vector<T>>& a)
{
   const int imax = a.size();
   const int jmax = a[0].size();
   
   vector<vector<T>> result(imax,vector<T>(jmax));
   const vector<T>   bvector(jmax,b);
   const vector<vector<T>> bmatrix(imax,bvector);
   result = a*bmatrix;

   return result;
}

// matrix a * double b 
// 
template <typename T>
vector<vector<T>> operator*(const vector<vector<T>>& a, const T b)
{
   const int imax = a.size();
   const int jmax = a[0].size();
   
   vector<vector<T>> result(imax,vector<T>(jmax));
   const vector<T>   bvector(jmax,b);
   const vector<vector<T>> bmatrix(imax,bvector);
   result = a*bmatrix;

   return result;
}

// matrix a / matrix b
//
template <typename T>
vector<vector<T>> operator/(const vector<vector<T>>& a, const vector<vector<T>>& b)
{
   const int imax = a.size();
   const int jmax = a[0].size();
   assert(a.size() == b.size());
   assert(a[0].size() == b[0].size());
   
   /*
   vector< vector<T> > result;
   vector<T> rvec;
   
   for (auto i=0; i<imax; i++) {
       rvec.reserve(jmax);
       transform(a.at(i).begin(), a.at(i).end(), b.at(i).begin(),
                 back_inserter(rvec), divides<T>());
       if(i==0) {
	   result.resize(1,rvec);
       }
       else { 
           result.push_back(rvec);
       }
       rvec.clear();
   }
   */

   vector<vector<T>> result(imax,vector<T>(jmax));
   
   for (auto i=0; i<imax; i++) {
      for (auto j=0; j<jmax; j++) {
         result[i][j] = a[i][j]/b[i][j];
      }
   }

   return result;
}

// matrix a / double b 
// 
template <typename T>
vector<vector<T>> operator/(const vector<vector<T>>& a, const T b)
{
   const int imax = a.size();
   const int jmax = a[0].size();
   
   vector<vector<T>> result(imax,vector<T>(jmax));
   const vector<T>   bvector(jmax,b);
   const vector<vector<T>> bmatrix(imax,bvector);
   result = a/bmatrix;

   return result;
}

// double b / matrix a 
// 
template <typename T>
vector<vector<T>> operator/(const T b, const vector<vector<T>>& a)
{
   const int imax = a.size();
   const int jmax = a[0].size();
   
   vector<vector<T>> result(imax,vector<T>(jmax));
   const vector<T>   bvector(jmax,b);
   vector<vector<T>> bmatrix(imax,bvector);
   result = bmatrix/a;

   return result;
}

//////////////////////////////////////////
//
//    vector operators
//


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
// exp(), sqrt(), tanh(), log(), cos(), sin() ...
//

vector<double> exp(const vector<double> &fin);

vector<double> sqrt(const vector<double> &fin);

vector<double> tanh(const vector<double> &fin);

vector<double> log(const vector<double> &fin);

vector<double> cos(const vector<double> &fin);

vector<double> sin(const vector<double> &fin);

vector<double> pow(const vector<double> &fin, const double exponent);

vector<double> abs(const vector<double> &fin);

double min(const vector<double> &fin);

double max(const vector<double> &fin);

double vanleer(const double a, const double b);
double minmod(const double a, const double b);
double superbee(const double a, const double b);


// operators for vector<vector>
//
vector<vector<double>> exp(const vector<vector<double>> &fin);

vector<vector<double>> sqrt(const vector<vector<double>> &fin);

vector<vector<double>> tanh(const vector<vector<double>> &fin);

vector<vector<double>> log(const vector<vector<double>> &fin);

vector<vector<double>> cos(const vector<vector<double>> &fin);

vector<vector<double>> sin(const vector<vector<double>> &fin);

vector<vector<double>> pow(const vector<vector<double>> &fin, const double exponent);

vector<vector<double>> abs(const vector<vector<double>> &fin);

double min(const vector<vector<double>> &fin);

double max(const vector<vector<double>> &fin);

#endif

