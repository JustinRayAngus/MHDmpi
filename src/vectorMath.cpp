/***
 *
 * template for doing vector math oporations
 *  such as +, -, *, and /
 *
***/


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


#include "vectorMath.h"

using namespace std;


////////////////////////////////////////////////////////////////
//
// functions for other common math operators
// exp(), tanh, log(), cos(), sin() ...
//

vector<double> exp(const vector<double> &fin) {

   vector<double> result;
   const int imax = fin.size();
   result.resize(imax);

   //cout << "fin SIZE IS " << fin.size() << endl;
   //cout << "RESULT SIZE IS " << result.size() << endl;

   for (auto i=0; i<imax; i++) {
      result.at(i) = exp(fin.at(i));
   }

   return result;
}

vector<double> tanh(const vector<double> &fin) {

   vector<double> result;
   const int imax = fin.size();
   result.resize(imax);

   for (auto i=0; i<imax; i++) {
      result.at(i) = tanh(fin.at(i));
   }

   return result;
}

vector<double> log(const vector<double> &fin) {

   vector<double> result;
   const int imax = fin.size();
   result.resize(imax);

   for (auto i=0; i<imax; i++) {
      result.at(i) = log(fin.at(i));
   }

   return result;
}

vector<double> cos(const vector<double> &fin) {

   vector<double> result;
   const int imax = fin.size();
   result.resize(imax);

   for (auto i=0; i<imax; i++) {
      result.at(i) = cos(fin.at(i));
   }

   return result;
}

vector<double> sin(const vector<double> &fin) {

   vector<double> result;
   const int imax = fin.size();
   result.resize(imax);

   for (auto i=0; i<imax; i++) {
      result.at(i) = sin(fin.at(i));
   }

   return result;
}

vector<double> abs(const vector<double> &fin) {

   vector<double> result;
   const int imax = fin.size();
   result.resize(imax);

   for (auto i=0; i<imax; i++) {
      result.at(i) = abs(fin.at(i));
   }

   return result;
}

double min(const vector<double> &fin) {

   double result, localresult;
   const int imax = fin.size();

   localresult = fin.at(0);
   for (auto i=1; i<imax; i++) {
      if (fin.at(i)<=localresult) {
         localresult = fin.at(i);
      }
   }

   MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

   return result;
}

double max(const vector<double> &fin) {

   double result, localresult;
   const int imax = fin.size();

   localresult = fin.at(0);
   for (auto i=1; i<imax; i++) {
      if (fin.at(i)>=localresult) {
         localresult = fin.at(i);
      }
   }

   MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

   return result;
}


