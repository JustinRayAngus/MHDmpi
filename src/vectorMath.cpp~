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
// exp(), sqrt(), tanh, log(), cos(), sin() ...
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

vector<double> sqrt(const vector<double> &fin) {

   vector<double> result;
   const int imax = fin.size();
   result.resize(imax);

   //cout << "fin SIZE IS " << fin.size() << endl;
   //cout << "RESULT SIZE IS " << result.size() << endl;

   for (auto i=0; i<imax; i++) {
      result.at(i) = sqrt(fin.at(i));
      //if(result.at(i) != result.at(i)) {
      //    cout << "ERROR: sqrt() on vector returning complex results !!! " << endl;
      //    exit (EXIT_FAILURE);
      //}    
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

vector<double> pow(const vector<double> &fin, const double exponent) {

   vector<double> result;
   const int imax = fin.size();
   result.resize(imax);

   for (auto i=0; i<imax; i++) {
      result.at(i) = pow(fin.at(i), exponent);
      if(result.at(i) != result.at(i)) {
          cout << "ERROR: pow() on vector returning complex results !!! " << endl;
          cout << "ERROR: exponent = " << exponent << endl;
          cout << "ERROR: fin.at(i) = " << fin.at(i) << endl;
          cout << "ERROR: fin.at(i+1) = " << fin.at(i+1) << endl;
          cout << "ERROR: fin.at(i+2) = " << fin.at(i+2) << endl;
          cout << "ERROR: i = " << i << endl;
          exit (EXIT_FAILURE);
      }    
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

////////////////////////////////////////////////////
//
//       operators for vector<vector>
//
////////////////////////////////////////////////////

vector<vector<double>> exp(const vector<vector<double>> &fin) {

   const int imax = fin.size();
   const int jmax = fin[0].size();
   vector<vector<double>> result;
   result.assign(imax,vector<double>(jmax));

   for (auto i=0; i<imax; i++) {
       for (auto j=0; j<jmax; j++) {
           result[i][j] = exp(fin[i][j]);
       }
   }

   return result;
}

vector<vector<double>> sqrt(const vector<vector<double>> &fin) {

   const int imax = fin.size();
   const int jmax = fin[0].size();
   vector<vector<double>> result(imax,vector<double>(jmax));

   for (auto i=0; i<imax; i++) {
       for (auto j=0; j<jmax; j++) {
           result[i][j] = sqrt(fin[i][j]);
	   /*
	   if(result[i][j] != result[i][j]) {
	       cout << "ERROR: sqrt() on matrix returning complex results !!! " << endl;
               exit (EXIT_FAILURE);
	   }
           */	   
       }
   }

   return result;
}

vector<vector<double>> log(const vector<vector<double>> &fin) {

   const int imax = fin.size();
   const int jmax = fin[0].size();
   vector<vector<double>> result(imax,vector<double>(jmax));

   for (auto i=0; i<imax; i++) {
       for (auto j=0; j<jmax; j++) {
           result[i][j] = log(fin[i][j]);
       }
   }

   return result;
}

vector<vector<double>> tanh(const vector<vector<double>> &fin) {

   const int imax = fin.size();
   const int jmax = fin[0].size();
   vector<vector<double>> result(imax,vector<double>(jmax));

   for (auto i=0; i<imax; i++) {
       for (auto j=0; j<jmax; j++) {
           result[i][j] = tanh(fin[i][j]);
       }
   }

   return result;
}

vector<vector<double>> cos(const vector<vector<double>> &fin) {

   const int imax = fin.size();
   const int jmax = fin[0].size();
   vector<vector<double>> result(imax,vector<double>(jmax));

   for (auto i=0; i<imax; i++) {
       for (auto j=0; j<jmax; j++) {
           result[i][j] = cos(fin[i][j]);
       }
   }

   return result;
}

vector<vector<double>> sin(const vector<vector<double>> &fin) {

   const int imax = fin.size();
   const int jmax = fin[0].size();
   vector<vector<double>> result(imax,vector<double>(jmax));

   for (auto i=0; i<imax; i++) {
       for (auto j=0; j<jmax; j++) {
           result[i][j] = sin(fin[i][j]);
       }
   }

   return result;
}

vector<vector<double>> pow(const vector<vector<double>> &fin, const double exponent) {

   const int imax = fin.size();
   const int jmax = fin[0].size();
   vector<vector<double>> result(imax,vector<double>(jmax));

   for (auto i=0; i<imax; i++) {
       for (auto j=0; j<jmax; j++) {
           result[i][j] = pow(fin[i][j],exponent);
	   // check for imaginary stuff
	   /*
	   if(result[i][j] != result[i][j]) {
	       cout << "ERROR: pow() on matrix returning complex results !!! " << endl;
               exit (EXIT_FAILURE);
	   } 
	   */   
       }
   }

   return result;
}

vector<vector<double>> abs(const vector<vector<double>> &fin) {

   const int imax = fin.size();
   const int jmax = fin[0].size();
   vector<vector<double>> result(imax,vector<double>(jmax));

   for (auto i=0; i<imax; i++) {
       for (auto j=0; j<jmax; j++) {
           result[i][j] = abs(fin[i][j]);
       }
   }

   return result;
}

double min(const vector<vector<double>> &fin) {

   double result, thisfin, localresult;
   const int imax = fin.size();
   const int jmax = fin[0].size();

   localresult = fin[0][0];
   for (auto i=0; i<imax; i++) {
       for (auto j=0; j<jmax; j++) {
	   thisfin = fin[i][j];
           if (thisfin<=localresult) {
               localresult = thisfin;
           }
       }
   }

   MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

   return result;
}

double max(const vector<vector<double>> &fin) {

   double result, thisfin, localresult;
   const int imax = fin.size();
   const int jmax = fin[0].size();

   localresult = fin[0][0];
   for (auto i=0; i<imax; i++) {
       for (auto j=0; j<jmax; j++) {
	   thisfin = fin[i][j];
           if (thisfin>=localresult) {
               localresult = thisfin;
           }
       }
   }

   MPI_Allreduce(&localresult, &result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

   return result;
}


double vanleer(const double a, const double b) {

   // vanleer flux limiter = 2*a*b/(a+b) (if a and b same sign)
   // zero otherwise
   //
   double result;
   
   result = 2.0*a*b/(a+b);

   if(a<=0 && b>=0) result = 0.0;
   if(a>=0 && b<=0) result = 0.0;

   return result;
}

double minmod(const double a, const double b) {

   // minmod flux limiter = 1/2(sign(a)+sign(b))*min(|a|,|b|)
   //
   double result, signa, signb;
   
   signa = 0.0;
   signb = 0.0;
   if(a>0) signa = 1.0;
   if(a<0) signa = -1.0;
   if(b>0) signb = 1.0;
   if(b<0) signb = -1.0;

   result = 0.5*(signa + signb)*min(abs(a),abs(b));

   return result;
}

double superbee(const double a, const double b) {

   // superbee flux limiter = minmod(a,2b) if |a|>=|b|
   //                       = minmod(2a,b) otherwise
   //
   double result;
   
   if(abs(a)>=abs(b)) result = minmod(a,2*b);
   else result = minmod(2*a,b);

   return result;
}





