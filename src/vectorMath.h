/***
 *
 * template for doing vector math oporations
 *  such as +, -, *, and /
 *
***/

#ifndef vectorMath_h
#define vectorMath_h

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

// double b + vector a
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

#endif

