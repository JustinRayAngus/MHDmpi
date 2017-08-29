/***
 *
 * HDF5 data file class
 *
***/

#ifndef HDF5dataFile_h
#define HDF5dataFile_h

#include "H5Cpp.h"

using namespace std;
#ifndef H5_NO_NAMESPACE
   using namespace H5;
#endif

class HDF5dataFile
{
public:
   void setOutputFile(const string&);
    
   void add(double &sclVar, const char *varName, int grow = 0);
   void add(vector<double> &vecVar, const char *varName, int grow = 0); 
   void add(vector<vector<double>> &vecVar, const char *varName, int grow = 0); 
   
   void writeAll(); // appends all growing variables in output file

private:
   string outputFile;
   
   bool varExists(const H5File&, const char*); // probably obsolete now
   bool varAdded(const string &name); 
   
   void writeScl(const double&, const char*, const bool&);   
   void writeVec(const vector<double>&, const char*, const bool&);  
   void writeVecVec(const vector<vector<double>>&, const char*, const bool&); 
   
   void appendVecInOutput(vector<double>*, const char*);
   void appendSclInOutput(double*, const char*);

   template<class T>
      struct VarStr {
         T *ptr;
         string name;
         bool grow;
      };
   vector< VarStr<double> > sclVar_arr;
   vector< VarStr< vector<double> > > vecVar_arr;
   vector< VarStr< vector<vector<double>> > > vecvecVar_arr;

};

void HDF5dataFile::setOutputFile(const string& thisOutputFile)
{
   outputFile = thisOutputFile;
   H5File file(outputFile.c_str(), H5F_ACC_TRUNC); // overwrites old
   file.close();
   cout << "\nOutput file " << outputFile.c_str() << " created " << endl;
}

void HDF5dataFile::add(double &varData, const char *varName, int grow) 
{   
   if(varAdded(string(varName))) {
      cout << "ERROR: Output variable " << string(varName) << " already added to datafile " << endl;
      exit (EXIT_FAILURE);
   }
   VarStr<double> scl;

   scl.ptr  = &varData;
   scl.name = string(varName);
   scl.grow = (grow>0) ? true:false;

   sclVar_arr.push_back(scl); // add variable to the list
   writeScl(varData, varName, grow); // add variable to the output file
}

void HDF5dataFile::add(vector<double> &varData, const char *varName, int grow) 
{
   if(varAdded(string(varName))) {
      cout << "ERROR: Output variable " << string(varName) << " already added to datafile " << endl;
      exit (EXIT_FAILURE);
   }
   VarStr< vector<double> > vec;

   vec.ptr  = &varData;
   vec.name = string(varName);
   vec.grow = (grow>0) ? true:false;

   vecVar_arr.push_back(vec); // add variable to the list
   writeVec(varData, varName, grow); // add variable to the output file
}

void HDF5dataFile::add(vector<vector<double>> &varData, const char *varName, int grow) 
{
   if(varAdded(string(varName))) {
      cout << "ERROR: Output variable " << string(varName) << " already added to datafile " << endl;
      exit (EXIT_FAILURE);
   }
   VarStr< vector<vector<double>> > vecvec;

   vecvec.ptr  = &varData;
   vecvec.name = string(varName);
   vecvec.grow = (grow>0) ? true:false;

   vecvecVar_arr.push_back(vecvec); // add variable to the list
   writeVecVec(varData, varName, grow); // add variable to the output file
}

void HDF5dataFile::writeAll()
{
   //cout << "Growing scalars are: " << endl;
   for(vector< VarStr<double> >::iterator it = sclVar_arr.begin(); it != sclVar_arr.end(); it++) {
      if(it->grow==1) {
         //cout << it->name << endl;
         appendSclInOutput(it->ptr, it->name.c_str());
      }
   }
   //cout << "Growing vectors are: " << endl;
   for(vector< VarStr< vector<double> > >::iterator it = vecVar_arr.begin(); it != vecVar_arr.end(); it++) {
      if(it->grow==1) {
         //cout << it->name << endl;
         appendVecInOutput(it->ptr, it->name.c_str());
      }
   }
}

/***
 *
 * HDF5 private functions
 *
***/

bool HDF5dataFile::varExists(const H5File& file2, const char* thisVar)
{
   try {
      Exception::dontPrint();
	  DataSet dataset0 = DataSet(file2.openDataSet(thisVar));
	  dataset0.close();
	  return 1; // dataset does not exist in output file yet
   }
   catch( ... ) {  
      return 0; // dataset already exists in output file
   }
}

bool HDF5dataFile::varAdded(const string &name)
{
   for(vector< VarStr<double> >::iterator it = sclVar_arr.begin(); it != sclVar_arr.end(); it++) {
      if(name == it->name)
         return true;
   }

   for(vector< VarStr< vector<double> > >::iterator it = vecVar_arr.begin(); it != vecVar_arr.end(); it++) {
      if(name == it->name)
         return true;
   }
   
   for(vector< VarStr< vector<vector<double>> > >::iterator it = vecvecVar_arr.begin(); it != vecvecVar_arr.end(); it++) {
      if(name == it->name)
         return true;
   }
   
   return false;
}

void HDF5dataFile::writeScl(const double& varData, const char* varName, const bool& growVar)
{     
   // Open file and check to see if output variable already exists
   //
   H5File file(outputFile.c_str(), H5F_ACC_RDWR);
   if(varExists(file,varName)) { // redundant, already checked in parent function
      cout << "ERROR: Output variable " << varName << " already exists " << endl;
	  exit (EXIT_FAILURE);
   }
   
   // Need to copy vectors to an array for writing purposes
   //
   int RANK = 2;
   PredType varType = PredType::NATIVE_DOUBLE;
   
   if(growVar) { // create extendable dataset
      
      hsize_t dimsf[RANK];
      dimsf[0] = 1; 
      dimsf[1] = 1;
      hsize_t mdimsf[2] = {1, H5S_UNLIMITED};
      DataSpace mdataspace(RANK, dimsf, mdimsf); // memory dataspace
      
      DSetCreatPropList cparms;
      hsize_t chunk_dims[2]={1,1}; // chunk size doesn't effect here
      cparms.setChunk(RANK,chunk_dims);
      
      DataType datatype(varType);
      DataSet dataset = file.createDataSet(varName, datatype, mdataspace, cparms);
      
      hsize_t size[2] = {1, 1};
      dataset.extend(size);
      DataSpace dataspace = dataset.getSpace();
      hsize_t offset[2] = {0,0};
      hsize_t dims[2] = {1, 1};
      dataspace.selectHyperslab(H5S_SELECT_SET,dims,offset); // data dataspace
               
      dataset.write(&varData, varType, mdataspace, dataspace);
      cout << "Extendable scalar " << varName << " added to " << outputFile << endl;    
           
      // close opened stuff
      //
      dataset.close();  
      mdataspace.close();    
      dataspace.close();
      file.close();   
   }
   else { // create non-extendable dataset
   
      // Write data to hdf5 file and try to catch errors
      //
      //H5S_class_t type = H5S_SCALAR;
      PredType varType = PredType::NATIVE_DOUBLE;
      //Exception::dontPrint();   
      DataSpace dataspace(H5S_SCALAR);
      DataType datatype(varType);
      DataSet dataset = file.createDataSet(varName, datatype, dataspace);
      dataset.write(&varData, varType); 
      cout << "Non-extendable scalar " << varName << " written to " << outputFile << endl; 
   
      // close opened stuff
      //
      dataset.close();      
      dataspace.close();
      file.close();    
   }
}

void HDF5dataFile::writeVec(const vector<double>& varData, const char* varName, const bool& growVar)
{     
   // Open file and check to see if output variable already exists
   //
   H5File file(outputFile.c_str(), H5F_ACC_RDWR);
   if(varExists(file,varName)) { // redundant, already checked in parent function
      cout << "ERROR: Output variable " << varName << " already exists " << endl;
	  exit (EXIT_FAILURE);
   }
   
   // Need to copy vectors to an array for writing purposes
   //
   int RANK;
   PredType varType = PredType::NATIVE_DOUBLE;
   const int Nvar = varData.size(); //cout << varData.size() << endl;
   double data[Nvar];
   for (int i=0; i<Nvar; i++) {
      data[i] = varData[i];
   } 
    
   if(growVar) { // create extendable dataset
      
      RANK = 2;
      hsize_t dimsf[RANK];
      dimsf[0] = Nvar; 
      dimsf[1] = 1;
      hsize_t mdimsf[2] = {H5S_UNLIMITED, varData.size()};
      DataSpace mdataspace(RANK, dimsf, mdimsf); // memory dataspace
      
      DSetCreatPropList cparms;
      hsize_t chunk_dims[2]={1,varData.size()}; // chunk size doesn't effect here
      cparms.setChunk(RANK,chunk_dims);
      
      DataType datatype(varType);
      DataSet dataset = file.createDataSet(varName, datatype, mdataspace, cparms);
      
      hsize_t size[2] = {1,varData.size()};
      dataset.extend(size);
      DataSpace dataspace = dataset.getSpace();
      hsize_t offset[2] = {0,0};
      hsize_t dims[2] = {1,varData.size()};
      dataspace.selectHyperslab(H5S_SELECT_SET,dims,offset); // data dataspace
               
      dataset.write(data, varType, mdataspace, dataspace);
      cout << "Extendable vector " << varName << " added to " << outputFile << endl;    
           
      // close opened stuff
      //
      dataset.close();
      mdataspace.close();      
      dataspace.close();
      file.close();   
   }
   else { // create non-extendable dataset
   
      RANK = 1;
      hsize_t dimsf[RANK];
      dimsf[0] = Nvar; //varData.size();
      DataSpace dataspace(RANK, dimsf);
      DataType datatype(varType);
      DataSet dataset = file.createDataSet(varName, datatype, dataspace);
      dataset.write(data, varType);
      cout << "Non-extendable vector " << varName << " added to " << outputFile << endl;   
           
      // close opened stuff
      //
      dataset.close();      
      dataspace.close();
      file.close();     
   }
}


void HDF5dataFile::writeVecVec(const vector<vector<double>>& varData, const char* varName, const bool& growVar)
{     
   // Open file and check to see if output variable already exists
   //
   H5File file(outputFile.c_str(), H5F_ACC_RDWR);
   if(varExists(file,varName)) { // redundant, already checked in parent function
      cout << "ERROR: Output variable " << varName << " already exists " << endl;
	  exit (EXIT_FAILURE);
   }
   
   // Need to copy vectors to an array for writing purposes
   //
   int RANK;
   PredType varType = PredType::NATIVE_DOUBLE;
   const int numVars = varData.size(); //cout << varData.size() << endl;
   const int Nvar = varData[0].size();
   double data[numVars][Nvar];
   for (auto j=0; j<numVars; j++) {
      vector<double> thisvarData = varData[j];
      for (int i=0; i<Nvar; i++) {
         data[j][i] = thisvarData[i];
      } 
   }
    
//    if(growVar) { // create extendable dataset
//       
//       RANK = 2;
//       hsize_t dimsf[RANK];
//       dimsf[0] = Nvar; 
//       dimsf[1] = numVars;
//       hsize_t mdimsf[RANK] = {H5S_UNLIMITED, varData.size()};
//       DataSpace mdataspace(RANK, dimsf, mdimsf); // memory dataspace
//       
//       DSetCreatPropList cparms;
//       hsize_t chunk_dims[2]={1,varData.size()}; // chunk size doesn't effect here
//       cparms.setChunk(RANK,chunk_dims);
//       
//       DataType datatype(varType);
//       DataSet dataset = file.createDataSet(varName, datatype, mdataspace, cparms);
//       
//       hsize_t size[RANK] = {1,varData.size()};
//       dataset.extend(size);
//       DataSpace dataspace = dataset.getSpace();
//       hsize_t offset[RANK] = {0,0};
//       hsize_t dims[RANK] = {1,varData.size()};
//       dataspace.selectHyperslab(H5S_SELECT_SET,dims,offset); // data dataspace
//                
//       dataset.write(data, varType, mdataspace, dataspace);
//       cout << "Extendable vector " << varName << " added to " << outputFile << endl;    
//            
//       // close opened stuff
//       //
//       dataset.close();
//       mdataspace.close();      
//       dataspace.close();
//       file.close();   
//    }
//    else { // create non-extendable dataset
   
      RANK = 2;
      hsize_t dimsf[RANK];
      dimsf[0] = numVars;
      dimsf[1] = Nvar;
      DataSpace dataspace(RANK, dimsf);
      DataType datatype(varType);
      DataSet dataset = file.createDataSet(varName, datatype, dataspace);
      dataset.write(data, varType);
      cout << "Non-extendable vector vector " << varName << " added to " << outputFile << endl;   
           
      // close opened stuff
      //
      dataset.close();      
      dataspace.close();
      file.close();     
//   }
}



void HDF5dataFile::appendSclInOutput(double* varData, const char* varName)
{  
   H5File file(outputFile.c_str(), H5F_ACC_RDWR);
   // Need to copy vectors to an array for writing purposes
   //
   int RANK = 2;
   PredType varType = PredType::NATIVE_DOUBLE;
   const int Nvar = 1; 
   
   // open extendable dataset and get ID's to stuff
   //
   DataSet dataset = file.openDataSet(varName);
   hid_t fileID, dsetID, dspaceID;
   fileID = H5Fopen(outputFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   dsetID = H5Dopen(fileID, varName, H5P_DEFAULT);
   dspaceID = H5Dget_space(dsetID);
   // get dimensions of extendable dataset
   //
   const int ndims = H5Sget_simple_extent_ndims(dspaceID);
   if(ndims == 1) {
      cout << "ERROR: Trying to extend non-extendable variable " 
		 << varName << endl;
	exit (EXIT_FAILURE);
   }
   hsize_t dims[ndims];
   H5Sget_simple_extent_dims(dspaceID,dims,NULL); // sets dims
   
   // extend dataset and write
   //
   hsize_t offset[RANK], size[RANK];
   offset[0] = 0;
   offset[1] = dims[1]; // equal to total number of previous time steps
   hsize_t dimsappend[RANK];
   dimsappend[0] = Nvar;
   dimsappend[1] = 1;
   dims[1] = dims[1] + dimsappend[1];
   size[0] = dims[0];
   size[1] = dims[1];
   dataset.extend(size);
   DataSpace fspace = dataset.getSpace();
   fspace.selectHyperslab(H5S_SELECT_SET, dimsappend, offset);
   DataSpace mspace(RANK,dimsappend);
   dataset.write(varData, varType, mspace, fspace);
   //cout << "Extendable scalar " << varName << " updated in " << outputFile << endl; 
	
   // close opened stuff
   //
   H5Fclose(fileID);  
   H5Fclose(dsetID);  
   dataset.close();      
   fspace.close();
   mspace.close();
   file.close();
}

void HDF5dataFile::appendVecInOutput(vector<double>* varData, const char* varName)
{  
   H5File file(outputFile.c_str(), H5F_ACC_RDWR);
   // Need to copy vectors to an array for writing purposes
   //
   int RANK = 2;
   PredType varType = PredType::NATIVE_DOUBLE;
   const int Nvar = varData->size(); //cout << varData.size() << endl;
   double data[Nvar];
   for (int i=0; i<Nvar; i++) {
      data[i] = (*varData)[i];
   } 
   
   // open extendable dataset and get ID's to stuff
   //
   DataSet dataset = file.openDataSet(varName);
   hid_t fileID, dsetID, dspaceID;
   fileID = H5Fopen(outputFile.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
   dsetID = H5Dopen(fileID, varName, H5P_DEFAULT);
   dspaceID = H5Dget_space(dsetID);
   
   // get dimensions of extendable dataset
   //
   const int ndims = H5Sget_simple_extent_ndims(dspaceID);
   if(ndims == 1) {
      cout << "ERROR: Trying to extend non-extendable variable " 
		 << varName << endl;
	exit (EXIT_FAILURE);
   }
   hsize_t dims[ndims];
   H5Sget_simple_extent_dims(dspaceID,dims,NULL); // sets dims
  
   // extend dataset and write
   //
   hsize_t offset[RANK], size[RANK];
   offset[0] = dims[0];
   offset[1] = 0; // equal to total number of previous time steps
   hsize_t dimsappend[RANK];
   dimsappend[0] = 1;
   dimsappend[1] = Nvar;
   dims[1] = dims[1]; //+ dimsappend[1];
   size[0] = dims[0] + dimsappend[0];
   size[1] = dims[1];
   dataset.extend(size);
   DataSpace fspace = dataset.getSpace();
   fspace.selectHyperslab(H5S_SELECT_SET, dimsappend, offset);
   DataSpace mspace(RANK,dimsappend);
   dataset.write(data, varType, mspace, fspace);
   //cout << "Extendable vector " << varName << " updated in " << outputFile << endl; 
	
   // close opened stuff
   //
   H5Fclose(fileID);  
   H5Fclose(dsetID);  
   dataset.close();      
   fspace.close();
   mspace.close();
   file.close();
}


#endif


