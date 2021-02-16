#ifndef NUEGENCONFIG_H
#define NUEGENCONFIG_H

/**********************************
  Reads in the ranges for the possible values to "scan"

  range in delta m23
  range in delta cp
  range in sin theta 13

  SetList doesn't actually work yet, but the system is in place for it
                                                                                
Input files:  .inp along the lines of skzp seems best
*****************************************************************/

#include <string>
#include <vector>

struct GCDataParam{
   bool isfixed;
   double start;
   double end;
   double step;
};

struct GCCompleteSet{
    GCDataParam EventClass[6];
};

struct place
{
  double num[6];
};

class NueGenConfig
{
 public:
   typedef enum GenMethod{
      kFileList = 1,
      kSetList  = 2
   } GCTypes_t;

   NueGenConfig();
   NueGenConfig(std::string input);
   void ReadInput(std::string input);
   bool CheckConfig();
   void Reset();

   //Accessing the data
   bool LoadNextNumberSet(double* num);
   double * GetErrors() { return fErrors;};
   double * GetOscPar() { return fOscPar;};

   //Input Data
   int GetDataMethod()  {return fDataMethod;}

   void SetDataRange(int nuClass, std::string line, GCCompleteSet &set);
   void AddInputFile(std::string file)  {files.push_back(file);}

   void SetOffset(int in) { fOffSet = in; };
   void SetReadNumber(int in) {fReadUntil = in; };
 

private:
   void SetErrorVal(int pos, std::string line);

   int fOffSet;
   int fReadUntil;
   int fDataMethod;
   bool found[6];      // given class type has data
   std::vector<std::string> files;   // File of the energy histogram
                                                                                
   double fErrors[6]; 
   std::vector<GCCompleteSet> EventGroups;

   double fOscPar[6];
};
 

#endif 
