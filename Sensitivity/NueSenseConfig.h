#ifndef NUESENSECONFIG_H
#define NUESENSECONFIG_H

/******************************************************************
                                                                                
Input:
  5 numbers -> can't make use of delta term
                                                                                
  4 numbers and 1 energy histogram (oscillated) + osc params used
  4 numbers oscillated to standard CC result
  POT for these entries
                                                                                
  4 numbers (oscillated) and 1 energy histogram (unoscillated)
  POT for these numberse
                                                                                
  file with the selection (the file should include a pottree anyway)
  these files will be far mc files
  TChain to a trimmed nue file with just the selection
                                                                                
  range in delta
  range in sin(2th13)
  range in delta m23
                                                                                
  value for delta m12
  value for sin theta 23
  value for sin theta 12
  value for the density
                                                                                
  Things to know:
    systematic error
    inflate NC par
    target POT range
                                                                                
Ways to oscillate spectrum
 - simple -> no matter effect 2 neutrino approx
 - matter effect
                                                                                
Output -> chi2 surface
                                                                                
Input files:  .inp along the lines of skzp seems best
*****************************************************************/

#include <string>
#include <vector>

struct NSCDataParam{
   bool isfixed;
   float start;
   float end;
};

class NSCErrorParam{
 public:
  NSCErrorParam();
  float bg_systematic;
  float sig_systematic;
  float scale[5];
};

class NueSenseConfig
{
 public:
   enum SenseMethod{
      kAllNumbers = 1,
      kNueHist = 2,
      kAnaNueFiles =3,
      kNumbers = 4
   };

   NueSenseConfig();
   NueSenseConfig(std::string input);
   void ReadInput(std::string input);
   bool CheckConfig();
   void Reset();

   //Retrieving Settings
   double GetPOT() {return fPOT;}

   //Input Data
   int GetDataMethod()  {return fDataMethod;}

   bool ShouldDeOsc() {return fMustDeOsc; }
   float GetOldUe3Square() { return fOscUe3Square;}
   float GetOldDeltaMSquare() { return fOscDeltaMS23; }

   float GetNumber(int i) {return number[i];};
   std::string GetNueHistFile() { return signuehfile; }
   std::string GetNueHistName() { return signuehname; }
   int GetNumFiles() {return datafiles.size(); } ;
   std::string GetFile(int i) {return datafiles[i];};

   //Oscillation Parameters
   NSCDataParam GetDeltaMS23(){ return DeltaMS23;}
   NSCDataParam GetDelta()    { return Delta;}
   NSCDataParam GetSinS2Th13(){ return SinS2Th13;}
   
   float GetDeltaMS12() { return DeltaMS12;};  //units of 1e-5 eV^2
   float GetSinS2Th12() { return SinS2Th12;};  //Sin^2(2Theta12)
   float GetSinS2Th23() { return SinS2Th23;};  //Sin^2(2Theta23)
   float GetDensity()   { return density; };

   //Error settings
   int GetNumberConfig()    {return errorsets.size(); }
   NSCErrorParam GetErrorConfig(int i) { return errorsets[i]; };

   std::vector<std::string> GetDataFiles() { return datafiles;};

private:
   void SetDataInfo(int nuClass, std::string line);
   void SetParamInfo(NSCDataParam &par, std::string line);

   int fDataMethod;
   double fPOT;        // 1e20
   bool found[5];      // given class type has data
   std::string signuehfile;  //File of the energy histogram
   std::string signuehname;  //Name of the energy histogram
                                                                                
   float fOscUe3Square;  // parameter used when preparing numbers
   float fOscDeltaMS23;  // parameter used when preparing numbers
   bool fMustDeOsc;           // Do I need to deoscillate?

   NSCDataParam DeltaMS23;   //units of 1e-3 eV^2
   NSCDataParam Delta;       //units of pi
   NSCDataParam SinS2Th13;   //Sin^2(2Theta13)

   float DeltaMS12;   //units of 1e-5 eV^2
   float SinS2Th12;   //Sin^2(2Theta12)
   float SinS2Th23;   //Sin^2(2Theta23)
   float density;

   float number[5];   // Information (number for given class)

   std::vector<std::string> datafiles;  //list of input files for data
   std::vector<NSCErrorParam> errorsets;
};
 

#endif 
