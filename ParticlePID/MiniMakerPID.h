#ifndef MiniMakerPIDa //must end in lower case on osx?!
#define MiniMakerPIDa

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
                                                       
#include <iostream>
#include <string>
#include <vector>                                                                                                                                        
#include "NueAna/NueRecord.h"
#include "NueAna/NuePOT.h"
#include "NueAna/NueAnalysisCuts.h"
#include "Registry/Registry.h"
#include "NueAna/ParticlePID/NueMiniPID.h"

using std::string;


class MiniMakerPID
{
  public:
     MiniMakerPID();
     void AddFiles(string files);
     void SetOutputFile(string name);
     void Config(Registry &r);
     void SetCuts(string type, int level = 0);
     bool EvaluateCuts(NueRecord *nr, NueMiniPID *nm);
     
     void SetIsMRCC(int ismrcc=0){fIsMRCC=ismrcc;};

     void RunMiniMakerPID();

     void SetDeltaMSquare(double dm2);
     void SetUe3Square(double dUe32);
     void SetTheta23(double t23);
     void SetBaseline(double bl);
     void RecalculatePOT() {fOverWritePOT = true;}; 
     void SetReweight(int e=0){  fReweight = e;};

  private:
     string outFile;
     bool outSet;
     vector<string> files;
     NueAnalysisCuts fCuts;

     double fDeltaMSquare;
     double fUe3Square;
     double fTheta23;
     double fBaseline;
     bool fReweight;

     string cutSet;
     int cutLevel;

     bool fOverWritePOT;
     int fIsMRCC;
};

#endif


