#ifndef NUEEXPBUILDER
#define NUEEXPBUILDER

/****************************************************************
This class generates experiment samples using an input set
 of MC,  Takes a target number of poissson fluctuates it
****************************************************************/

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
                                                       
#include <iostream>
#include <string>
#include <vector>                                                                                                                                        
#include "NueAna/NueRecord.h"
#include "NueAna/NuePOT.h"
#include "OscProb/OscCalc.h"

#include "Registry/Registry.h"

class NueExpBuilder
{
  public:
     NueExpBuilder();

     void AddFiles(std::string files);

     void SetOutputFile(std::string name);

     void SetCuts(std::string type, int level = 0);
     bool EvaluateCuts(NueRecord *nr);

     //Oscillation Probabilities 
     void SetDeltaMSquare(double dm2);
     void SetUe3Square(double dUe32);
     void SetTheta23(double t23);
     void SetBaseline(double bl);

     void SetTargetPOT(double num) {fTargetPOT = num;};

     void SetMeanNumberOfEvents(double num);
     void SetNumberOfEvents(int num);

     void GenerateExperiment(int nExp);

     void SetSeed(int in) {fSeed = in;};
     void SetMaxProb(double max);
     double GetMaxProb();


  private:
     std::string outFile;
     bool outSet;
     std::vector<std::string> files;

     double fDeltaMSquare;
     double fUe3Square;
     double fTheta23;
     double fBaseline;
     bool fReweight;

     std::string cutSet;
     int cutLevel;

     OscCalc fOsc;

     int fNumberOfEvents;
     double fMeanNumber;
     double fMaxProb;
     double fTargetPOT;
     int fSeed;
};

#endif
