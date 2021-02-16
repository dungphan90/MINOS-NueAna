#ifndef NUEEXPGENERATOR_H
#define NUEEXPGENERATOR_H

#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TH1D.h"
                                                       
#include <iostream>
#include <string>
#include <vector>                                                                                                                                        
#include "NueAna/NueRecord.h"
#include "NueAna/NuePOT.h"
#include "NueAna/NueAnalysisCuts.h"
#include "Registry/Registry.h"

/***********************************
Takes in a point with systematics and rolls a set of
 experiments corresponding to that point.  At the moment only
 works for single count rolls but could be reworked to do spectrum 
 without too much difficulty

Stores the 68 and 90 CL hits as well 
***********************************/


struct Position{
   TDirectory *fDirectory;
   string id;
   double fEventNo[5];
   TH1D* nevent;
   TH1D* deltaLog;
//   TH1D* sysevents[6];
//   TH1D* statevents[5];
   double sixtyeight;
   double ninety;
   double threesigma;
}; 

class NueExpGenerator
{
  public:
     NueExpGenerator();
     void Run(string errfile, string data, string output);
     void SetOutputFile(string name) { outFile = name;};
     Position* GenerateExperimentSet(double* num, double*err);
     double CalculateDeltaLog(double nexp, double nobs, double b = -1);
     void WriteToFile(Position *pos);

     void SetSeed(double in) {fSeed = in;};
     void SetNumberOfExp(int num) {fNumberOfExp = num;};
     void SetScale(double in) {fScale = in;};
     void ReScale(double *num);
     void SetMinTotal(double num) {fMinTotal = num;};

     void CleanPos(Position* pos);

     void SetOffset(int in) { fOffSet = in; };
     void SetReadNumber(int in) {fReadUntil = in; };

     bool DetermineCuts(Position* pos);
     void SetMethod(int met) {fMethod = met;};
     void SetFitMethod(int met) {fFitMethod = met;};
     void SetObservation(int obs) { fObservation = obs; };

     double CalculateChi2withBestFit(double b, double s, double k, double errBg, double errK, int n);

  private:
     string outFile;
     vector<string> files;
     double fSeed;

     TRandom3 fRand;
     int fNumberOfExp;
     double fErrors[6];
     double fNumbers[6];
     vector<Position*> posList;
     double fMinTotal;

     double fOscPar[6];

     int fOffSet;
     int fReadUntil;

     double fScale;
     int fMethod;
     int fFitMethod;
     double fObservation;
};

#endif
