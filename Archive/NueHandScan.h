////////////////////////////////////////////////////////////////////////
// $Id: NueHandScan.h,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef NUEHANDSCAN_H
#define NUEHANDSCAN_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#include "CandNtupleSR/NtpSREvent.h"
#include <string>
#include <map>
#endif

class TFile;

class NueHandScan : public JobCModule
{
 public:
  NueHandScan();
  ~NueHandScan();

  JobCResult Ana(const MomNavigator* mom);
  void EndJob();
  Bool_t PassCuts(NtpSREvent *,Int_t);

  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);
  
private:

  std::map<Int_t,Int_t> runMap;
  std::map<Int_t,Int_t> subrunMap;
  std::map<Int_t,Int_t> snarlMap;
  std::map<Int_t,Int_t> eventMap;
  std::map<Int_t,Int_t> passMap;
  std::vector<std::string> ntupleFileNames;

  Int_t allEntries;
  Int_t randomSeed;
  Double_t fracPassed;

  Double_t preScaleFactorNue;
  Double_t preScaleFactorNuMu;
  Double_t preScaleFactorNuTau;
  Double_t preScaleFactorBNue;
  Double_t preScaleFactorNC;

  Float_t nNueCC;
  Float_t nNueNC;
  Float_t nNuMuCC;
  Float_t nNuMuNC;
  Float_t nNuTauCC;
  Float_t nNuTauNC;
  Float_t nBeamNueCC;
  Float_t nBeamNueNC;

  Float_t nPassNueCC;
  Float_t nPassNueNC;
  Float_t nPassNuMuCC;
  Float_t nPassNuMuNC;
  Float_t nPassNuTauCC;
  Float_t nPassNuTauNC;
  Float_t nPassBeamNueCC;
  Float_t nPassBeamNueNC;

  std::string fileTag;
  Bool_t firstPass;

};
#endif // NUEHANDSCAN_H
////////////////////////////////////////////////////////////////////////
