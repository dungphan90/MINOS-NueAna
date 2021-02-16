////////////////////////////////////////////////////////////////////////
// $Id: NueSensitivity.h,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef NUESENSITIVITY_H
#define NUESENSITIVITY_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif

class TH1F;
class TFile;
class TF1;
class TH2F;

class NueSensitivity : public JobCModule
{
public:
  NueSensitivity();
  ~NueSensitivity();

public:
  // Analysis and Reconstruction methods
  void BeginJob();
  JobCResult Ana(const MomNavigator* mom);
  void Analysis();
  void EndJob();
  
  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);
  
private:
  // Module member data
  
  void SetPOT(); //translate #files -> #POT

  TF1 *nueAppear;
  TF1 *numuSurvive;
  
  TFile *systematicFile;
  TH2F *systematicHist_allSys;
  TH2F *systematicHist_oscSys;
  Double_t systematicHistNorm;

  Double_t MDCNearToFar;

  //for normalisation:
  Double_t MDCChallengePOT;
  
  Int_t nNuMuFiles;
  Int_t nNueFiles;
  Int_t nNuTauFiles;
  Int_t nNearFiles;
  Int_t nChallengeNearFiles;
  
  Double_t NuMuFilesPOT;
  Double_t NueFilesPOT;
  Double_t NuTauFilesPOT;
  Double_t NearFilesPOT;
  Double_t ChallengeNearFilesPOT;
  
  Int_t nNearUnknownEvents;
  Int_t nNearNueEvents;
  Int_t nNearNuMuEvents;
  Int_t nNearNuTauEvents;
  Int_t nNearBeamNueEvents;
  Int_t nNearNCEvents;

  Int_t nFarUnknownEvents;
  Int_t nFarNueEvents;
  Int_t nFarNuMuEvents;
  Int_t nFarNuTauEvents;
  Int_t nFarBeamNueEvents;
  Int_t nFarNCEvents;

  Int_t nThetaPoints;
  Double_t *theta13;

  Int_t nDeltaPoints;
  Double_t *delta23;

  Int_t currentRun;

};
#endif // NUESENSITIVITY_H
////////////////////////////////////////////////////////////////////////
