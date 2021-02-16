////////////////////////////////////////////////////////////////////////
//
// FILL_IN: Comparisson package
//
// boehm@physics.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef COMPAREALL_H
#define COMPAREALL_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include <vector>
#include "TTree.h"

using namespace std;


class TH1F;
class TFile;
class HistMan;
class NueRecord;


class CompareAll : public JobCModule
{
public:
  CompareAll();
  ~CompareAll();

public:
  // Analysis and Reconstruction methods
  JobCResult Ana(const MomNavigator* mom);
  void EndJob();
  void BeginJob();

    bool PassesCuts(NueRecord* nr);
    bool PassesBeamCuts(NueRecord* nr);
    TString MakeIdString(NueRecord *nr);
    void FillFromList(NueRecord* nr, TString id, Float_t weight);
    bool NeedsSpecialAttention(TString name, NueRecord *nr, Float_t &value);

  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);
private:
  int counter;
    vector<string> varName;
    vector<Float_t> beg;
    vector<Float_t> end;
    vector<Int_t> nbins;
    vector<TString> gtype;
  Int_t kHiPlaneTrackCut;
  Int_t kLoPlaneEventCut;
  Int_t kHiTrackLikeCut;
  Int_t kDPlaneCut;
  Int_t kLoPhNStripCut;
  Int_t kLoPhNPlaneCut;
  Float_t kHiEnergyCut;
  Float_t kLoEnergyCut;
  Float_t kHiEnergyShowerCut;
  Float_t kLoEnergyShowerCut;
  Float_t kPhStripCut;
  Float_t kPhPlaneCut;
  Float_t  kLoCurrentCut;
  Float_t  kLoHorBeamWidth;
  Float_t  kHiHorBeamWidth;
  Float_t  kLoVertBeamWidth;
  Float_t  kHiVertBeamWidth;
  Float_t  kLoNuTarZ;
  Float_t  kHiNuTarZ;
  Int_t kOscillate;
  std::string kOutputFile;
  
  Float_t fPOT[4]; //{far_data, far_mc, near_data, near_mc}
  Float_t fOscParams[4];   //{L, dm2, theta23, UE32}
  TTree * fHistRecord;

  enum POTID {far_mc, far_data, near_mc, near_data};
  
};
#endif // COMPAREALL_H
////////////////////////////////////////////////////////////////////////
