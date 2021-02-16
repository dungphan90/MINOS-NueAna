////////////////////////////////////////////////////////////////////////
//
// FILL_IN: Comparisson package
//
// boehm@physics.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef COMPAREMD_H
#define COMPAREMD_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include <vector>

using namespace std;


class TH1F;
class TFile;
class HistMan;
class NueRecord;


class CompareMD : public JobCModule
{
public:
  CompareMD();
  ~CompareMD();

public:
  // Analysis and Reconstruction methods
  JobCResult Ana(const MomNavigator* mom);
  void EndJob();
  void BeginJob();

  bool PassesBeamCuts(NueRecord* nr);
  bool PassesCuts(NueRecord* nr);

  TString MakeIdString(NueRecord *nr);
  void FillFromList(NueRecord* nr, TString id, HistMan* hm, Float_t weight);
  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);

private:
  int counter;
  vector<TString> varName;
  vector<Float_t> beg;
  vector<Float_t> end;
  vector<Int_t> nbins;
  vector<TString> gtype;

  Int_t kNMCFiles;

};
#endif // COMPAREMD_H
////////////////////////////////////////////////////////////////////////
