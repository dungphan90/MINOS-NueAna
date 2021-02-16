////////////////////////////////////////////////////////////////////////
//
// FILL_IN: Comparisson package
//
// boehm@physics.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef TRIMMOD_H
#define TRIMMOD_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include <vector>
#include "TTree.h"
#include "NueAna/NueAnalysisCuts.h"


using namespace std;


class TH1F;
class TFile;
class NueRecord;


class TrimModule : public JobCModule
{
public:
  TrimModule();
  ~TrimModule();

public:
  // Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);
  void EndJob();
  void BeginJob();

    bool PassesCuts(NueRecord* nr);
    bool PassesBeamCuts(NueRecord* nr);

  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);

private:
  int counter;
  int kept;

  NueAnalysisCuts fCuts;

  std::string kOutputFile;

  Int_t kReWeight;
  Float_t kTheta23;
  Float_t kUe3Square;
  Float_t kDeltaMSquare;
};
#endif // COMPAREALL_H
////////////////////////////////////////////////////////////////////////
