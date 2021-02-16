////////////////////////////////////////////////////////////////////////
//
// FILL_IN: Comparisson package
//
// boehm@physics.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef SPILLTYPETRIMMER_H
#define SPILLTYPETRIMMER_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include <vector>
#include "TTree.h"
#include "NueAna/NueAnalysisCuts.h"

using namespace std;


class TH1F;
class TFile;
class NtpStRecord;

class SpillTypeFilter : public JobCModule
{
public:
  SpillTypeFilter();
  ~SpillTypeFilter();

public:
  // Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);
  void EndJob();
  void BeginJob();

  bool PassesCuts();
//  bool PassesBeamCuts(NueRecord* nr);

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
