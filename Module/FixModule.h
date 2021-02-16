////////////////////////////////////////////////////////////////////////
//
// FILL_IN: Comparisson package
//
// boehm@physics.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef FIXMOD_H
#define FIXMOD_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include <vector>
#include "TTree.h"
#include "NueAna/NueAnalysisCuts.h"
#include <string>

using namespace std;


class TH1F;
class TFile;
class NueRecord;


class FixModule : public JobCModule
{
public:
  FixModule();
  ~FixModule();

public:
  // Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);
  void EndJob();
  void BeginJob();

  const Registry& DefaultConfig() const;
  void Config(const Registry& r);

private:
  int counter;
  int kept;
  std::string fInFile;
  std::string fOutFile;

};
#endif // COMPAREALL_H
////////////////////////////////////////////////////////////////////////
