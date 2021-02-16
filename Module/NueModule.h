////////////////////////////////////////////////////////////////////////
// $Id: NueModule.h,v 1.17 2009/06/30 23:53:42 pawloski Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef NUEMODULE_H
#define NUEMODULE_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif

#include <string>
#include "NueAna/NueAnalysisCuts.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "Conventions/ReleaseType.h"
#include "Conventions/BeamType.h"
#include "CalDetDST/UberRecord.h"

#include "MuonRemoval/NtpMRRecord.h"

class TProfile2D;
class BeamSummary;
class MuParentHelper;
class Zbeam;
class Zfluk;
class Kfluk;
class MCReweight;
class SKZPWeightCalculator;

class NueModule : public JobCModule
{
public:
  NueModule();
  ~NueModule();

public:
  // Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);
  JobCResult Analyze(MomNavigator* mom, NtpStRecord* str, NtpMRRecord* mr, UberRecord* ur, NtpStRecord* oldst);
  void BeginJob();
  void EndJob();
  // Module configuration
  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);

  void LoadBeamSummaryData(const char* bd);
  void SetMRCCRunning(bool val = true) {mrccRun = val;};  
  BeamType::BeamType_t DetermineBeamType(VldContext vc);
  BeamType::BeamType_t DetermineBeamType(string file);

private:
  bool PassesBlindingCuts(NueRecord *nr);

  std::string tmpltfile;
  std::string beamstring;

  NueAnalysisCuts fCut;

  Int_t kDPlaneCut;
  Int_t kLoPhNStripCut;
  Int_t kLoPhNPlaneCut;
  Float_t kPhStripCut;
  Float_t kPhPlaneCut;
  Float_t  kCPhPlaneCut;

//  Int_t kBeamType;
  BeamType::BeamType_t kBeamType;
  string kMuPiDir;
  bool kFixMuParents;
  bool kOpenedMuPiFile;

  int counter;
  int passcounter;
  int passcutcounter;
  int failcounter;
  int writecounter;
  bool foundmeu;

  float MSTminsig;
  float MSTmaxmetric;
  float MSTminfarsig;
  float MSTmaxmetriclowz;

  float SIGMAPMEU;
  
  TProfile2D *MSTetemplate; //!
  TProfile2D *MSTbtemplate; //!
  TProfile2D *MSTemtemplate; //!
  TProfile2D *MSTbmtemplate; //!

  Double_t threshCut;
  std::string sasFileName;

  float pot;
  float MEUPERGEV;

  bool foundRelease;
  int release;
  string skzpcfg;
  string title;
  string filename;

  bool mrccRun;

  MuParentHelper *mupar;
  Zbeam *zbeam;
  Zfluk *zfluk;
  Kfluk *kfluk;
  SKZPWeightCalculator *skzpCalc;

  MCReweight *mcr;
  Registry *xsecreweightreg;

  string kPidName;
  double kPIDHighCut;
  double kPIDLowCut;
  double kCCPIDCut;
  double kHighECut;
  double kLowECut;
  string kCCPidName;

};
#endif // NUEMODULE_H
////////////////////////////////////////////////////////////////////////
