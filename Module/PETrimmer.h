////////////////////////////////////////////////////////////////////////
//
// FILL_IN: Comparisson package
//
// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef PETRIMMER_H
#define PETRIMMER_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include <vector>
#include "TTree.h"

#include "StandardNtuple/NtpStRecord.h"

#include "NueAna/NueAnalysisCuts.h"

using namespace std;


class TH1F;
class TFile;
class NtpStRecord;

class PETrimmer : public JobCModule
{
public:
  PETrimmer();
  ~PETrimmer();

public:
  // Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);
  void TrimRecord(NtpStRecord *str);

  void EndJob();
  void BeginJob();
  void Config(const Registry&);
  const Registry& DefaultConfig() const;

  void DumpIt(NtpStRecord * str);

int CleanEvent(NtpSREvent * evt,std::vector< NtpSRShower* > shw, std::vector< NtpSRTrack* > trk, std::vector<const NtpSRStrip* > stp, int *pepass, float *stripe_trl, float* stripe_shw, int * goodshower, int*goodtrack);

int CleanShower(NtpSRShower * evt, std::vector<const NtpSRStrip* > stp, int *pepass, float *stripe);
int CleanTrack(NtpSRTrack * evt, std::vector<const NtpSRStrip* > stp, int *pepass, float *stripe);


private:

 
  float pecut;
  float updateEventEnergy;
};
#endif // COMPAREALL_H
////////////////////////////////////////////////////////////////////////
