////////////////////////////////////////////////////////////////////////
// $Id: NueReweight.h,v 1.2 2005/05/21 16:26:09 vahle Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle
////////////////////////////////////////////////////////////////////////
#ifndef NUEREWEIGHT_H
#define NUEREWEIGHT_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include "MCReweight/MCReweight.h"
#include "NueAna/Reweight/NueRW.h"
#include "Registry/Registry.h"
class NtpMCTruth;


class NueReweight : public JobCModule
{
public:
  NueReweight();
  ~NueReweight();

public:
  // Handle job status changes
  void EndJob();

  int ReadRandom();
  virtual void Reset();

  // Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);

  // Module configuration
  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);


  NueRW *rw;
private:
  // Module member data
  MCReweight mcr;
  int counter;
  int rand;
  bool firstevent;
  Registry rwtconfig;

  ClassDef(NueReweight,0)
};
#endif // NUEREWEIGHT_H
////////////////////////////////////////////////////////////////////////
