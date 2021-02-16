////////////////////////////////////////////////////////////////////////
// $Id: FilterPID.h,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef FILTERPID_H
#define FILTERPID_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif

class TH1F;
class TFile;

class FilterPID : public JobCModule
{
public:
  FilterPID();
  ~FilterPID();

public:
  // Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);

private:
  // Module member data

};
#endif // FILTERPID_H
////////////////////////////////////////////////////////////////////////
