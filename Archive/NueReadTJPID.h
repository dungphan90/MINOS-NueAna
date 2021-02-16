////////////////////////////////////////////////////////////////////////
// $Id: NueReadTJPID.h,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef NUEREADTJPID_H
#define NUEREADTJPID_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif

class TH1F;
class TFile;

class NueReadTJPID : public JobCModule
{
public:
  NueReadTJPID();
  ~NueReadTJPID();

public:
  // Analysis and Reconstruction methods
  JobCResult Ana(const MomNavigator* mom);

private:
  // Module member data

};
#endif // NUEREADTJPID_H
////////////////////////////////////////////////////////////////////////
