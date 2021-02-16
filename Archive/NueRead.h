////////////////////////////////////////////////////////////////////////
// $Id: NueRead.h,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef NUEREAD_H
#define NUEREAD_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif

class TH1F;
class TFile;

class NueRead : public JobCModule
{
public:
  NueRead();
  ~NueRead();

public:
  // Analysis and Reconstruction methods
  JobCResult Ana(const MomNavigator* mom);
  void EndJob();

private:
  TFile *f;
  TH1F *hpar0;
  int counter;
  int passcounter;
  // Module member data

};
#endif // NUEREAD_H
////////////////////////////////////////////////////////////////////////
