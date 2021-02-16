////////////////////////////////////////////////////////////////////////
// $Id: CompareMST.h,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef COMPAREMST_H
#define COMPAREMST_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif

class TH1F;
class TFile;

class CompareMST : public JobCModule
{
public:
  CompareMST();
  ~CompareMST();

public:
  // Analysis and Reconstruction methods
  JobCResult Ana(const MomNavigator* mom);
  void EndJob();
  void BeginJob();

private:
  int counter;
};
#endif // COMPAREMST_H
////////////////////////////////////////////////////////////////////////
