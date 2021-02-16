////////////////////////////////////////////////////////////////////////
// $Id: MSTTemplate.h,v 1.1 2005/04/25 19:57:26 minoscvs Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#ifndef MSTTEMPLATE_H
#define MSTTEMPLATE_H
#ifndef JOBCMODULE_H
#include "JobControl/JobCModule.h"
#endif
#include<string>

class TProfile2D;
class TFile;

class MSTTemplate : public JobCModule
{
public:
  MSTTemplate();
  ~MSTTemplate();

public:
  // Analysis and Reconstruction methods
  JobCResult Ana(const MomNavigator* mom);

  void EndJob();
  void BeginJob();


  // Module configuration
  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);

private:
  int counter;
  std::string fname;
  TFile *fout;
  TProfile2D *nlambdanele;
  TProfile2D *nlambdanother;
  TProfile2D *mipdistele;
  TProfile2D *mipdistother;
  double *xbins;
  double *ybins;
  double *ymbins;
  double *zbins;

};
#endif // MSTTEMPLATE_H
////////////////////////////////////////////////////////////////////////
