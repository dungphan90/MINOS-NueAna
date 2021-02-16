////////////////////////////////////////////////////////////////////////
/// $Id: AnaTrim.h,v 1.1 2009/06/19 14:32:36 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef ANATRIM_H
#define ANATRIM_H
#include "JobControl/JobCModule.h"




class AnaTrim : public JobCModule
{
public:
  AnaTrim();
  ~AnaTrim();

public:
  // Handle job status changes
  void BeginJob();
  void EndJob();

  // Analysis and Reconstruction methods
  JobCResult Reco(MomNavigator* mom);

  // Module configuration
  const Registry& DefaultConfig() const;
  void            Config(const Registry& r);

  // User interface methods
  void Reset();



private:
  // Module configuration

  // Module member data


};
#endif // ANATRIM_H
////////////////////////////////////////////////////////////////////////

