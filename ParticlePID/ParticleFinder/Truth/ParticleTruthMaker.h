////////////////////////////////////////////////////////////////////////
/// $Id: ParticleTruthMaker.h,v 1.1 2009/06/19 14:32:40 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef PARTICLETRUTHMAKER_H
#define PARTICLETRUTHMAKER_H
#include "JobControl/JobCModule.h"

#include "NueAna/ParticlePID/ParticleFinder/Truth/ParticleTruthHelper.h"

class ParticleTruthMaker : public JobCModule
{
public:
  ParticleTruthMaker();
  ~ParticleTruthMaker();

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

  ParticleTruthHelper * helper;






};
#endif // PARTICLETRUTHMAKER_H
////////////////////////////////////////////////////////////////////////
