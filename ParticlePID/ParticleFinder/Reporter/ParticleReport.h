////////////////////////////////////////////////////////////////////////
/// $Id: ParticleReport.h,v 1.1 2009/06/19 14:32:39 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef PARTICLEREPORT_H
#define PARTICLEREPORT_H
#include "JobControl/JobCModule.h"

#include "NueAna/ParticlePID/ParticleFinder/Reporter/ParticleReportHelper.h"

class ParticleReport : public JobCModule
{
public:
  ParticleReport();
  ~ParticleReport();

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

   ParticleReportHelper * reporter;


  // Module configuration

  // Module member data

};
#endif // PARTICLEREPORT_H
////////////////////////////////////////////////////////////////////////
