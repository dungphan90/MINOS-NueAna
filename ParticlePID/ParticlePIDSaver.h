////////////////////////////////////////////////////////////////////////
/// $Id: ParticlePIDSaver.h,v 1.3 2009/08/13 11:27:38 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef PARTICLEPIDSAVER_H
#define PARTICLEPIDSAVER_H
#include "JobControl/JobCModule.h"

#include "NueAna/ParticlePID/ParticleAna/PRecord.h"
#include "NueAna/NueRecord.h"

#include <string>

using std::string;

class ParticlePIDSaver : public JobCModule
{
public:
  ParticlePIDSaver();
  ~ParticlePIDSaver();

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
  
	int isMatched(NueRecord *nr, PRecord *pr);


	string kPOTTreeName;
	string kPOTTreeNameIn;




};
#endif // PARTICLEANA_H
////////////////////////////////////////////////////////////////////////

