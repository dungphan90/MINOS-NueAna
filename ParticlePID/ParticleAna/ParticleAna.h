////////////////////////////////////////////////////////////////////////
/// $Id: ParticleAna.h,v 1.1 2009/06/19 14:32:37 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef PARTICLEANA_H
#define PARTICLEANA_H
#include "JobControl/JobCModule.h"


//#include "MCReweight/SKZPWeightCalculator.h"
#include "MCNtuple/NtpMCTruth.h"
//#include "MCReweight/Zbeam.h"
//#include "MCReweight/Zfluk.h"
//#include "MCReweight/Kfluk.h"

//#include "NueAna/ParticlePID/ParticleFinder/POT.h"

#include "NueAna/ParticlePID/ParticleFinder/ParticleBeamMon.h"

class ParticleAna : public JobCModule
{
public:
  ParticleAna();
  ~ParticleAna();

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

	void FillPot(ParticleBeamMon *bmon);

private:
  // Module configuration

  // Module member data
  
 // static SKZPWeightCalculator *skzpCalc;


    std::string kPOTTreeName;
   // POT *pot; 

	int lastrun;

	int DoMRCC;
	std::string outObjectName;

};
#endif // PARTICLEANA_H
////////////////////////////////////////////////////////////////////////

