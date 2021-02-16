////////////////////////////////////////////////////////////////////////
/// $Id: ParticleFinder.h,v 1.2 2009/06/23 18:08:31 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////
#ifndef PARTICLEFINDER_H
#define PARTICLEFINDER_H
#include "JobControl/JobCModule.h"

#include "NueAna/ParticlePID/ParticleFinder/Finder.h"


#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "MCReweight/SKZPWeightCalculator.h"
#include "MCNtuple/NtpMCTruth.h"
#include "MCReweight/Zbeam.h"
#include "MCReweight/Zfluk.h"
#include "MCReweight/Kfluk.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleBeamMon.h"
#include "TF1.h"
#include "NueAna/ParticlePID/ParticleFinder/POT.h"


#include "NueAna/ParticlePID/ParticleFinder/ParticleBeamMonAna.h"

	const int nfront=5;  //number of hits to use for near det fully instrumentented front containment


class ParticleFinder : public JobCModule
{
public:
  ParticleFinder();
  ~ParticleFinder();

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

double OscillationProb(TF1* f2, int ntype, double NuE, double sinth23, double sin2th13);

  Finder finder;

	void CheckContainment(ParticleObjectHolder *p, double  front_z[2][nfront], double front_t[2][nfront]);


		 static TF1 * osceq;
  // Module configuration

  // Module member data
  
      std::string kPOTTreeName;
  
  static SKZPWeightCalculator *skzpCalc;

    //static POT *pot; 

	int lastrun;
	
	int IsInsideNearFiducial_Nue_Standard(double x, double y, double z, bool /*isMC*/);
	int IsInsideFarFiducial_Nue_Standard(double x, double y, double z, bool isMC);

		ParticleBeamMonAna bma;
	int recoCnt;

	std::string outObjectName;
	int DoMRCC;

	double pecutlevel;
	double MIPperPE;
	double MEUperGeV;
};
#endif // PARTICLEFINDER_H
////////////////////////////////////////////////////////////////////////

