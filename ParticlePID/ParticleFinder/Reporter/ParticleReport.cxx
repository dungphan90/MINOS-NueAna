////////////////////////////////////////////////////////////////////////
///
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////
#include "NueAna/ParticlePID/ParticleFinder/Reporter/ParticleReport.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro

#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "CandNtupleSR/NtpSRPulseHeight.h"

#include "StandardNtuple/NtpStRecord.h"


#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObject.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle.h"



JOBMODULE(ParticleReport, "ParticleReport",
          "finds particles");

//......................................................................

ParticleReport::ParticleReport():reporter(0)  
{


reporter = new ParticleReportHelper();


  ///
  /// (Document me!)
  ///
}
//......................................................................

ParticleReport::~ParticleReport()
{
  ///
  /// (Document me!)
  ///

if(reporter) delete reporter;

}

//......................................................................

void ParticleReport::BeginJob()
{
  ///
  /// (Document me!)
  ///
}

//......................................................................

void ParticleReport::EndJob()
{
  ///
  /// (Document me!)
  ///
}

//......................................................................

JobCResult ParticleReport::Reco(MomNavigator* mom)
{
  ///
  /// (Document me!)
  ///

//MSG("ParticleReport",Msg::kDebug) << "reco called\n";

  NtpStRecord *str=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord"));


  //std::vector<ParticleObjectHolder *>* hf;
  std:vector<TObject* > hft =( mom->GetFragmentList("ParticleObjectHolder"));



  int make_new_holder=1;


  if(hft.size()>0)
  {
      make_new_holder=0;
  }
 

printf("size of found array %d\n",hft.size()); 

  std::vector<const NtpSREvent* > evts = str->GetEvents();
  std::vector<const NtpSRStrip* > stp = str->GetStrips();

  printf("Events %d \n",evts.size());
  if(evts.size()<1)
	{

		if(make_new_holder)
		{
                	ParticleObjectHolder * h = new ParticleObjectHolder(); 
                	mom->AdoptFragment(h);
		}

		return JobCResult::kPassed;  
	}


	printf("!!! %d events, %d  currently there \n", evts.size(), hft.size());

  //for(int ievt=0;ievt<evts.size();ievt++)
 // {
int ievt=0;

	reporter->Reset();


        ParticleObjectHolder *n=0;


        if(make_new_holder || hft.size() < ievt+1)
        {
                ParticleObjectHolder * h = new ParticleObjectHolder();
                mom->AdoptFragment(h);
                n=h;
        }else{
                n=dynamic_cast<ParticleObjectHolder *>(hft[ievt]);
        }


printf("partices in current evt %d array %d \n", ievt, n->particles->GetEntries());

	if(n->particlematch->GetEntries()>0)return JobCResult::kPassed;//continue; //dont rerun it


	//loop over found particles
	

		//loop over planes, strips in found particles

		for(int i=0;i<n->particletruth[0].GetEntries();i++)
			reporter->addtruth((ParticleTruthObject*)n->particletruth[0][i]);

	        for(int i=0;i<n->particles[0].GetEntries();i++)
			reporter->addfound((ParticleObject *)n->particles[0][i]);

		//record 
		reporter->Process(*n);


 // }



  return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................

const Registry& ParticleReport::DefaultConfig() const
{
  ///
  /// Supply the default configuration for the module
  ///
  static Registry r; // Default configuration for module

  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  // Set values in configuration
  r.UnLockValues();
  r.LockValues();

  return r;
}

//......................................................................

void ParticleReport::Config(const Registry& r)
{
  ///
  /// Configure the module given the Registry r
  ///
}

//......................................................................

void ParticleReport::Reset()
{
  ///
  /// (Document me!)
  ///
}

////////////////////////////////////////////////////////////////////////
