////////////////////////////////////////////////////////////////////////
/// $Id: ParticleTruthMaker.cxx,v 1.1 2009/06/19 14:32:40 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////
#include "NueAna/ParticlePID/ParticleFinder/Truth/ParticleTruthMaker.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro

#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "CandNtupleSR/NtpSRPulseHeight.h"

#include "StandardNtuple/NtpStRecord.h"

#include "MCNtuple/NtpMCDigiScintHit.h"
#include "MCNtuple/NtpMCStdHep.h"

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObject.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle.h"



JOBMODULE(ParticleTruthMaker, "ParticleTruthMaker",
          "finds particles");

//......................................................................

ParticleTruthMaker::ParticleTruthMaker() :helper(0)  
{

helper = new ParticleTruthHelper();



  ///
  /// (Document me!)
  ///
}
//......................................................................

ParticleTruthMaker::~ParticleTruthMaker()
{
  ///
  /// (Document me!)
  ///

if(helper) delete helper;

}

//......................................................................

void ParticleTruthMaker::BeginJob()
{
  ///
  /// (Document me!)
  ///
}

//......................................................................

void ParticleTruthMaker::EndJob()
{
  ///
  /// (Document me!)
  ///
}

//......................................................................

JobCResult ParticleTruthMaker::Reco(MomNavigator* mom)
{
  ///
  /// (Document me!)
  ///

//MSG("ParticleTruthMaker",Msg::kDebug) << "reco called\n";

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
  

  //TClonesArray * digihit = str->digihit;

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

//for now multiple events are in the same TObjArray...

int ievt=0;
 // for(int ievt=0;ievt<evts.size();ievt++)
 // {

   //    std::vector<const NtpSRStrip* > stp = str->GetStrips(ievt);
//	TClonesArray * digihit = &( str->digihit[ievt]); 

	TClonesArray * digihit = str->digihit;


	TClonesArray * stdhep = str->stdhep;
	
//	TClonesArray *digihit = &dh[ievt];

	helper->Reset();


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

	if(n->particletruth->GetEntries()>0)return JobCResult::kPassed; ///continue; //dont rerun it




	printf("-event %d  with %d strips\n",ievt,evts[ievt]->nstrip);

	//Load processor with event strip data
	for(int i=0;i<digihit[ievt].GetEntries();i++)
	{
	//	helper->AddStrip(stp[evts[ievt]->stp[i]]->plane, stp[evts[ievt]->stp[i]]->strip, (stp[evts[ievt]->stp[i]]->ph0.sigcor+stp[evts[ievt]->stp[i]]->ph1.sigcor)/565.0, stp[evts[ievt]->stp[i]]->tpos, stp[evts[ievt]->stp[i]]->z);
		NtpMCDigiScintHit * dh =  ( NtpMCDigiScintHit *)digihit[ievt][i];
		

		int trkid = dh->trkId;
		trkid=trkid<0?-trkid:trkid; //abs it

		NtpMCStdHep * sh = (NtpMCStdHep *)stdhep[ievt][trkid];

		

		helper->AddStrip(sh->IdHEP, trkid, dh->plane,dh->strip, dh->dE, dh->planeview);


//		printf("%d %d %d %f\n", dh->pId, dh->plane,dh->strip, dh->dE);

	}

	helper->Process(*n);
	

	//process event, returns particleobject	

//dont need digihits anymore, so reset them
digihit->Clear();

  //}





  return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................

const Registry& ParticleTruthMaker::DefaultConfig() const
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

void ParticleTruthMaker::Config(const Registry& r)
{
  ///
  /// Configure the module given the Registry r
  ///
}

//......................................................................

void ParticleTruthMaker::Reset()
{
  ///
  /// (Document me!)
  ///
}

////////////////////////////////////////////////////////////////////////
