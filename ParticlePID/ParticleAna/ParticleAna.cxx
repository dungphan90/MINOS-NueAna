////////////////////////////////////////////////////////////////////////
/// $Id: ParticleAna.cxx,v 1.5 2009/06/30 16:45:12 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TROOT.h"

#include "NueAna/ParticlePID/ParticleAna/ParticleAna.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro

#include "NueAna/ParticlePID/ParticleAna/PRecord.h"
#include "NueAna/ParticlePID/ParticleAna/PRecordAna.h"


#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObject.h"


#include "NueAna/ParticlePID/ParticleFinder/ParticleBeamMonAna.h"

#include "DataUtil/MCInfo.h"

JOBMODULE(ParticleAna, "ParticleAna",
          "does the ana for ParticleFinder");
CVSID("$Id: ParticleAna.cxx,v 1.5 2009/06/30 16:45:12 scavan Exp $");

//SKZPWeightCalculator *ParticleAna::skzpCalc=0;


//......................................................................

ParticleAna::ParticleAna() 
{

	lastrun=-1;
DoMRCC=0;
outObjectName="normal";
//	pot = new POT();
  //if(!skzpCalc)skzpCalc=new SKZPWeightCalculator("DetXs",true);
	
  ///
  /// (Document me!)
  ///
}
//......................................................................

ParticleAna::~ParticleAna()
{
  ///
  /// (Document me!)
  ///


//  if(pot)delete pot;
//  pot=0;
  
  //we can reususe it!
  //if(skzpCalc)delete skzpCalc;
  //skzpCalc=0;

}

//......................................................................

void ParticleAna::BeginJob()
{
  ///
  /// (Document me!)
  ///
}

//......................................................................

void ParticleAna::EndJob()
{
  ///
  /// (Document me!)
  ///

/*   MSG("ParticleAna",Msg::kInfo)<<"Number of POT in this job: "<<pot->pot<<endl;
        
   TDirectory *savedir = gDirectory;
        
   TFile *fpf = dynamic_cast<TFile *>(gROOT->GetListOfFiles()->FindObject(kPOTTreeName.c_str()));
   if(fpf){
     fpf->cd();
     TTree *pottree = new TTree("pottree","pottree");
     pottree->Branch("ParticlePOT",&pot);
     pottree->Fill();
     pottree->Write();
     savedir->cd();
   }
   else{
     MSG("ParticleAna",Msg::kError)<<"Could not find the file to write the pottree to("<<kPOTTreeName<<"), there will be no pottree"<<endl;
   }
*/


}

//......................................................................

JobCResult ParticleAna::Reco(MomNavigator* mom)
{
  ///
  /// (Document me!)
  ///

//MSG("ParticleAna",Msg::kDebug) << "reco called\n";

  //std::vector<ParticleObjectHolder *>* hf;
  std::vector<TObject* > hft =( mom->GetFragmentList("ParticleObjectHolder",outObjectName.c_str()));


 
//MSG("ParticleAna",Msg::kDebug) << "size of found array " << hft.size()<<endl; 

 	ParticleBeamMon *bmon=0;

for(unsigned int s =0;s<hft.size();s++)
{


	ParticleObjectHolder *n=0;
	n=dynamic_cast<ParticleObjectHolder *>(hft[s]);


//MSG("ParticleAna",Msg::kDebug) << &hft[s]<<"\n";

//MSG("ParticleAna",Msg::kDebug) << "snarl "<<n->GetHeader().GetSnarl() <<"evt "<< n->GetHeader().GetEvent()<<"\n";
 	
	RecCandHeader ntphdr(n->GetHeader().GetVldContext(),n->GetHeader().GetRun(),
         				n->GetHeader().GetSubRun(),n->GetHeader().GetRunType(),n->GetHeader().GetErrorCode(),
         				n->GetHeader().GetSnarl(),n->GetHeader().GetTrigSrc(),n->GetHeader().GetTimeFrame(),
         				n->GetHeader().GetRemoteSpillType(),n->GetHeader().GetEvent());


	
 	//get beam info
 	if(s==0)
 	{
 	
 		RecCandHeader ntphdr(n->GetHeader().GetVldContext(),n->GetHeader().GetRun(),
         				n->GetHeader().GetSubRun(),n->GetHeader().GetRunType(),n->GetHeader().GetErrorCode(),
         				n->GetHeader().GetSnarl(),n->GetHeader().GetTrigSrc(),n->GetHeader().GetTimeFrame(),
         				n->GetHeader().GetRemoteSpillType(),-1);
 	
		bmon = new ParticleBeamMon(ntphdr);

/*	
		ParticleBeamMonAna bma;
		bma.ana(n,bmon);
		
		//FillPot(bmon);
		mom->AdoptFragment(bmon);			
 */	}
 	
 //	if(!bmon)bmon=(ParticleBeamMon *)mom->GetFragment("ParticleBeamMon");
 	
 	
	

	
	PRecord * r = new PRecord(ntphdr);


	PRecordAna pra;

	
	pra.ana(n,r,bmon);
	r->SetName(outObjectName.c_str());
	mom->AdoptFragment(r);	

	delete bmon;
	bmon=0;
}

		


  return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................

const Registry& ParticleAna::DefaultConfig() const
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
  //r.Set("POTTreeFileName","pottree.root");
  r.Set("DoMRCC",0);
  r.Set("OutObjectName","Normal");
  r.LockValues();


  return r;
}

//......................................................................

void ParticleAna::Config(const Registry& r)
{
  ///
  /// Configure the module given the Registry r
  ///
  
 
    const char* tmps;
    int tmpi;
 //   if(r.Get("POTTreeFileName",tmps)){kPOTTreeName=tmps;}
    if(r.Get("DoMRCC",tmpi)){DoMRCC=tmpi;}
    if(r.Get("OutObjectName",tmps)){outObjectName=tmps;}
}

//......................................................................

void ParticleAna::Reset()
{
  ///
  /// (Document me!)
  ///
}

////////////////////////////////////////////////////////////////////////


void ParticleAna::FillPot(ParticleBeamMon *bmon)
{
		//int ismc=bmon->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC;
		VldContext evt_vldc = bmon->GetHeader().GetVldContext();



/*
		if(!ismc)
		{
			pot->pot+=bmon->GetPot();
		}else{
		
			std::string relName = bmon->GetTitle();

  			std::string mcinfo = "";
  			if(bmon->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC){
     			mcinfo = "Carrot";
     	//		string temp = mchdr.geninfo.codename;
     	//		if(temp.size() != 0){   mcinfo = temp;  }
  			}
  					
			ReleaseType::Release_t rel = ReleaseType::MakeReleaseType(relName, mcinfo);
		
			
			//near for every snarl
			if(evt_vldc.GetDetector()==Detector::kNear){
              pot->pot+=MCInfo::GetMCPoT(Detector::kNear, bmon->beamtype, rel);
              
      //        printf("near %d %d\n",bmon->beamtype,rel);
            }
			
			//far only once per run
			if(evt_vldc.GetDetector()==Detector::kFar && 
             	bmon->GetHeader().GetRun()!=lastrun){
            		pot->pot+=MCInfo::GetMCPoT(Detector::kFar, bmon->beamtype, rel);
          	}
   
		}
		
		/ *
	 / *	
		pot->nsnarls++;
		if(bmon->GetHeader().GetRun()!=lastrun)
		{
			pot->nruns++;
			


			pot->beamtype = bmon->beamtype;
 			
			
			
		}
		
	*/	
		lastrun=bmon->GetHeader().GetRun();
		

}



