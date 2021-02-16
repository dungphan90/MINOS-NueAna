////////////////////////////////////////////////////////////////////////
// $Id: NueReweight.cxx,v 1.8 2008/11/19 18:22:51 rhatcher Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle
////////////////////////////////////////////////////////////////////////
#include <iostream>
#include <fstream>
#include <cmath>
#include "TMath.h"
#include "NueAna/NueRecord.h"
#include "NueAna/Reweight/NueReweight.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "Conventions/Detector.h"
#include "MCReweight/NeugenWeightCalculator.h"
#include "MCReweight/MCReweight.h"
#include "MCReweight/ReweightHelpers.h"
#include "Registry/Registry.h"

JOBMODULE(NueReweight, "NueReweight",
          "nue reweighting");
CVSID("$Id: NueReweight.cxx,v 1.8 2008/11/19 18:22:51 rhatcher Exp $");
ClassImp(NueReweight)

//......................................................................

NueReweight::NueReweight():
   mcr(MCReweight::Instance()),
   counter(0),
   rand(0),
   firstevent(true)
{
   MSG("NueReweight",Msg::kDebug)<<"Making a NueReweight"<<endl;
   rw = new NueRW();

   // make a WeightCalculator and add it to the MCReweight singleton
   NeugenWeightCalculator *n=new NeugenWeightCalculator();
   mcr.AddWeightCalculator(n);
}

//......................................................................

NueReweight::~NueReweight()
{
   MSG("NueReweight",Msg::kDebug)<<"Killing a NueReweight"<<endl;
   if(rw){
      delete rw;
      rw=0;
   }
}

//......................................................................

void NueReweight::EndJob()
{
   MSG("NueReweight",Msg::kDebug)<<"End Job"<<endl;
}

//......................................................................

JobCResult NueReweight::Reco(MomNavigator* mom)
{
  //  if(counter%1000==0){
  //      MSG("NueReweight",Msg::kInfo)<<"On entry "<<counter<<endl;
  //  }
  rw->randrow=rand;
   counter++;
   TObject *obj=0;
   TIter objiter = mom->FragmentIter();
   int objcount=0;
   while((obj=objiter.Next())){
      const char *cn=obj->ClassName();
      //      MSG("NueReweight",Msg::kDebug)<<"Found a "<<cn
      //	      		    <<". Objcount "<<objcount<<endl;
      if(strcmp(cn,"NueRecord")!=0){
	continue;
      }
      objcount++;
      NueRecord *nr = dynamic_cast<NueRecord *>(obj);
      if(firstevent){
	MSG("NueReweight",Msg::kDebug)<<"this is the first event"<<endl;
	rw->Reset();
	rw->fRun = nr->GetHeader().GetRun();
	rw->fSubRun = nr->GetHeader().GetSubRun();
	rw->fDet=nr->GetHeader().GetVldContext().GetDetector();
	int iflv = (int)((rw->fRun%100000-rw->fRun%1000)/10000);
	MSG("NueReweight",Msg::kDebug)<<"run is "<<rw->fRun	  
				      <<" run%100000 is "<<rw->fRun%100000
				      <<" run%1000 is "<<rw->fRun%1000<<endl;
	
	MSG("NueReweight",Msg::kDebug)<<"iflv is "<<iflv<<endl;
	if(iflv==0){
	  rw->fFileType = NueRW::kBEAM;
	}
	else if(iflv==1){
	  rw->fFileType = NueRW::kNUE;
	}
	else if(iflv==3){
	  rw->fFileType = NueRW::kTAU;
	}
	else{
	  rw->fFileType = NueRW::kUnknown;
	}
	rw->randrow=rand;
	ReadRandom();      
      // to calculate a reweight factor, two registry objects are used
      // the first contains the list of parameters in the WeightCalculator to
      // vary 
      // the list of parameters that can be varied can be seen by looking in the
      // appropriate WeightCalculator::Config() methods

	rwtconfig.UnLockValues();
	rwtconfig.UnLockKeys();
	rwtconfig.Set("neugen:use_scale_factors",1);
	rwtconfig.Set("neugen:ma_qe",rw->qel_ma);
	rwtconfig.Set("neugen:ma_res",rw->res_ma);
	//      rwtconfig.Set("neugen:ma_coh",rw->coh_ma);
	//      rwtconfig.Set("neugen:qel_fa0",rw->qel_fa0);
	//      rwtconfig.Set("neugen:qek_eta",rw->qel_eta);
	//      rwtconfig.Set("neugen:res_qe",rw->res_omega);
	//      rwtconfig.Set("neugen:res_z",rw->res_z);
	//      rwtconfig.Set("neugen:coh_r0",rw->coh_r0);
	//      rwtconfig.Set("neugen:coh_rei",rw->coh_rei);
	//      rwtconfig.Set("neugen:kno_a1",rw->kno_a1);
	//      rwtconfig.Set("neugen:kno_a2",rw->kno_a2);
	//      rwtconfig.Set("neugen:kno_a3",rw->kno_a3);
	//      rwtconfig.Set("neugen:kno_a4",rw->kno_a4);
	//      rwtconfig.Set("neugen:kno_b",rw->kno_b);
	rwtconfig.Set("neugen:kno_r112",rw->kno_r112);
	rwtconfig.Set("neugen:kno_r122",rw->kno_r122);
	rwtconfig.Set("neugen:kno_r132",rw->kno_r132);
	rwtconfig.Set("neugen:kno_r142",rw->kno_r142);
	rwtconfig.Set("neugen:kno_r113",rw->kno_r113);
	rwtconfig.Set("neugen:kno_r123",rw->kno_r123);
	rwtconfig.Set("neugen:kno_r133",rw->kno_r133);
	rwtconfig.Set("neugen:kno_r143",rw->kno_r143);
	rwtconfig.LockValues();
	rwtconfig.LockKeys();

	firstevent=false;
	MSG("NueReweight",Msg::kDebug)<<"Set rw->randrow to "<<rw->randrow<<endl;
      }
      if(objcount==1){
	rw->nsnarls++;
      }
      rw->nevents++;   
      rw->nacc++;
      //do reweighting--

      //make event regsitry
      Registry evreg;
      ReweightHelpers::EventRegistryFilla(&(nr->mctrue),evreg);

      // now all that's left is to compute the weight as follows:
      float reweight = mcr.ComputeWeight(&evreg,&rwtconfig);

      MSG("NueReweight",Msg::kDebug)<<"Reweight is "<<reweight<<endl;
       
      MSG("NueReweight",Msg::kDebug)<<"*********************"<<endl;      
      //      evreg.Print();
      MSG("NueReweight",Msg::kDebug)<<"*********************"<<endl;      
      //      rwtconfig.Print();


      //oscillate
      float oscprob = 1.;
      if(nr->GetHeader().GetVldContext().GetDetector()==Detector::kFar){
	float theta23 = TMath::ASin(sqrt(rw->ss2th))/2.;
	oscprob = NueConvention::Oscillate(&(nr->mctrue),735.,
					  rw->dm2,theta23,rw->UE32);
	MSG("NueReweight",Msg::kDebug)<<"dm2 is "<<rw->dm2<<" theta23 is "<<theta23
				      <<" ue32 is "<<rw->UE32<<endl;
      }
      MSG("NueReweight",Msg::kDebug)<<"oscprob is "<<oscprob<<endl;
      int ebin=rw->FindEBin(nr->srevent.energyGeV);
      float scale=reweight*oscprob;

      int iaction = nr->mctrue.interactionType;
      int inu = nr->mctrue.nuFlavor;
      int inunoosc = nr->mctrue.nonOscNuFlavor;
      
      if(iaction==0){
	 //neutral current;
	rw->nnc+=scale;
	rw->nbg+=scale;
	rw->nncE[ebin]+=scale;
	rw->nbgE[ebin]+=scale;
      }
      else{
	if(abs(inu)==12){
	  if(abs(inunoosc)==12){
	     rw->nnueb+=scale;
	     rw->nbg+=scale;
	     rw->nnuebE[ebin]+=scale;
	     rw->nbgE[ebin]+=scale;
	  }
	  else if(abs(inunoosc)==14){
	    rw->nsig+=scale;
	    rw->nsigE[ebin]+=scale;
	  }
	}
	else if(abs(inu)==14){
	  rw->nnumu+=scale;
	  rw->nbg+=scale;
	  rw->nnumuE[ebin]+=scale;
	  rw->nbgE[ebin]+=scale;
	}
	else if(abs(inu)==16){
	  rw->nnutau+=scale;
	  rw->nbg+=scale;
	  rw->nnutauE[ebin]+=scale;
	  rw->nbgE[ebin]+=scale;
	}
      }
   }          
//   MSG("NueReweight",Msg::kDebug)<<"InReco"<<endl;
   

   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................

const Registry& NueReweight::DefaultConfig() const
{
  static Registry r; // Default configuration for module

  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  // Set values in configuration
  r.UnLockValues();
  r.Set("RandNumber", 0);
  r.LockValues();

  return r;
}

//......................................................................

void NueReweight::Config(const Registry& r)
{
  int    tmpi;

  if (r.Get("RandNumber",tmpi)) { rand = tmpi; }
}


//......................................................................

int NueReweight::ReadRandom()
{
   const Int_t NCOLS=25;

   string fname(getenv("SRT_PRIVATE_CONTEXT"));
   //   fname+=("/NueAna/data/randnumbers.txt");
   //   fname+=("/NueAna/data/uniformnumbers.txt");
   fname+=("/NueAna/data/biguniformnumbers.txt");
   ifstream in(fname.c_str());
   if(!in){
      MSG("NueReweight",Msg::kError)<<"Could not open random number file in priviate area"<<endl
				    <<"Tried: "<<fname<<endl;
      fname=getenv("SRT_PUBLIC_CONTEXT");
      //      fname+=("/NueReweight/data/randnumbers.txt");
      fname+=("/NueAna/data/biguniformnumbers.txt");
      MSG("NueReweight",Msg::kError)<<"Will now try: "<<fname<<endl;
      in.open(fname.c_str());
      if(!in){
	 MSG("NueReweight",Msg::kError)<<"Could not open a random number file"<<endl;
	 return -1;
      }
   }

   float rands[NCOLS];
   int counter=0;
   in>>rands[0];
   while(!in.eof()){
      for(int i=1;i<NCOLS;i++){
	 in>>rands[i];
      }
      if(counter==rw->randrow){
	 break;
      }
      counter++;
      in>>rands[0];
   }

   if(counter<rw->randrow){
      MSG("NueReweight",Msg::kError)<<"Couldn't get row "<<rw->randrow<<" from "<<fname<<endl;
      return -1;
   }
   
   rw->qel_ma=fabs(rands[0]);
   rw->res_ma=fabs(rands[1]);
   rw->coh_ma=fabs(rands[2]);
   rw->qel_fa0=-1.*fabs(rands[3]);
   rw->qel_eta=fabs(rands[4]);
   rw->res_omega=fabs(rands[5]);
   rw->res_z=fabs(rands[6]);
   rw->coh_r0=fabs(rands[7]);
   rw->coh_rei=fabs(rands[8]);
   rw->kno_a1=fabs(rands[9]);
   rw->kno_a2=fabs(rands[10]);
   rw->kno_a3=fabs(rands[11]);
   rw->kno_a4=fabs(rands[12]);
   rw->kno_b=fabs(rands[13]);

   //use same scale for all kno_r* variables
   rw->kno_r112=fabs(rands[14]);
   rw->kno_r122=fabs(rands[14]);
   rw->kno_r132=fabs(rands[14]);
   rw->kno_r142=fabs(rands[14]);   
   rw->kno_r113=fabs(rands[14]);
   rw->kno_r123=fabs(rands[14]);
   rw->kno_r133=fabs(rands[14]);
   rw->kno_r143=fabs(rands[14]);
   rw->dm2=fabs(rands[22]);
   if(rands[23]>1){
     rw->ss2th=1.-(rands[23]-1.);
   }
   else{
     rw->ss2th=fabs(rands[23]);
   }
   rw->UE32=fabs(rands[24]);
   return 0;
}

//......................................................................

void NueReweight::Reset()
{
   counter=0;
   rand=0;
   firstevent=true;
   rw->Reset();
   rwtconfig.UnLockKeys();
   rwtconfig.UnLockValues();
   rwtconfig.Clear();
   rwtconfig.LockKeys();
   rwtconfig.LockValues();
   return;
}


////////////////////////////////////////////////////////////////////////

