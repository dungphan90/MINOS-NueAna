////////////////////////////////////////////////////////////////////////
// $Id: NueReadwPID.cxx,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TFile.h"
#include "NueAna/NueReadwPID.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NuePID.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "HistMan/HistMan.h"
JOBMODULE(NueReadwPID, "NueReadwPID",
          "Reads in ana_nue ntuple and pid ntuple");
CVSID("$Id:");
//......................................................................

NueReadwPID::NueReadwPID()
 {}

//......................................................................

NueReadwPID::~NueReadwPID()
{}

//......................................................................

JobCResult NueReadwPID::Ana(const MomNavigator* mom)
{
   //get all NueRecords from mom 
   //may have more than one per go since mom reads in a snarl's worth of data
   //so, this is a little more complicated than just asking for a NueRecord
   MSG("NueReadwPID",Msg::kDebug)<<"***********IN ANA*****************"<<endl;
   
   TObject *obj=0;
   TIter objiter = mom->FragmentIter();
   vector<NueRecord *> vr;
   vector<NuePID *> vpid;
   while((obj=objiter.Next())){
      const char *cn=obj->ClassName();
      if(strcmp(cn,"NueRecord")==0){
	 NueRecord *nr = dynamic_cast<NueRecord *>(obj);
	 MSG("NueReadwPID",Msg::kDebug)<<"Found a NueRecord in MOM"
				       <<" Snarl "<<nr->GetHeader().GetSnarl()
				       <<" Event "<<nr->GetHeader().GetEventNo()<<endl;
	 vr.push_back(nr);
      }
      else if(strcmp(cn,"NuePID")==0){
	 NuePID *npid  = dynamic_cast<NuePID *>(obj);
	 MSG("NueReadwPID",Msg::kDebug)<<"Found a NuePID in MOM"
				       <<" Snarl "<<npid->GetHeader().GetSnarl()
				       <<" Event "<<npid->GetHeader().GetEventNo()<<endl;
	 vpid.push_back(npid);
      }
      else{
	 continue;
      }
   }


   //use HistMan to plot something for events that pass/fail

   static HistMan *hm = new HistMan("nueana");
   static TH1F *h1=0;
   static TH1F *h2=0;
   static TH1F *h3=0;
   static TH1F *h4=0;
   static TH1F *h5=0;
   static TH1F *h6=0;
   static TH1F *h7=0;
   static TH1F *h8=0;
   static TH1F *h9=0;

   h1=hm->Book<TH1F>("par_a_pass","par_a_pass",200,0,50);
   h2=hm->Book<TH1F>("par_b_pass","par_b_pass",200,0,50);
   h3=hm->Book<TH1F>("par_e0_pass","par_e0_pass",200,0,50);
   h4=hm->Book<TH1F>("par_a_fail","par_a_fail",200,0,50);
   h5=hm->Book<TH1F>("par_b_fail","par_b_fail",200,0,50);
   h6=hm->Book<TH1F>("par_e0_fail","par_e0_fail",200,0,50);
   h7=hm->Book<TH1F>("par_a_hmmm","par_a_hmmm",200,0,50);
   h8=hm->Book<TH1F>("par_b_hmmm","par_b_hmmm",200,0,50);
   h9=hm->Book<TH1F>("par_e0_hmmm","par_e0_hmmm",200,0,50);

   
   //so, mom will match up snarls for us,
   //but, we have to match up events for ourselves.
   for(unsigned int i=0;i<vr.size();i++){
      int event = vr[i]->GetHeader().GetEventNo();
      bool foundmatch=false;
      for(unsigned int j=0;j<vpid.size();j++){
	 if(vpid[j]->GetHeader().GetEventNo()==event){
	    int pass = vpid[j]->IsNue;
	    MSG("NueReadwPID",Msg::kDebug)<<"Found match!"<<endl
					  <<" Record snarl: "<<vr[i]->GetHeader().GetSnarl()
					  <<" event: "<<vr[i]->GetHeader().GetEventNo()<<endl
					  <<" PID snarl: "<<vpid[j]->GetHeader().GetSnarl()
					  <<" event: "<<vpid[j]->GetHeader().GetEventNo()<<endl;
	    MSG("NueReadwPID",Msg::kDebug)<<"pass "<<pass<<" j "<<j<<" "<<vpid[j]->IsNue<<endl;

	    if(pass==1){
	       hm->Fill1d("par_a_pass",vr[i]->shwfit.par_a);
	       hm->Fill1d("par_b_pass",vr[i]->shwfit.par_b);
	       hm->Fill1d("par_e0_pass",vr[i]->shwfit.par_e0);
	    }
	    else if(pass==-1){
	       hm->Fill1d("par_a_fail",vr[i]->shwfit.par_a);
	       hm->Fill1d("par_b_fail",vr[i]->shwfit.par_b);
	       hm->Fill1d("par_e0_fail",vr[i]->shwfit.par_e0);
	    }
	    else{
	       hm->Fill1d("par_a_hmmm",vr[i]->shwfit.par_a);
	       hm->Fill1d("par_b_hmmm",vr[i]->shwfit.par_b);
	       hm->Fill1d("par_e0_hmmm",vr[i]->shwfit.par_e0);
	    }
	 }
	 //delete pid from vector so we don't loop over it next time
	 vector<NuePID *>::iterator vi(&(vpid[j]));
	 vpid.erase(vi);
	 foundmatch=true;
	 break;
      }
      if(!foundmatch){
	 MSG("NueReadwPID",Msg::kError)<<"Could not find PID match for"
				       <<" Snarl "<<vr[i]->GetHeader().GetSnarl()
				       <<" Event "<<vr[i]->GetHeader().GetEventNo()<<endl;
      }
   }

   MSG("NueReadwPID",Msg::kDebug)<<"**********************************"<<endl;

   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
