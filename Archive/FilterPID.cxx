////////////////////////////////////////////////////////////////////////
// $Id: FilterPID.cxx,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TFile.h"
#include "NueAna/FilterPID.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NuePID.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "HistMan/HistMan.h"
JOBMODULE(FilterPID, "FilterPID",
          "Reads ana_nue ntuple and pid ntuple, writes out only those events that pass pid");
CVSID("$Id:");
//......................................................................

FilterPID::FilterPID()
 {}

//......................................................................

FilterPID::~FilterPID()
{}

//......................................................................

JobCResult FilterPID::Reco(MomNavigator* mom)
{
   //get all NueRecords from mom 
   //may have more than one per go since mom reads in a snarl's worth of data
   //so, this is a little more complicated than just asking for a NueRecord
   MSG("FilterPID",Msg::kDebug)<<"***********IN ANA*****************"<<endl;
   
   TObject *obj=0;
   TIter objiter = mom->FragmentIter();
   vector<NueRecord *> vr;
   vector<NuePID *> vpid;
   while((obj=objiter.Next())){
      const char *cn=obj->ClassName();
      MSG("FilterPID",Msg::kDebug)<<"Found a "<<cn<<endl;
      if(strcmp(cn,"NueRecord")==0){
	 NueRecord *nr = dynamic_cast<NueRecord *>(obj);
	 MSG("FilterPID",Msg::kDebug)<<"Found a NueRecord in MOM"
				       <<" Snarl "<<nr->GetHeader().GetSnarl()
				       <<" Event "<<nr->GetHeader().GetEventNo()<<endl;
	 vr.push_back(nr);
      }
      else if(strcmp(cn,"NuePID")==0){
	 NuePID *npid  = dynamic_cast<NuePID *>(obj);
	 MSG("FilterPID",Msg::kDebug)<<"Found a NuePID in MOM"
				       <<" Snarl "<<npid->GetHeader().GetSnarl()
				       <<" Event "<<npid->GetHeader().GetEventNo()<<endl;
	 vpid.push_back(npid);
      }
      else{
	 continue;
      }
   }

   
   //so, mom will match up snarls for us,
   //but, we have to match up events for ourselves.
   for(unsigned int i=0;i<vr.size();i++){
      int event = vr[i]->GetHeader().GetEventNo();
      bool foundmatch=false;
      for(unsigned int j=0;j<vpid.size();j++){
	 if(vpid[j]->GetHeader().GetEventNo()==event){
	    int pass = vpid[j]->IsNue;
	    MSG("FilterPID",Msg::kDebug)<<"Found match!"<<endl
					<<" Record snarl: "
					<<vr[i]->GetHeader().GetSnarl()
					<<" event: "
					<<vr[i]->GetHeader().GetEventNo()<<endl
					<<" PID snarl: "
					<<vpid[j]->GetHeader().GetSnarl()
					<<" event: "
					<<vpid[j]->GetHeader().GetEventNo()<<endl;
	    MSG("FilterPID",Msg::kDebug)<<"pass "<<pass<<" j "<<j
					<<" "<<vpid[j]->IsNue<<endl;
	    if(pass==1){
	      MSG("FilterPID",Msg::kDebug)<<"Gave to mom"<<endl;
	      NueRecord *passnr = new NueRecord(*(vr[i]));
	      //	      NueHeader passh(vr[i]->GetHeader().GetVldContext());
	      //	      NueRecord *passnr = new NueRecord(passh);
	      //	      cout<<"passnr "<<passnr<<" vr[i] "<<vr[i]<<endl;
	      //	      passnr->Print();
	      passnr->SetName("NueRecord_filt");
	      mom->AdoptFragment(passnr);
	      //	      mom->Print();
	    }

	 }
	 //delete pid from vector so we don't loop over it next time
	 vector<NuePID *>::iterator vi(&(vpid[j]));
	 vpid.erase(vi);
	 foundmatch=true;
	 break;
      }
      if(!foundmatch){
	 MSG("FilterPID",Msg::kError)<<"Could not find PID match for"
				       <<" Snarl "<<vr[i]->GetHeader().GetSnarl()
				       <<" Event "<<vr[i]->GetHeader().GetEventNo()<<endl;
      }
   }

   MSG("FilterPID",Msg::kDebug)<<"**********************************"<<endl;

   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
