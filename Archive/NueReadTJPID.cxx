////////////////////////////////////////////////////////////////////////
// $Id: NueReadTJPID.cxx,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TFile.h"
#include "NueAna/NueReadTJPID.h"
#include "NueAna/NuePID.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "HistMan/HistMan.h"
#include "MCNtuple/NtpMCRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include "TruthHelperNtuple/NtpTHRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "NueAna/SntpHelpers.h"


JOBMODULE(NueReadTJPID, "NueReadTJPID",
          "Reads in ana_nue ntuple and pid ntuple");
CVSID("$Id:");
//......................................................................

NueReadTJPID::NueReadTJPID()
{}

//......................................................................

NueReadTJPID::~NueReadTJPID()
{}

//......................................................................

JobCResult NueReadTJPID::Ana(const MomNavigator* mom)
{
   //get all NueRecords from mom 
   //may have more than one per go since mom reads in a snarl's worth of data
   //so, this is a little more complicated than just asking for a NueRecord
   MSG("NueReadTJPID",Msg::kDebug)<<"***********IN ANA*****************"<<endl;
   
   TObject *obj=0;
   TIter objiter = mom->FragmentIter();
   NtpSRRecord *sr=0;
   NtpMCRecord *mc=0;
   NtpTHRecord *th=0;
   vector<NuePID *> vpid;
   while((obj=objiter.Next())){
      const char *cn=obj->ClassName();
//      MSG("NueReadTJPID",Msg::kDebug)<<"Found a "<<cn<<"."<<endl;
      if(strcmp(cn,"NtpSRRecord")==0){
	 sr = dynamic_cast<NtpSRRecord *>(obj);
      }
      else if(strcmp(cn,"NtpMCRecord")==0){
	 mc = dynamic_cast<NtpMCRecord *>(obj);
      }
      else if(strcmp(cn,"NtpTHRecord")==0){
	 th = dynamic_cast<NtpTHRecord *>(obj);
      }
      else if(strcmp(cn,"NuePID")==0){
	 NuePID *npid  = dynamic_cast<NuePID *>(obj);
//	 MSG("NueReadTJPID",Msg::kDebug)<<"Found a NuePID in MOM"
//				       <<" Snarl "<<npid->GetHeader().GetSnarl()
//				       <<" Event "<<npid->GetHeader().GetEventNo()<<endl;
	 vpid.push_back(npid);
      }
      else{
	 continue;
      }
   }
   if(!sr){
      MSG("NueReadTJPID",Msg::kError)<<"No NtpSRRecord in Mom"<<endl;
      return JobCResult::kFailed;
   }
   
   if(!mc){
      MSG("NueReadTJPID",Msg::kError)<<"No NtpMCRecord in Mom"<<endl;
      return JobCResult::kFailed;
   }   
   if(!th){
      MSG("NueReadTJPID",Msg::kError)<<"No NtpTHRecord in Mom"<<endl;
      return JobCResult::kFailed;
   }
   
//   MSG("NueReadTJPID",Msg::kDebug)<<"SR snarl: "<<sr->GetHeader().GetSnarl()<<endl;

   if(vpid.size()==0){
      MSG("NueReadTJPID",Msg::kError)<<"No NuePID Records in Mom"<<endl;
      return JobCResult::kFailed;
   }
      
   //use HistMan to plot something for events that pass/fail

   static HistMan *hm = new HistMan("tjpid");
   static TH1F *htEall=0;
   static TH1F *htEsig=0;
   static TH1F *htEnumu=0;
   static TH1F *htEnutau=0;
   static TH1F *htEbe=0;
   static TH1F *htEnc=0;

   static TH1F *hphall=0;
   static TH1F *hphsig=0;
   static TH1F *hphnumu=0;
   static TH1F *hphnutau=0;
   static TH1F *hphbe=0;
   static TH1F *hphnc=0;

   static TH1F *htrklnall=0;
   static TH1F *htrklnsig=0;
   static TH1F *htrklnnumu=0;
   static TH1F *htrklnnutau=0;
   static TH1F *htrklnbe=0;
   static TH1F *htrklnnc=0;

   static TH1F *hshwlnall=0;
   static TH1F *hshwlnsig=0;
   static TH1F *hshwlnnumu=0;
   static TH1F *hshwlnnutau=0;
   static TH1F *hshwlnbe=0;
   static TH1F *hshwlnnc=0;



   htEall = hm->Book<TH1F>("htEall","htEall",50,0,50);
   htEsig = hm->Book<TH1F>("htEsig","htEsig",50,0,50);
   htEnumu = hm->Book<TH1F>("htEnumu","htEnumu",50,0,50);
   htEnutau = hm->Book<TH1F>("htEnutau","htEnutau",50,0,50);
   htEbe = hm->Book<TH1F>("htEbe","htEbe",50,0,50);
   htEnc = hm->Book<TH1F>("htEnc","htEnc",50,0,50);

   hphall = hm->Book<TH1F>("hphall","hphall",50,0,50);
   hphsig = hm->Book<TH1F>("hphsig","hphsig",50,0,50);
   hphnumu = hm->Book<TH1F>("hphnumu","hphnumu",50,0,50);
   hphnutau = hm->Book<TH1F>("hphnutau","hphnutau",50,0,50);
   hphbe = hm->Book<TH1F>("hphbe","hphbe",50,0,50);
   hphnc = hm->Book<TH1F>("hphnc","hphnc",50,0,50);

   htrklnall = hm->Book<TH1F>("htrklnall","htrklnall",50,0,50);
   htrklnsig = hm->Book<TH1F>("htrklnsig","htrklnsig",50,0,50);
   htrklnnumu = hm->Book<TH1F>("htrklnnumu","htrklnnumu",50,0,50);
   htrklnnutau = hm->Book<TH1F>("htrklnnutau","htrklnnutau",50,0,50);
   htrklnbe = hm->Book<TH1F>("htrklnbe","htrklnbe",50,0,50);
   htrklnnc = hm->Book<TH1F>("htrklnnc","htrklnnc",50,0,50);

   hshwlnall = hm->Book<TH1F>("hshwlnall","hshwlnall",50,0,50);
   hshwlnsig = hm->Book<TH1F>("hshwlnsig","hshwlnsig",50,0,50);
   hshwlnnumu = hm->Book<TH1F>("hshwlnnumu","hshwlnnumu",50,0,50);
   hshwlnnutau = hm->Book<TH1F>("hshwlnnutau","hshwlnnutau",50,0,50);
   hshwlnbe = hm->Book<TH1F>("hshwlnbe","hshwlnbe",50,0,50);
   hshwlnnc = hm->Book<TH1F>("hshwlnnc","hshwlnnc",50,0,50);


   
   //so, mom will match up snarls for us,
   //but, we have to match up events for ourselves.
   for(unsigned int i=0;i<vpid.size();i++){
      int evtno = vpid[i]->GetHeader().GetEventNo();
      int mcindex=0;
      if(evtno<0){
	 MSG("NueReadTJPID",Msg::kDebug)<<"can not get mctruth for event "<<evtno<<endl;
      }
      else{
	 MSG("NueReadTJPID",Msg::kDebug)<<"Trying to get mc index "<<evtno<<endl;
	 mcindex = SntpHelpers::GetEvent2MCIndex(evtno,th);
	 MSG("NueReadTJPID",Msg::kDebug)<<"got mc index "<<mcindex<<endl;
      }
      NtpMCTruth *mcth = SntpHelpers::GetMCTruth(mcindex,mc);
      if(mcth==0){
	 MSG("NueReadTJPID",Msg::kError)<<"can not get mctruth for event "<<evtno<<endl;
	 continue;
      }

      NtpSREvent *event = 0;
      if(evtno>=0){
	 event = SntpHelpers::GetEvent(evtno,sr);
      }
      //loop over tracks in this event, find longest
      int longtrack=0;
      if(event!=0){
	 for(int j=0;j<event->ntrack;j++){
	    int tindex = SntpHelpers::GetTrackIndex(j,event);
	    NtpSRTrack *track = SntpHelpers::GetTrack(tindex,sr);
	    if(longtrack<track->plane.n){
	       longtrack = track->plane.n;
	    }
	 }
      }
      //loop over showers in this event, find longest
      int longshower=0;
      if(event!=0){
	 for(int j=0;j<event->nshower;j++){
	    int sindex = SntpHelpers::GetShowerIndex(j,event);
	    NtpSRShower *shower = SntpHelpers::GetShower(sindex,sr);
	    if(longshower<shower->plane.n){
	       longshower = shower->plane.n;
	    }
	 }
      }

      hm->Fill1d("htEall",mcth->p4neu[3]);
      if(event!=0){
	 hm->Fill1d("hphall",event->ph.pe);
      }
      hm->Fill1d("htrklnall",longtrack);
      hm->Fill1d("hshwlnall",longshower);
      
      int pass = vpid[i]->IsNue;
      if(pass==1){
	 int cls=0;
	 if(mcth->iaction==0){
	    cls=5;
	 }
	 else if(abs(mcth->inu)==12){
	    if(abs(mcth->inunoosc)==12){
	       cls=4;
	    }
	    else if(abs(mcth->inu)==14){
	       cls=1;
	    }
	 }
	 else if(abs(mcth->inu)==14){
	    cls=2;
	 }
	 else if(abs(mcth->inu)==16){
	    cls=3;
	 }
	 if(cls==1){
	    hm->Fill1d("htEsig",mcth->p4neu[3]);
	    if(event!=0) hm->Fill1d("hphsig",event->ph.pe);
	    hm->Fill1d("htrklnsig",longtrack);
	    hm->Fill1d("hshwlnsig",longshower);
	 }
	 else if(cls==2){
	    hm->Fill1d("htEnumu",mcth->p4neu[3]);
	    if(event!=0) hm->Fill1d("hphnumu",event->ph.pe);
	    hm->Fill1d("htrklnnumu",longtrack);
	    hm->Fill1d("hshwlnnumu",longshower);
	 }
	 else if(cls==3){
	    hm->Fill1d("htEnutau",mcth->p4neu[3]);
	    if(event!=0) hm->Fill1d("hphnutau",event->ph.pe);
	    hm->Fill1d("htrklnnutau",longtrack);
	    hm->Fill1d("hshwlnnutau",longshower);
	 }
	 else if(cls==4){
	    hm->Fill1d("htEbe",mcth->p4neu[3]);
	    if(event!=0) hm->Fill1d("hphbe",event->ph.pe);
	    hm->Fill1d("htrklnbe",longtrack);
	    hm->Fill1d("hshwlnbe",longshower);
	 }
	 else if(cls==5){
	    hm->Fill1d("htEnc",mcth->p4neu[3]);
	    if(event!=0) hm->Fill1d("hphnc",event->ph.pe);
	    hm->Fill1d("htrklnnc",longtrack);
	    hm->Fill1d("hshwlnnc",longshower);
	 }
      }
   }

   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
