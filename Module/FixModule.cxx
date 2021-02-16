////////////////////////////////////////////////////////////////////////
//
// FILL_IN: [Document your code!!]
//
// boehm@physics.harvard.edu
////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TFile.h"
#include "Conventions/Detector.h"
#include "NueAna/NueRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "HistMan/HistMan.h"
#include <fstream>
#include "AnalysisNtuples/ANtpDefaultValue.h"                                   
#include "TClass.h"
#include "TDataType.h"
#include "TDataMember.h"
#include "TRealData.h"
#include <string>
#include "FixModule.h"
#include "TChain.h"
#include "NueAna/NuePOT.h"
#include "TDirectory.h"
#include "TROOT.h"


#include "NueAna/Module/SetKNNModule.h"
#include "PhysicsNtuple/Handle.h"
#include "PhysicsNtuple/Factory.h"

const int CALend=110;
const int SM1end=230;
const int SM2end=484;

JOBMODULE(FixModule, "FixModule",
          "Reduce the file size of AnaNue files by filtering out events");
CVSID("$Id: FixModule.cxx,v 1.2 2008/11/19 18:33:13 rhatcher Exp $");
//......................................................................

FixModule::FixModule():
  counter(0)
{}

//......................................................................

FixModule::~FixModule()
{}

//......................................................................
void FixModule::BeginJob()
{
  kept = 0;
//  fCuts.Reset();
  if(counter%1000==0){
    cout<<"On entry "<<counter<<endl;
  }
}


JobCResult FixModule::Reco(MomNavigator* mom)
{
   //get all NueRecords from mom 
   //may have more than one per go since mom reads in a snarl's worth of data
   //so, this is a little more complicated than just asking for a NueRecord
   TObject *obj=0;
//  static Float_t total_pot = 0;

   vector<NueRecord *> records;

   Anp::Handle<StorekNNData> data = Anp::Factory<StorekNNData>::Instance().Get("kNNData");
   if(!data.valid())
   {
      MAXMSG("NueModule", Msg::kError, 1)
           << "NueModule::Reco - Handle<StorekNNData> is invalid, assuming no Rustem variable to run"
           << endl;
   }

   bool foundOne = false;
   TIter objiter = mom->FragmentIter();
   while((obj=objiter.Next())){
      NueRecord *nr = dynamic_cast<NueRecord *>(obj);
      if(nr){
         foundOne = true;
	 MSG("FixModule",Msg::kDebug)<<"Found a NueRecord in MOM"<<endl;
      }
      else{
	 MSG("FixModule",Msg::kDebug)<<"Didn't find a NueRecord in MOM"<<endl;
	 continue;
      }
      MSG("FixModule",Msg::kDebug)<<"Found a NueRecord in MOM"
                       <<" Snarl "<<nr->GetHeader().GetSnarl()
                       <<" Event "<<nr->GetHeader().GetEventNo()<<endl;
 

      if(nr->mri.orig_event > -10){
        float knn_pid;

        data -> SetPrefix("OldSNTP");
        data -> Get(nr->mri.orig_event, "knn_pid", knn_pid);
        nr->mri.orig_roCCPID = knn_pid; 
      }

      if(counter%1000==0){
         cout<<"On entry "<<counter<<endl;
      }
      counter++;
   }

   data->Clear();
   if(!foundOne) return JobCResult::kFailed;

   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
void FixModule::EndJob()
{
  cout<<counter<<" events processed and "<<kept<<" were output"<<endl;

  TChain* ch = new TChain("pottree");
  ch->Add(fInFile.c_str());
  
  NuePOT* pt = new NuePOT();
  
  ch->SetBranchAddress("NuePOT", &pt);
  ch->GetEntry(0);  

  TDirectory *savedir = gDirectory;

  cout<<"Trying to access: "<<fOutFile<<endl;

  TFile *fpf = dynamic_cast<TFile *>(gROOT->GetListOfFiles()->FindObject(fOutFile.c_str()));
  if(fpf){
     cout<<"Got the file and rolling"<<endl;
     fpf->cd();
     TTree *pottree = new TTree("pottree","pottree");
     pottree->Branch("NuePOT",&pt);
     pottree->Fill();
     pottree->Write();
     savedir->cd();
   }else{
    cout<<"failed!"<<endl;
  }
}

const Registry& FixModule::DefaultConfig() const
{
    static Registry r;
    r.UnLockValues();
    r.Set("InFileName","blank.root");
    r.Set("FileName","pottree.root");
    r.LockValues();
    
    return r;
}

void FixModule::Config(const Registry& r)
{
    const char* tmps;
    if(r.Get("InFileName",tmps)){ fInFile=tmps;}
    if(r.Get("FileName",tmps)){fOutFile=tmps;}
}


