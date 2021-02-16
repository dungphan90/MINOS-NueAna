////////////////////////////////////////////////////////////////////////
// $Id: NueRead.cxx,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TFile.h"
#include "NueAna/NueRead.h"
#include "NueAna/NueRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro

JOBMODULE(NueRead, "NueRead",
          "Reads in ana_nue ntuples");
CVSID("$Id:");
//......................................................................

NueRead::NueRead():
   counter(0),
   passcounter(0)
{
   f=new TFile("testfile.root","RECREATE");
   f->cd();
   hpar0=new TH1F("hpar0","hpar0",100,0,100);
}

//......................................................................

NueRead::~NueRead()
{
   if(f){
      if(f->IsOpen()){
	 f->Write();
	 f->Close();
      }
   }

//======================================================================
// FILL_IN: [Document your code!!]
//======================================================================
}

//......................................................................

JobCResult NueRead::Ana(const MomNavigator* mom)
{
   //get all NueRecords from mom 
   //may have more than one per go since mom reads in a snarl's worth of data
   //so, this is a little more complicated than just asking for a NueRecord
   TObject *obj=0;
   TIter objiter = mom->FragmentIter();
   while((obj=objiter.Next())){
      NueRecord *nr = static_cast<NueRecord *>(obj);
      if(nr){
	 MSG("NueRead",Msg::kDebug)<<"Found a NueRecord in MOM"<<endl;
      }
      else{
	 MSG("NueRead",Msg::kError)<<"Didn't find a NueRecord in MOM"<<endl;
	 counter++;
	 continue;
      }
      MSG("NueRead",Msg::kDebug)<<"Found a NueRecord "<<nr<<endl;
      counter++;
      
      //fill histogram
      hpar0->Fill(nr->shwfit.par_a);
      passcounter++;
   }
   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
void NueRead::EndJob()
{
   MSG("NueRead",Msg::kInfo)<<"Counter "<<counter<<" passcounter "<<passcounter<<endl;
   MSG("NueRead",Msg::kInfo)<<"Entries in histogram "<<hpar0->GetEntries()<<endl;
}
