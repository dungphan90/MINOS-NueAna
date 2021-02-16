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
#include "NtpStTrimmer.h"
#include "NueAna/NueAnalysisCuts.h"
#include "CandNtupleSR/NtpSRRecord.h"


const int CALend=110;
const int SM1end=230;
const int SM2end=484;

JOBMODULE(NtpStTrimmer, "NtpStTrimmer",
          "Reduce the file size of AnaNue files by filtering out events");
CVSID("$Id: NtpStTrimmer.cxx,v 1.3 2008/11/19 18:22:51 rhatcher Exp $");
//......................................................................

NtpStTrimmer::NtpStTrimmer():
  counter(0),
  kOutputFile("TrimmedOut.root"),
  kReWeight(0),
  kTheta23(1.0),
  kUe3Square(0.01),
  kDeltaMSquare(0.0025)
{}

//......................................................................

NtpStTrimmer::~NtpStTrimmer()
{}

//......................................................................
void NtpStTrimmer::BeginJob()
{
  kept = 0;
//  fCuts.Reset();
  if(counter%1000==0){
    cout<<"On entry "<<counter<<endl;
  }
}


JobCResult NtpStTrimmer::Reco(MomNavigator* mom)
{
   bool foundST=false;
                                                                                
   VldContext vc;
   //first ask mom for a NtpStRecord
   NtpStRecord *str=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord","Primary"));
   if(str){
     foundST=true;
     vc=str->GetHeader().GetVldContext();
   }

   cout<<"whats up"<<endl; 
   if(!foundST){
       MSG("NtpStTrimmer",Msg::kError)<<"This code ONLY runs on NtpSt"<<endl;
   }

   Int_t evtn = 0;
   if(foundST){
     evtn=str->evthdr.nevent;
   }

   str->SetName("NtpSt_junk");
   
   for(int i = 0; i < evtn; i++)
   {
      cout<<"looking at event "<<i<<endl;
      fCuts.SetInfoObject(i, str);
      if(PassesCuts()) {
          MSG("NtpStTrimmer",Msg::kDebug)<<"Excellent a NueRecord survives"
                    <<" Snarl "<<str->GetHeader().GetSnarl()
                    <<" Event "<<i<<endl;
                                                                                
          str->SetName("NtpSt_trim");
          kept++;
      }
   }
   if(counter%1000==0){
      cout<<"On entry "<<counter<<endl;
   }
   counter++;
   
   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
void NtpStTrimmer::EndJob()
{
  cout<<counter<<" events processed and "<<kept<<" were output"<<endl;

  //fCuts.Report();
//Now i have to output the tree
  // Here is where all of the writeout work will be done, first the tree and then the histos

  //If we are oscillating the files, then the effective exposure is only one third per file
  //    if you are not using an equal number of files from each type I accept no responsibility
  //    for the nature of your results

}

const Registry& NtpStTrimmer::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("NtpStTrimmer",Msg::kDebug)<<"In Trimmer::DefaultConfig"<<endl;

  static Registry r = fCuts.DefaultConfig();
 
  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());
                                                                                
  // Set values in configuration
  r.UnLockValues();
  r.Set("OutputFile", "TrimmedNtuple.root");

  r.Set("DeltaMSquare", 0.0025);
  r.Set("Theta23", TMath::Pi()/4);
  r.Set("Ue3Square", 0.01);
  r.Set("ReWeight", 0);
  
  r.LockValues();
                                                                                
  return r;
}

void NtpStTrimmer::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("NtpStTrimmer",Msg::kDebug)<<"In Trimmer::Config"<<endl;
  
  fCuts.Config(r);
                                                                               
  //  const char* tmps;
  int imps;
  if(r.Get("ReWeight", imps)) {kReWeight = imps;}
  
  double fmps;
  if(r.Get("DeltaMSquare", fmps)) {kDeltaMSquare = fmps;}
  if(r.Get("Ue3Square", fmps))    {kUe3Square = fmps;}
  if(r.Get("Theta23", fmps))  {kTheta23 = fmps;}

  const char* tmps;
  if(r.Get("OutputFile", tmps)) {kOutputFile = tmps;}
  
}

bool NtpStTrimmer::PassesCuts()
{
    bool passes = true;   
    if(!fCuts.PassesFiducialVolume())    passes = false;
    if(!fCuts.PassesFullContainment())  passes = false;   
     
    if(!fCuts.PassesAllCuts()) passes = false;
//    if(!fCuts.PassesFileCut())  passes = false;

     // Cut on min Total pulse height per prong (sigcor)
//     if((nr->srevent.planes<1.05||nr->srevent.planes>16) ||
//         (nr->hitcalc.fHitLongEnergy<30||nr->hitcalc.fHitLongEnergy>5000)||
//         (nr->shwfit.uv_molrad_vert<0||nr->shwfit.uv_molrad_vert>5.51)||
//        (nr->shwfit.par_b<0.315||nr->shwfit.par_b>1.47)) passes=false;

   return passes;
}

