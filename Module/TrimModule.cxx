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
#include "TrimModule.h"
#include "NueAna/NueAnalysisCuts.h"

const int CALend=110;
const int SM1end=230;
const int SM2end=484;

JOBMODULE(TrimModule, "TrimModule",
          "Reduce the file size of AnaNue files by filtering out events");
CVSID("$Id: TrimModule.cxx,v 1.2 2008/11/19 18:22:51 rhatcher Exp $");
//......................................................................

TrimModule::TrimModule():
  counter(0),
  kOutputFile("TrimmedOut.root"),
  kReWeight(0),
  kTheta23(1.0),
  kUe3Square(0.01),
  kDeltaMSquare(0.0025)
{}

//......................................................................

TrimModule::~TrimModule()
{}

//......................................................................
void TrimModule::BeginJob()
{
  kept = 0;
//  fCuts.Reset();
  if(counter%1000==0){
    cout<<"On entry "<<counter<<endl;
  }
}


JobCResult TrimModule::Reco(MomNavigator* mom)
{
   //get all NueRecords from mom 
   //may have more than one per go since mom reads in a snarl's worth of data
   //so, this is a little more complicated than just asking for a NueRecord
   TObject *obj=0;
//  static Float_t total_pot = 0;

   vector<NueRecord *> records;

   TIter objiter = mom->FragmentIter();
   while((obj=objiter.Next())){
      NueRecord *nr = dynamic_cast<NueRecord *>(obj);
      if(nr){
	 MSG("TrimModule",Msg::kDebug)<<"Found a NueRecord in MOM"<<endl;
      }
      else{
	 MSG("TrimModule",Msg::kError)<<"Didn't find a NueRecord in MOM"<<endl;
	 continue;
      }
      nr->SetName("NueRecord_junk");
      MSG("TrimModule",Msg::kDebug)<<"Found a NueRecord in MOM"
                       <<" Snarl "<<nr->GetHeader().GetSnarl()
                       <<" Event "<<nr->GetHeader().GetEventNo()<<endl;
 
      if(nr->GetHeader().GetEventNo()<0){
	continue;
      }

      if(kReWeight)
      {
          int nuFlavor = nr->mctrue.nuFlavor;
          int  nonOsc = nr->mctrue.nonOscNuFlavor;
          float energy = nr->mctrue.nuEnergy;
                                                                                
          Float_t newWeight = NueConvention::Oscillate(nuFlavor, nonOsc, energy,
                                735, kDeltaMSquare, kTheta23, kUe3Square);
                                                                  
          nr->mctrue.Ue3Squared = kUe3Square;
          nr->mctrue.DeltamSquared23 = kDeltaMSquare;
          nr->mctrue.Theta23  = kTheta23;
          nr->mctrue.fOscProb = newWeight;
      }

      fCuts.SetInfoObject(nr);
      if(PassesCuts(nr) && PassesBeamCuts(nr)){
           MSG("TrimModule",Msg::kDebug)<<"Excellent a NueRecord survives"
                     <<" Snarl "<<nr->GetHeader().GetSnarl()
                     <<" Event "<<nr->GetHeader().GetEventNo()<<endl;

           nr->SetName("NueRecord_trim");
           kept++;
      }
      if(counter%1000==0){
         cout<<"On entry "<<counter<<endl;
      }
      counter++;
   }

   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
void TrimModule::EndJob()
{
  cout<<counter<<" events processed and "<<kept<<" were output"<<endl;

  //fCuts.Report();
//Now i have to output the tree
  // Here is where all of the writeout work will be done, first the tree and then the histos

  //If we are oscillating the files, then the effective exposure is only one third per file
  //    if you are not using an equal number of files from each type I accept no responsibility
  //    for the nature of your results

}

const Registry& TrimModule::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("TrimModule",Msg::kDebug)<<"In TrimModule::DefaultConfig"<<endl;

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

void TrimModule::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("TrimModule",Msg::kDebug)<<"In TrimModule::Config"<<endl;
  
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

bool TrimModule::PassesBeamCuts(NueRecord* nr)
{
  //bool passes = true;
  if(nr->GetHeader().GetVldContext().GetSimFlag()!=SimFlag::kData) return true;
  if(fCuts.PassesBeamCut()) return true;;

  return false;
}



bool TrimModule::PassesCuts(NueRecord* nr)
{
    bool passes = true;   
    if(nr == 0) return false;
    if(!fCuts.PassesFiducialVolume())    passes = false;
    if(!fCuts.PassesFullContainment())  passes = false;   
     
    if(!fCuts.PassesAllCuts())           passes = false;

     // Cut on min Total pulse height per prong (sigcor)
//     if((nr->srevent.planes<1.05||nr->srevent.planes>16) ||
//         (nr->hitcalc.fHitLongEnergy<30||nr->hitcalc.fHitLongEnergy>5000)||
//         (nr->shwfit.uv_molrad_vert<0||nr->shwfit.uv_molrad_vert>5.51)||
//        (nr->shwfit.par_b<0.315||nr->shwfit.par_b>1.47)) passes=false;

   return passes;
}

