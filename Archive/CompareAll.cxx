////////////////////////////////////////////////////////////////////////
//
// FILL_IN: [Document your code!!]
//
// boehm@physics.harvard.edu
////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TFile.h"
#include "Conventions/Detector.h"
#include "NueAna/CompareAll.h"
#include "NueAna/NueRecord.h"
#include "NueAna/EventFilter.h"
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
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include <string>

const int CALend=110;
const int SM1end=230;
const int SM2end=484;

JOBMODULE(CompareAll, "CompareAll",
          "Makes hists to All variables between detectors and truth classes");
CVSID("$Id:");
//......................................................................

CompareAll::CompareAll():
    counter(0),
   kHiPlaneTrackCut(25),
   kLoPlaneEventCut(-1),
   kHiTrackLikeCut(-1),
   kDPlaneCut(-1),
   kLoPhNStripCut(-1),
   kLoPhNPlaneCut(-1),
   kHiEnergyCut(-1),
   kLoEnergyCut(-1),
   kHiEnergyShowerCut(-1),
   kLoEnergyShowerCut(-1),
   kPhStripCut(-1),
   kPhPlaneCut(-1),
   kLoCurrentCut(0.1),
   kLoHorBeamWidth(0.0),
   kHiHorBeamWidth(2.9),
   kLoVertBeamWidth(0.0),
   kHiVertBeamWidth(2.9),
   kLoNuTarZ(-1),
   kHiNuTarZ(1000),
   kOscillate(0),   
   kOutputFile("HistManOut.root")
{}

//......................................................................

CompareAll::~CompareAll()
{}

//......................................................................
void CompareAll::BeginJob()
{

  if(counter%1000==0){
    cout<<"On entry "<<counter<<endl;
  }
  counter++;

  string dumName;
  TString  dumtype;
  Float_t dumstart, dumend;
  Int_t dumbins, nvar;
  ifstream ins;
  ins.open("AllParam.txt");

  nvar = 0;

  fPOT[0] = fPOT[1] = fPOT[2] = fPOT[3] = 0.0;
  for(int i = 0 ; i < 4; i++)
    fOscParams[i] = ANtpDefVal::kFloat;

  //read in the file
  while(!ins.eof()) {
      ins>>dumName>>dumstart>>dumend>>dumbins>>dumtype;
      if(!ins.eof()){
          varName.push_back(dumName);
          beg.push_back(dumstart);  end.push_back(dumend);
          nbins.push_back(dumbins);   gtype.push_back(dumtype);
          nvar++;
         }
    }
  cout<<nvar<<" variables read in"<<endl;

// static HistMan *hm = new HistMan("allcomp");
	
  //create all possible histograms
 
  vector<TString> names;
  names.push_back("_far_mc_nue");
  names.push_back("_far_mc_numu");
  names.push_back("_far_mc_bnue");
  names.push_back("_far_mc_nutau");
  names.push_back("_far_mc_nc");
  names.push_back("_near_mc_bnue");
  names.push_back("_near_mc_numu");
  names.push_back("_near_mc_nc");
  names.push_back("_far_data");
  names.push_back("_near_data");
  names.push_back("_unknown");
 
  
  for(UInt_t i = 0; i < names.size(); i++){
     for(UInt_t l = 0; l < varName.size(); l++){
	   string temp = (varName[l]);
           string dirstring = "allcomp/" + temp.substr(0,temp.find_first_of("."));
	   HistMan hm2(dirstring.c_str());
	   TString param = varName[l] + names[i];
           hm2.Book<TH1F>(param,param,nbins[l],beg[l],end[l]);	   
    }
  }
  
  fHistRecord = new TTree("histRecord", "Historgram Filling Information");

  fHistRecord->Branch("POT", &fPOT, "farMCPOT/F:farDataPOT/F:nearMCPOT/F:nearDataPOT/F");
  fHistRecord->Branch("OscParams", &fOscParams, "Baseline/F:deltam23/F:theta23/F:Ue3/F");

}


JobCResult CompareAll::Ana(const MomNavigator* mom)
{
   //get all NueRecords from mom 
   //may have more than one per go since mom reads in a snarl's worth of data
   //so, this is a little more complicated than just asking for a NueRecord
   TObject *obj=0;
//  static Float_t total_pot = 0;
  static Int_t runno = -1;
  static Int_t snarlno = -1;

//  if its a brand new far det mc file add a far mc number of POT
//  if its a brand new near det mc file add near mc number of POT
// if mc also chekc htat the nueosc params are the same

//  if its data -> add up the number of POT recorded in the bmon nutple for each snarl


//  if(

//  static Int_t printpot = 1;


   TIter objiter = mom->FragmentIter();
   while((obj=objiter.Next())){
      NueRecord *nr = static_cast<NueRecord *>(obj);
      if(nr){
	 MSG("CompareAll",Msg::kDebug)<<"Found a NueRecord in MOM"<<endl;
      }
      else{
	 MSG("CompareAll",Msg::kError)<<"Didn't find a NueRecord in MOM"<<endl;
	 continue;
      }
      MSG("CompareAll",Msg::kDebug)<<"Found a NueRecord "<<nr<<endl;
      SimFlag::SimFlag_t s = nr->GetHeader().GetVldContext().GetSimFlag();
      Detector::Detector_t d = 
	      nr->GetHeader().GetVldContext().GetDetector();
  
//  need at least 1 sucessfully reconstructed event
       
     if(counter%1000==0){
         cout<<"On entry "<<counter<<endl;
     }
     counter++;
      
     if(nr->GetHeader().GetEventNo()<0){
	continue;
     }

     //  if its a brand new far det mc file add a far mc number of POT
     //  //  if its a brand new near det mc file add near mc number of POT
     //  // if mc also chekc htat the nueosc params are the same
     //  


     //if it is data
     if(s==SimFlag::kData){
       if(PassesBeamCuts(nr)){
         if (runno != nr->GetHeader().GetRun() ||
                runno == nr->GetHeader().GetRun() &&
                snarlno != nr->GetHeader().GetSnarl())
	 {
             runno = nr->GetHeader().GetRun();
             snarlno = nr->GetHeader().GetSnarl();

             if(d == Detector::kFar)
	       fPOT[far_data] += nr->bmon.bI;
             if(d == Detector::kNear)
               fPOT[near_data] += nr->bmon.bI;			  
        }
     }                                                                         
   }	   
     if(s==SimFlag::kMC){
       if (runno != nr->GetHeader().GetRun() ) 
       {
           runno = nr->GetHeader().GetRun();
                                                                                
           if(d == Detector::kFar)
               fPOT[far_mc] += 6.5e8;
           if(d == Detector::kNear)
               fPOT[near_mc] += 550*24.;
      }
      if(fOscParams[0] < 0)
      {
        fOscParams[0] = nr->mctrue.Baseline;
	fOscParams[1] = nr->mctrue.DeltamSquared23;
        fOscParams[2] = nr->mctrue.Theta23;
        fOscParams[3] = nr->mctrue.Ue3Squared;		
      }
       
    }
           
     if(PassesCuts(nr)){
       if (s==SimFlag::kData){
         if(PassesBeamCuts(nr)){
           TString id = MakeIdString(nr);
           FillFromList(nr,id,1.);
         }
       }else{
           TString id = MakeIdString(nr);
           Float_t weight  = 1.0;
           if(kOscillate) weight = nr->mctrue.fOscProb*nr->mctrue.fNueWeight;
           FillFromList(nr,id,weight);
        }
     }
   }

   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
void CompareAll::EndJob()
{
//Now i have to output the tree
  // Here is where all of the writeout work will be done, first the tree and then the histos

  //If we are oscillating the files, then the effective exposure is only one third per file
  //    if you are not using an equal number of files from each type I accept no responsibility
  //    for the nature of your results
  if(kOscillate){
    fPOT[far_mc] =  fPOT[far_mc]/3.0;
  }

  TFile* file = new TFile(kOutputFile.c_str(), "update");
  file->cd();
  fHistRecord->Fill();
  fHistRecord->Write();
//  delete file;

  vector<TString> branches;  
  for(UInt_t i = 0; i < varName.size(); i++)
  {
     string temp = (varName[i]);
     string dirstring = "allcomp/" + temp.substr(0,temp.find_first_of("."));
     bool newBranch = true;
     
     for(UInt_t j = 0; j < branches.size() && newBranch; j++)
     {
	 if(dirstring == branches[j])
		newBranch = false;
     }
     if(newBranch) branches.push_back(dirstring);
  }

  for(UInt_t i = 0; i < branches.size(); i++)
  {  
    HistMan *hm2 = new HistMan(branches[i]);
    hm2->WriteOut(*file);
  }


  delete file;  
}

const Registry& CompareAll::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("CompareAll",Msg::kDebug)<<"In CompareAll::DefaultConfig"<<endl;
                                                                                
  static Registry r; // Default configuration for module
                                                                                
  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());
                                                                                
  // Set values in configuration
  r.UnLockValues();
  r.Set("HiPlaneTrackCut",25);
  r.Set("LoPlaneEventCut",-1);
  r.Set("HiTrackLikeCut",-1);
  r.Set("DPlaneCut",-1);
  r.Set("LoPhNStripCut",-1);
  r.Set("LoPhNPlaneCut",-1);
  r.Set("HiEnergyCut",-1);
  r.Set("LoEnergyCut",-1);
  r.Set("HiEnergyShowerCut",-1);
  r.Set("LoEnergyShowerCut",-1);
  r.Set("PhStripCut",-1);
  r.Set("PhPlaneCut",-1);
  r.Set("LoCurrentCut", 0.1);
  r.Set("LoHorBeamWidth", 0.0);
  r.Set("HiHorBeamWidth", 2.9);
  r.Set("LoVertBeamWidth", 0.0);
  r.Set("HiVertBeamWidth", 2.9);
  r.Set("LoNuTarZ", -1);
  r.Set("HiNuTarZ", 1000);
  r.Set("Oscillate", 0);
  r.Set("OutputFile", "HistManInfo.root");
  
  r.LockValues();
                                                                                
  return r;
}

void CompareAll::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("CompareAll",Msg::kDebug)<<"In CompareAll::Config"<<endl;
                                                                                
  //  const char* tmps;
  int imps;
  if(r.Get("HiPlaneTrackCut",imps)) { kHiPlaneTrackCut=imps;}
  if(r.Get("LoPlaneEventCut",imps)) { kLoPlaneEventCut=imps;}
  if(r.Get("HiTrackLikeCut",imps)) { kHiTrackLikeCut=imps;}
  if(r.Get("DPlaneCut",imps)) { kDPlaneCut=imps;}
  if(r.Get("LoPhNStripCut",imps)) { kLoPhNStripCut=imps;}
  if(r.Get("LoPhNPlaneCut",imps)) { kLoPhNPlaneCut=imps;}

  if(r.Get("LoNuTarZ", imps)) {kLoNuTarZ = imps;}
  if(r.Get("HiNuTarZ", imps)) {kHiNuTarZ = imps;}
  if(r.Get("Oscillate", imps)) {kOscillate = imps;}     
  
  double fmps;
  if(r.Get("HiEnergyCut",fmps)) { kHiEnergyCut=fmps;}
  if(r.Get("LoEnergyCut",fmps)) { kLoEnergyCut=fmps;}
  if(r.Get("HiEnergyShowerCut",fmps)) { kHiEnergyShowerCut=fmps;}
  if(r.Get("LoEnergyShowerCut",fmps)) { kLoEnergyShowerCut=fmps;}
  if(r.Get("PhStripCut",fmps)) { kPhStripCut=fmps;}
  if(r.Get("PhPlaneCut",fmps)) { kPhPlaneCut=fmps;}

  if(r.Get("LoCurrentCut", fmps)) {kLoCurrentCut = fmps;}
  if(r.Get("LoHorBeamWidth", fmps)) {kLoHorBeamWidth = fmps;}
  if(r.Get("HiHorBeamWidth", fmps)) {kHiHorBeamWidth = fmps;}
  if(r.Get("LoVertBeamWidth", fmps)) {kLoVertBeamWidth = fmps;}
  if(r.Get("HiVertBeamWidth", fmps)) {kHiVertBeamWidth = fmps;}

  const char* tmps;
  if(r.Get("OutputFile", tmps)) {kOutputFile = tmps;}
  
}




bool CompareAll::PassesBeamCuts(NueRecord* nr)
{
  //bool passes = true;
  if (nr->GetHeader().GetVldContext().GetSimFlag()!=SimFlag::kData) return true;

  if(kLoCurrentCut > 0)
    if (nr->bmon.bI < kLoCurrentCut) return false;      //beam intensity (in 1e12)
  if(kLoHorBeamWidth > 0)    
    if (nr->bmon.hbw*1e3 < kLoHorBeamWidth) return false;     //horizontal beam width
  if(kHiHorBeamWidth > 0)    
    if (nr->bmon.hbw*1e3 > kHiHorBeamWidth) return false;
  if(kLoVertBeamWidth > 0)	  
    if (nr->bmon.vbw*1e3 < kLoVertBeamWidth) return false;     //vertical beam width
  if(kHiVertBeamWidth > 0)
    if (nr->bmon.vbw*1e3 > kHiVertBeamWidth) return false;
  if(kLoNuTarZ > 0)
    if (nr->bmon.nuTarZ < kLoNuTarZ) return false;  //low energy position
  if(kHiNuTarZ > 0)
    if (nr->bmon.nuTarZ > kHiNuTarZ) return false;
    
  return true;
}



bool CompareAll::PassesCuts(NueRecord* nr)
{
    bool passes = true;   

    //dont even bother if there is a long track
    if(nr->srtrack.planes > 25) passes = false;
    if(nr->anainfo.inFiducialVolume != 1)
        passes = false;
    //dont even bother if there is a long track
    if (nr->srtrack.planes > 25) passes = false;
    if (nr->srevent.pulseHeight<2e4) passes = false;
                                                                               
    if (nr->anainfo.isFullyContained != 1 && nr->anainfo.isFullyContained!= -2)
	    passes = false;
                                                                                
    if (nr->srevent.phMeu>150) passes = false;
    if (nr->srevent.hotch==1) passes = false;
    if((TMath::Max(nr->srtrack.pulseHeight,nr->srshower.pulseHeight)<5000))
        passes = false;

   //only look at events for which mst actually gets calculated
   if(nr->GetHeader().GetTrackLength()>25)
       passes = false;
     
   if(!EventFilter::PassesAllCuts(nr,kHiPlaneTrackCut,kHiTrackLikeCut,
			   kLoPlaneEventCut, kHiEnergyCut,
			   kLoEnergyCut, kHiEnergyShowerCut,
			   kLoEnergyShowerCut))
     passes = false;

   return passes;
}


TString CompareAll::MakeIdString(NueRecord *nr)
{
  Detector::Detector_t d = nr->GetHeader().GetVldContext().GetDetector();
  SimFlag::SimFlag_t s = nr->GetHeader().GetVldContext().GetSimFlag();
  
  TString det, dm;
  TString type;

  if(d==Detector::kFar){
    det = "_far";
  }
  else if(d==Detector::kNear){
    det = "_near";
  }
  else{
    return "_unknown";
  }
                                                                        
  if(s==SimFlag::kData){
    dm = "_data";
    TString id = det+dm;
    return id;
  }
  else 
     if (s==SimFlag::kMC){
          dm = "_mc";
                                                                                
       if(nr->mctrue.interactionType==1){
         if(abs(nr->mctrue.nuFlavor)==12 &&
            abs(nr->mctrue.nonOscNuFlavor)==12){
              type = "bnue";
         }
         else if(abs(nr->mctrue.nuFlavor)==12){
               type = "nue";
         }
         else if(abs(nr->mctrue.nuFlavor)==14){
	        type = "numu";
         }
         else if(abs(nr->mctrue.nuFlavor)==16){
	        type = "nutau";
         }
         else{
	        return "unknown";
         }
       }
       else{
	      type = "nc";
       }

       TString id = det+dm+ "_" + type;
       return id;
  }
  else {
	 return "unknown";
       }
  
}


void CompareAll::FillFromList(NueRecord* nr, TString id, Float_t weight)
{
    if(varName.size() == 0) return;
    TString hname;
    UInt_t count = 0;
    
    TClass *cl;
    TRealData *rd;
    string vName;
    TDataMember *member;
    TDataType *membertype;
    Float_t value = 0.0;
    
    cl=nr->IsA();
    TIter  next(cl->GetListOfRealData());                                                                                 
    while ((rd =dynamic_cast<TRealData*>(next()))) {
      member = rd->GetDataMember();
      membertype = member->GetDataType();
      vName=rd->GetName();
                                                                         
      Int_t offset = rd->GetThisOffset();
      char *pointer = (char*)nr  + offset;

      for(UInt_t i = 0; i < varName.size();i++){
        if(vName == varName[i]){
            value = -9999;
            if(!NeedsSpecialAttention(vName, nr, value))
                value=atof(membertype->AsString(pointer));
             MSG("CompareAll",Msg::kDebug)<<"Found variable "
		                        <<vName<<" with value "<<value;
             MSG("CompareAll",Msg::kDebug)<<"Storing it w/ id "<<id<<endl;
	     
 	    
 	    if(!ANtpDefVal::IsDefault(value) && 
	        !ANtpDefVal::IsDefault(static_cast<Double_t> (value)) && 
		!ANtpDefVal::IsDefault(static_cast<Int_t> (value))){

                string dirstring = "allcomp/" + vName.substr(0,vName.find_first_of("."));
                HistMan hm2(dirstring.c_str());			 
		    
   	        hname = varName[i]+id;
	        TH1F* hist = hm2.Get<TH1F>(hname);	
		hist->Fill(value, weight);
             }
             MSG("CompareAll",Msg::kDebug)<<"Found variable "
		   <<vName<<" with value "<<value;

	  count++;
	  i = varName.size();
	  }
      }
      if(count == varName.size()) break;
    }
          
   return;     
}     
           

bool CompareAll::NeedsSpecialAttention(TString name, NueRecord *nr, Float_t &value)
{
   
   //All the fHeaders and four of hte MST vars require special effort
     if(name == "fHeader.fSnarl") {
         value = nr->GetHeader().GetSnarl();
     }if(name == "fHeader.fRun") {
         value = nr->GetHeader().GetRun();
     }if(name == "fHeader.fSubRun") {
         value = nr->GetHeader().GetSubRun();
     }if(name == "fHeader.fEvtNo") {
         value = nr->GetHeader().GetEventNo();
     }if(name == "fHeader.fEvents") {
         value = nr->GetHeader().GetEvents();
     }if(name == "fHeader.fTrackLength") {
         value = nr->GetHeader().GetTrackLength();
     }

     if(name == "mstvars.eallw1") {
        if(nr->mstvars.enn1 > 0) value = 0.0;
        for(int i=0;i<nr->mstvars.enn1;i++){
          value += nr->mstvars.eallw1[i];
        }
     }
     if(name == "mstvars.oallw1") {
        if(nr->mstvars.onn1 > 0) value = 0.0;
        for(int i=0;i<nr->mstvars.onn1;i++){
          value += nr->mstvars.oallw1[i];
        }
     } 
     if(name == "mstvars.eallm1") {
        if(nr->mstvars.enn1 > 0) value = 0.0;
        for(int i=0;i<nr->mstvars.enn1;i++){
          value += nr->mstvars.eallm1[i];
        }
     }
     if(name == "mstvars.oallm1") {
         if(nr->mstvars.onn1 > 0) value = 0.0;
        for(int i=0;i<nr->mstvars.onn1;i++){
          value += nr->mstvars.oallm1[i];
        }
     }    

     if(value > -9999) return true;
     return false;
}
