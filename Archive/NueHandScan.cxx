////////////////////////////////////////////////////////////////////////
// $Id: NueHandScan.cxx,v 1.2 2008/11/19 18:22:51 rhatcher Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#include "TFile.h"
#include "TRandom.h"
#include "TString.h"
#include "NueAna/NueHandScan.h"
#include "NueAna/NuePID.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "JobControl/JobCEnv.h"
#include "MCNtuple/NtpMCRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include "TruthHelperNtuple/NtpTHRecord.h"
#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "NueAna/SntpHelpers.h"
#include <iostream>
#include <fstream>

JOBMODULE(NueHandScan, "NueHandScan",
          "Makes Random Hand Scan file from ntuples");
CVSID("$Id:");

//......................................................................
NueHandScan::NueHandScan() : 
  allEntries(0),randomSeed(0),fracPassed(0.8),
  preScaleFactorNue(1),preScaleFactorNuMu(1),preScaleFactorNuTau(1),
  preScaleFactorBNue(1),preScaleFactorNC(1),
  nNueCC(0),nNueNC(0),nNuMuCC(0),nNuMuNC(0),
  nNuTauCC(0),nNuTauNC(0),nBeamNueCC(0),nBeamNueNC(0),nPassNueCC(0),
  nPassNueNC(0),nPassNuMuCC(0),nPassNuMuNC(0),nPassNuTauCC(0),
  nPassNuTauNC(0),nPassBeamNueCC(0),nPassBeamNueNC(0)
{
  firstPass = true;
}

//......................................................................
NueHandScan::~NueHandScan()
{}

//......................................................................
JobCResult NueHandScan::Ana(const MomNavigator* mom)
{

  std::string curFileName = this->GetCurrentFile();
  if(firstPass) {
    ntupleFileNames.push_back(curFileName);
    MSG("NueHandScan",Msg::kInfo) << curFileName << endl;
    gRandom->SetSeed(randomSeed);
    randomSeed = gRandom->GetSeed();
    firstPass = false;  
  }

  std::vector<std::string>::iterator iter = ntupleFileNames.begin();
  std::vector<std::string>::iterator end = ntupleFileNames.end();
  Bool_t fileFlag = false;
  while(iter!=end){
    if(strcmp(iter->c_str(),curFileName.c_str())==0) fileFlag = true;
    iter++;
  }
  if(!fileFlag) {
    ntupleFileNames.push_back(curFileName);
    MSG("NueHandScan",Msg::kInfo) << curFileName << endl;  
  }

  TObject *obj=0;
  TIter objiter = mom->FragmentIter();
  NtpSRRecord *sr=0;
  NtpStRecord *st=0;
  NtpMCRecord *mc=0;
  NtpTHRecord *th=0;
  vector<NuePID *> vpid;
  while((obj=objiter.Next())){
    const char *cn=obj->ClassName();
    if(strcmp(cn,"NtpSRRecord")==0){
      sr = dynamic_cast<NtpSRRecord *>(obj);
    }
    else if(strcmp(cn,"NtpStRecord")==0){
      st = dynamic_cast<NtpStRecord *>(obj);
      }
    else if(strcmp(cn,"NtpMCRecord")==0){
      mc = dynamic_cast<NtpMCRecord *>(obj);
    }
    else if(strcmp(cn,"NtpTHRecord")==0){
      th = dynamic_cast<NtpTHRecord *>(obj);
    }
    else if(strcmp(cn,"NuePID")==0){
      NuePID *npid  = dynamic_cast<NuePID *>(obj);
      vpid.push_back(npid);
    }
    else{
      continue;
    }
  }
  if(!sr&&!st){
    MSG("NueHandScan",Msg::kError)<<"No NtpSR(t)Record in Mom"<<endl;
    return JobCResult::kFailed;
  }
  if(!st&&!mc){
    MSG("NueHandScan",Msg::kError)<<"No NtpMCRecord in Mom"<<endl;
    return JobCResult::kFailed;
  }
  if(!st&&!th){
    MSG("NueHandScan",Msg::kError)<<"No NtpTHRecord in Mom"<<endl;
    return JobCResult::kFailed;
  }
  //if(vpid.size()==0){
  //MSG("NueHandScan",Msg::kError)<<"No NuePID Records in Mom"<<endl;     
  //}
  
  Int_t nEvent = 0;
  if(st) nEvent = st->evthdr.nevent;
  else if(sr) nEvent = sr->evthdr.nevent;
  
  for(int i=0;i<nEvent;i++){
    
    //in case want to apply cuts based on pid later
    Int_t pidevno = -1;
    if(vpid.size()>0) {
      for(unsigned int j=0;j<vpid.size();j++){
	if(vpid[j]->GetHeader().GetEventNo()==i){
	  pidevno = j;
	}
      }
    }
    
    NtpSREvent *event = 0;     
    Int_t mcindex = -1;
    NtpMCTruth *mcth = 0;
    Int_t detector = 0;
    if(st) {
      if(st->GetHeader().GetVldContext().
	 GetDetector()==Detector::kNear) detector = 1;
      else if(st->GetHeader().GetVldContext().
	      GetDetector()==Detector::kFar) detector = 2;
      event = SntpHelpers::GetEvent(i,st);
      mcindex = SntpHelpers::GetEvent2MCIndex(i,st);
      mcth = SntpHelpers::GetMCTruth(mcindex,st);
      runMap[allEntries] = st->GetHeader().GetRun();
      subrunMap[allEntries] = st->GetHeader().GetSubRun();
      snarlMap[allEntries] = st->GetHeader().GetSnarl();
      eventMap[allEntries] = i;
    }
    else if(sr){ 
      if(sr->GetHeader().GetVldContext().
	 GetDetector()==Detector::kNear) detector = 1;
      else if(sr->GetHeader().GetVldContext().
	      GetDetector()==Detector::kFar) detector = 2;
      event = SntpHelpers::GetEvent(i,sr);
      mcindex = SntpHelpers::GetEvent2MCIndex(i,th);
      mcth = SntpHelpers::GetMCTruth(mcindex,mc);
      runMap[allEntries] = sr->GetHeader().GetRun();
      subrunMap[allEntries] = sr->GetHeader().GetSubRun();
      snarlMap[allEntries] = sr->GetHeader().GetSnarl();
      eventMap[allEntries] = i;
    }

    Double_t preScale = 1;

    if(mcth->iaction == 1){
      if(TMath::Abs(mcth->inu)==14) { nNuMuCC++; preScale = preScaleFactorNuMu; }
      else if(TMath::Abs(mcth->inu)==16) { nNuTauCC++; preScale = preScaleFactorNuTau; }
      else if(TMath::Abs(mcth->inu)==12) {
	if(TMath::Abs(mcth->inunoosc)==12) { nBeamNueCC++; preScale = preScaleFactorBNue; }
	else { nNueCC++; preScale = preScaleFactorNue; }
      }
    }
    else {
      preScale = preScaleFactorNC;
      if(TMath::Abs(mcth->inu)==14) nNuMuNC++;
      else if(TMath::Abs(mcth->inu)==16) nNuTauNC++;
      else if(TMath::Abs(mcth->inu)==12) {
	if(TMath::Abs(mcth->inunoosc)==12) nBeamNueNC++;
	else nNueNC++;
      }
    }
    
    if(gRandom->Uniform()<=preScale && 
       this->PassCuts(event,detector)) {
      passMap[allEntries] = 1;
      if(mcth->iaction == 1){
	if(TMath::Abs(mcth->inu)==14) nPassNuMuCC++;
	else if(TMath::Abs(mcth->inu)==16) nPassNuTauCC++;
	else if(TMath::Abs(mcth->inu)==12) {
	  if(TMath::Abs(mcth->inunoosc)==12) nPassBeamNueCC++;
	  else nPassNueCC++;
	}
      }
      else {
	if(TMath::Abs(mcth->inu)==14) nPassNuMuNC++;
	else if(TMath::Abs(mcth->inu)==16) nPassNuTauNC++;
	else if(TMath::Abs(mcth->inu)==12) {
	  if(TMath::Abs(mcth->inunoosc)==12) nPassBeamNueNC++;
	  else nPassNueNC++;
	}
      }
    }
    else passMap[allEntries] = 0;
    allEntries+=1;
  }
  return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................
Bool_t NueHandScan::PassCuts(NtpSREvent *event,Int_t detector)
{
  if(!event) return false;
  if(event->plane.n>16) return false;
  if(detector==1){
    if(event->vtx.z<1||event->vtx.z>5) return false;
    if(TMath::Sqrt(TMath::Power(event->vtx.x-1.4885,2) +
		   TMath::Power(event->vtx.y-0.1397,2)) > 1) return false;
  }
  else if(detector==2){
    if(event->vtx.z<1||event->vtx.z>29) return false;
    if(event->vtx.z>14&&event->vtx.z<17) return false;
    if(TMath::Sqrt(TMath::Power(event->vtx.x,2) + 
		   TMath::Power(event->vtx.y,2)) > 3.87) return false;
  }
  return true;
}

//......................................................................
void NueHandScan::EndJob()
{

  char filename[256];
  sprintf(filename,"randomEventFile_%s.txt",fileTag.c_str());
  ofstream outFile(filename);
  
  std::vector<Int_t> alreadyGot(allEntries,0);
  for(int i=0;i<allEntries;i++) alreadyGot[i] = 0;

  Int_t counter = 0;
  Int_t targetNumber = Int_t(fracPassed*(nPassNueCC + nPassNueNC + 
					 nPassNuMuCC + nPassNuMuNC + 
					 nPassNuTauCC + nPassNuTauNC + 
					 nPassBeamNueCC + nPassBeamNueNC));
  while(counter<targetNumber){
    Int_t ev_no = Int_t(Double_t(allEntries)*gRandom->Uniform());
    if(passMap[ev_no]==1&&alreadyGot[ev_no]==0) {
      outFile << runMap[ev_no] << " "
	      << subrunMap[ev_no] << " "
	      << snarlMap[ev_no] << " "
	      << eventMap[ev_no] << " 0 0" << endl;
      alreadyGot[ev_no] = 1;
      counter++;
    }
  }
  outFile.close();

  sprintf(filename,"EfficiencySummary_%s.txt",fileTag.c_str());
  ofstream outSum(filename);
  outSum << "The following ntuple files are required" << endl;
  outSum << "to use the associated random event file:" << endl;
  std::vector<std::string>::iterator iter = ntupleFileNames.begin();
  std::vector<std::string>::iterator end = ntupleFileNames.end();
  while(iter!=end) {
    TString ts(iter->c_str());
    char c[2] = "/";
    ts.Remove(0,ts.Last(c[0])+1);
    outSum << ts.Data() << endl;
    iter++;
  }
  outSum << "------------------------------------------------------------------" << endl;
  outSum << "Total number of reconstructed events in files = " << allEntries << endl; 
  outSum << "Total number of events passing cuts = " 
	 << (nPassNueCC + nPassNueNC + nPassNuMuCC + nPassNuMuNC + 
	     nPassNuTauCC + nPassNuTauNC + nPassBeamNueCC + nPassBeamNueNC) << endl; 
  outSum << "Fraction of events passing cuts in random file = " << fracPassed << endl;
  outSum << "RandomSeed = " << randomSeed << endl;
  outSum << "Pre-Scale Factors (Nue,NuMu,NuTau,BNue,NC): " 
	 << preScaleFactorNue << ", " << preScaleFactorNuMu << ", "
	 << preScaleFactorNuTau << ", " << preScaleFactorBNue << ", "
	 << preScaleFactorNC << endl;
  outSum << "(pre-scale factors are applied in addition to cuts in table below)" << endl;
  outSum << "------------------------------------------------------------------" << endl;  
  outSum << "          \t" << "Total\tPassed\tEff.(%)\tError(%)" << endl;
  outSum << "NueCC:    \t" << nNueCC << "\t" << nPassNueCC << "\t";
  if(nNueCC>0) outSum << 100.*nPassNueCC/nNueCC
		   << "\t" << 100.*sqrt(nPassNueCC)/nNueCC << endl;
  else outSum << "0\t0" << endl;
  outSum << "NuMuCC:   \t" << nNuMuCC << "\t" << nPassNuMuCC << "\t";
  if(nNuMuCC>0) outSum << 100.*nPassNuMuCC/nNuMuCC 
		  << "\t" << 100.*sqrt(nPassNuMuCC)/nNuMuCC << endl;
  else outSum << "0\t0" << endl;
  outSum << "NuTauCC:  \t" << nNuTauCC << "\t" << nPassNuTauCC << "\t";
  if(nNuTauCC>0) outSum << 100.*nPassNuTauCC/nNuTauCC 
		     << "\t" << 100.*sqrt(nPassNuTauCC)/nNuTauCC << endl;
  else outSum << "0\t0" << endl;
  outSum << "BeamNueCC:\t" << nBeamNueCC << "\t" << nPassBeamNueCC << "\t";
  if(nBeamNueCC>0) outSum << 100.*nPassBeamNueCC/nBeamNueCC
		       <<"\t"<<100.*sqrt(nPassBeamNueCC)/nBeamNueCC << endl;  
  else outSum << "0\t0" << endl;  
  outSum << "NueNC:    \t" << nNueNC << "\t" << nPassNueNC << "\t";
  if(nNueNC>0) outSum << 100.*nPassNueNC/nNueNC 
		   << "\t" << 100.*sqrt(nPassNueNC)/nNueNC << endl;  
  else outSum << "0\t0" << endl;
  outSum << "NuMuNC:   \t" << nNuMuNC << "\t" << nPassNuMuNC << "\t";
  if(nNuMuNC>0) outSum << 100.*nPassNuMuNC/nNuMuNC
		    << "\t" << 100.*sqrt(nPassNuMuNC)/nNuMuNC << endl;  
  else outSum << "0\t0" << endl;
  outSum << "NuTauNC:  \t" << nNuTauNC << "\t" << nPassNuTauNC << "\t";
  if(nNuTauNC>0) outSum << 100.*nPassNuTauNC/nNuTauNC 
		     << "\t" << 100.*sqrt(nPassNuTauNC)/nNuTauNC << endl;  
  else outSum << "0\t0" << endl;
  outSum << "BeamNueNC:\t" << nBeamNueNC << "\t" << nPassBeamNueNC << "\t";
  if(nBeamNueNC>0) outSum << 100.*nPassBeamNueNC/nBeamNueNC
		       << "\t" << 100.*sqrt(nPassBeamNueNC)/nBeamNueNC<<endl;  
  else outSum << "0\t0" << endl;
  outSum << "------------------------------------------------------------------" << endl;
  outSum.close();
}

//......................................................................
const Registry& NueHandScan::DefaultConfig() const
{
  static Registry r; // Default configuration for module 
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());  
  r.UnLockValues();
  r.Set("RandomSeed",0);
  r.Set("FracPassed",0.8);
  r.Set("FileTag","test");
  r.Set("PreScaleFactorNue",1);
  r.Set("PreScaleFactorNuMu",1);
  r.Set("PreScaleFactorNuTau",1);
  r.Set("PreScaleFactorBNue",1);
  r.Set("PreScaleFactorNC",1);
  r.LockValues();
  return r;
}

//......................................................................
void NueHandScan::Config(const Registry& r)
{
  Int_t tmpi = 0;
  Double_t tmpd = 0;
  const char *tmpc = 0;
  if(r.Get("RandomSeed",tmpi))          { randomSeed          = tmpi; }
  if(r.Get("FracPassed",tmpd))          { fracPassed          = tmpd; }
  if(r.Get("FileTag",tmpc))             { fileTag             = tmpc; }
  if(r.Get("PreScaleFactorNue",tmpd))   { preScaleFactorNue   = tmpd; }
  if(r.Get("PreScaleFactorNuMu",tmpd))  { preScaleFactorNuMu  = tmpd; }
  if(r.Get("PreScaleFactorNuTau",tmpd)) { preScaleFactorNuTau = tmpd; }
  if(r.Get("PreScaleFactorBNue",tmpd))  { preScaleFactorBNue  = tmpd; }
  if(r.Get("PreScaleFactorNC",tmpd))    { preScaleFactorNC    = tmpd; }
}

////////////////////////////////////////////////////////////////////////
