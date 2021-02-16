////////////////////////////////////////////////////////////////////////
// $Id: NueSensitivity.cxx,v 1.1 2006/09/18 21:41:22 boehm Exp $
//
// FILL_IN: [Document your code!!]
//
// vahle@hep.ucl.ac.uk
////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TFile.h"
#include "TClonesArray.h"
#include "TF1.h"
#include "TH2F.h"
#include "TGraph.h"
#include "NueAna/NueSensitivity.h"
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
#include "NueAna/NueRWHelpers.h"
#include "NueAna/OscProb.h"

using namespace NueRWHelpers;

#define POTPERFARFILE 6.5e20
#define POTPERNEARFILE 550*2.4e13

JOBMODULE(NueSensitivity, "NueSensitivity",
          "Reads in ana_nue ntuple and pid ntuple");
CVSID("$Id:");

//......................................................................
NueSensitivity::NueSensitivity()
{

  nueAppear = NULL;
  numuSurvive = NULL;

  systematicHist_oscSys = NULL;
  systematicHist_allSys = NULL;
  systematicFile = NULL;
  systematicHistNorm = 7.4e20; //Trish's normalisation

  MDCChallengePOT = 7.4e20;
  MDCNearToFar = 1.1125e-03; //ratio of #near to #far based on MDC

  nNuMuFiles = 0;
  nNueFiles = 0;
  nNuTauFiles = 0;
  nNearFiles = 0;
  nChallengeNearFiles = 0;

  nNearUnknownEvents = 0;
  nNearNueEvents = 0;
  nNearNuMuEvents = 0;
  nNearNuTauEvents = 0;
  nNearBeamNueEvents = 0;
  nNearNCEvents = 0;

  nFarUnknownEvents = 0;
  nFarNueEvents = 0;
  nFarNuMuEvents = 0;
  nFarNuTauEvents = 0;
  nFarBeamNueEvents = 0;
  nFarNCEvents = 0;

  nDeltaPoints = 60;
  nThetaPoints = 40;
  theta13 = NULL;
  delta23 = NULL;

  currentRun = 0;

}

//......................................................................
NueSensitivity::~NueSensitivity()
{
  delete [] theta13;
  delete [] delta23;
  delete nueAppear;
  delete numuSurvive;
}

//......................................................................
const Registry& NueSensitivity::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
  static Registry r; // Default configuration for module
 
  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());
 
  //Set values in configuration
  r.UnLockValues();
  r.Set("nNuMuFiles",0);
  r.Set("nNueFiles",0);
  r.Set("nNuTauFiles",0);
  r.Set("nNearFiles",0);
  r.Set("nChallengeNearFiles",0);
  r.LockValues();
 
  return r;
}

//......................................................................
void NueSensitivity::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  
  int tmpi;
  if(r.Get("nNuMuFiles",tmpi))           { nNuMuFiles  = tmpi;}
  if(r.Get("nNueFiles",tmpi))            { nNueFiles   = tmpi;}
  if(r.Get("nNuTauFiles",tmpi))          { nNuTauFiles = tmpi;}
  if(r.Get("nNearFiles",tmpi))           { nNearFiles  = tmpi;}
  if(r.Get("nChallengeNearFiles",tmpi))  { nChallengeNearFiles  = tmpi;}
}


//......................................................................
void NueSensitivity::SetPOT(){

  NuMuFilesPOT          = float(nNuMuFiles)*POTPERFARFILE;
  NueFilesPOT           = float(nNueFiles)*POTPERFARFILE;
  NuTauFilesPOT         = float(nNuTauFiles)*POTPERFARFILE;
  NearFilesPOT          = float(nNearFiles)*POTPERNEARFILE;  
  ChallengeNearFilesPOT = float(nChallengeNearFiles)*POTPERNEARFILE;

  MSG("NueSensitivity",Msg::kInfo) << "POT: " << NuMuFilesPOT << " " 
				   << NueFilesPOT << " "
				   << NuTauFilesPOT << " " 
				   << NearFilesPOT << " " 
				   << ChallengeNearFilesPOT << endl;

}

//......................................................................
void NueSensitivity::BeginJob(){

  Double_t baseline = 735.0;
  nueAppear = new TF1("nueAppear",ElecAppear,0.05,100,10);
  nueAppear->SetParameter(0,baseline); //baseline (km)
  nueAppear->SetParameter(1,0.6); //sinsq_2th23
  nueAppear->SetParameter(2,0.554); //sinsq_2th12
  nueAppear->SetParameter(3,0.0); //sinsq_2th13
  nueAppear->SetParameter(4,0.002); //dmsq23
  nueAppear->SetParameter(5,8.2e-5); //dmsq12
  nueAppear->SetParameter(6,0); //density
  nueAppear->SetParameter(7,0); //cp phase
  nueAppear->SetParameter(8,1); //anti-nu

  numuSurvive = new TF1("numuSurvive",MuSurvive,0.05,100,10);
  numuSurvive->SetParameter(0,baseline); //baseline (km)
  numuSurvive->SetParameter(1,0.6); //sinsq_2th23
  numuSurvive->SetParameter(2,0.554); //sinsq_2th12
  numuSurvive->SetParameter(3,0.0); //sinsq_2th13
  numuSurvive->SetParameter(4,0.002); //dmsq23
  numuSurvive->SetParameter(5,8.2e-5); //dmsq12
  numuSurvive->SetParameter(6,0); //density
  numuSurvive->SetParameter(7,0); //cp phase
  numuSurvive->SetParameter(8,1); //anti-nu

  systematicFile = new TFile("sysHys/moneyplot_newest.root","READ");
  if(systematicFile->IsOpen() && !systematicFile->IsZombie()){
    systematicHist_oscSys = (TH2F*) systematicFile->Get("hfdvnd");
    systematicHist_allSys = (TH2F*) systematicFile->Get("hfdvnd");
  }

  theta13 = new Double_t[nThetaPoints];
  delta23 = new Double_t[nDeltaPoints];

  for(int i=0;i<nThetaPoints;i++){
    theta13[i] = 0.000 + 0.005*Double_t(i);
  }

  for(int i=0;i<nDeltaPoints;i++){
    delta23[i] = 0.001 + 0.0001*Double_t(i);
  }

  HistMan man("sensitivity");
  man.Book<TH2F>("nearHist","nearHist",6,0,6,200,0,50);
  for(int i=0;i<nThetaPoints;i++){
    for(int j=0;j<nDeltaPoints;j++){
      char name[256];
      sprintf(name,"farHist%.0f_%.0f",1000.*theta13[i],10000.*delta23[j]);
      man.Book<TH2F>(name,name,6,0,6,200,0,50);
    }
  }

  man.Book<TH2F>("sensitivityHist1","sensitivityHist1",
		 nThetaPoints,theta13[0]-0.0025,
		 theta13[nThetaPoints-1]+0.0025,
		 nDeltaPoints,delta23[0]-0.00005,
		 delta23[nDeltaPoints-1]+0.00005);
  man.Book<TH2F>("sensitivityHist2","sensitivityHist2",
		 nThetaPoints,theta13[0]-0.0025,
		 theta13[nThetaPoints-1]+0.0025,
		 nDeltaPoints,delta23[0]-0.00005,
		 delta23[nDeltaPoints-1]+0.00005);


  SetPOT();
}

//......................................................................
JobCResult NueSensitivity::Ana(const MomNavigator* mom)
{
  //get all NueRecords from mom 
  //may have more than one per go since mom reads in a snarl's worth of data
  //so, this is a little more complicated than just asking for a NueRecord
  MSG("NueSensitivity",Msg::kVerbose) << "**********IN ANA**********" << endl;

  TObject *obj=0;
  TIter objiter = mom->FragmentIter();
  NtpSRRecord *sr=0;
  NtpMCRecord *mc=0;
  NtpTHRecord *th=0;
  vector<NuePID *> vpid;
  while((obj=objiter.Next())){
    const char *cn=obj->ClassName();
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
      vpid.push_back(npid);
    }
    else{
      continue;
    }
  }

  bool isMC = true;
  bool isTH = true;
  if(!sr){
    MSG("NueSensitivity",Msg::kError)<<"No NtpSRRecord in Mom"<<endl;
    return JobCResult::kFailed;
  }  
  if(!mc) isMC = false;
  if(!th) isTH = false;
  if(vpid.size()==0){
    MSG("NueSensitivity",Msg::kError)<<"No NuePID Records in Mom"<<endl;
    return JobCResult::kFailed;
  }

  if(sr->GetHeader().GetRun()!=currentRun) {
    currentRun = sr->GetHeader().GetRun();
    MSG("NueSensitivity",Msg::kInfo)<< "Current Run: " 
				    << currentRun << endl;
  }

  //use HistMan to plot something for events that pass/fail
  HistMan man("sensitivity");
  Int_t det = 0;
  if(sr->GetHeader().GetVldContext().GetDetector()==Detector::kNear) 
    {
      nueAppear->SetParameter(0,1.0);
      numuSurvive->SetParameter(0,1.0);
      det = 1;
    }
  else if(sr->GetHeader().GetVldContext().GetDetector()==Detector::kFar) 
    {
      nueAppear->SetParameter(0,735.);
      numuSurvive->SetParameter(0,735.);
      det = 2;
    }

  //so, mom will match up snarls for us,
  //but, we have to match up events for ourselves.  
  Int_t *nmc = NULL; 
  if(isMC){ 
    TClonesArray& mcArray = *(mc->mc);
    nmc = new Int_t[mcArray.GetEntries()];
    for(int i=0;i<mcArray.GetEntries();i++) nmc[i] = 0;
  }
  for(unsigned int i=0;i<vpid.size();i++){
    int evtno = vpid[i]->GetHeader().GetEventNo();
    int mcindex=0;
    NtpSREvent *evt = NULL;
    NtpMCTruth *mcth = NULL;
    if(evtno<0){
      MSG("NueSensitivity",Msg::kDebug)<< "can't get event "
				       << evtno << endl;
      if(isMC && isTH) {
	mcth = SntpHelpers::GetMCTruth(mcindex,mc);
	if(mcth==0) continue;
      }
    }
    else{
      evt = SntpHelpers::GetEvent(evtno,sr);
      if(isMC && isTH){
	mcindex = SntpHelpers::GetEvent2MCIndex(evtno,th);
	if(mcindex>=0) nmc[mcindex] += 1;
	mcth = SntpHelpers::GetMCTruth(mcindex,mc);
	if(mcth==0){
	  MSG("NueSensitivity",Msg::kError)<< "can't get mctruth for event "
					   << evtno << endl;
	  continue;
	}
      }
    }

    //nuIntType:
    //0=unknown; 1=nue from oscillated numu; 2=numuCC; 
    //3=nutauCC; 4=beam nue; 5=nc      
    Int_t nuIntType = 0;
    if(isMC) {
      if(mcth->iaction==0){
	nuIntType=5;
	//if(nmc[mcindex]==1) {
	if(det==1) nNearNCEvents +=1;
	else if(det==2) nFarNCEvents +=1;
	//}
      }
      else if(abs(mcth->inu)==12){
	if(abs(mcth->inunoosc)==12){
	  nuIntType=4;
	  //if(nmc[mcindex]==1) {
	  if(det==1) nNearBeamNueEvents +=1;
	  else if(det==2) nFarBeamNueEvents +=1;
	  //}
	}
	else if(abs(mcth->inunoosc)==14){
	  nuIntType=1;
	  //if(nmc[mcindex]==1) {
	  if(det==1) nNearNueEvents +=1;
	  else if(det==2) nFarNueEvents +=1;
	  //}
	}
      }
      else if(abs(mcth->inu)==14){
	nuIntType=2;
	//if(nmc[mcindex]==1) {
	if(det==1) nNearNuMuEvents +=1;
	else if(det==2) nFarNuMuEvents +=1;
	//}
      }
      else if(abs(mcth->inu)==16){
	nuIntType=3;
	//if(nmc[mcindex]==1) {
	if(det==1) nNearNuTauEvents +=1;
	else if(det==2) nFarNuTauEvents +=1;
	//}
      }
    }

    if(nuIntType==0) {
      if(!isMC) {
	if(det==1) nNearUnknownEvents +=1;
	else if(det==2) nFarUnknownEvents +=1;
      }
    }
    
    int pass = vpid[i]->IsNue;
    if(pass==1){
      
      if(evt==0) {
	MSG("NueSensitivity",Msg::kError)<< "Have PID==pass but evt==0!"
					 << endl;
	continue;
      }

      //POT norms:
      Double_t totPOT = 0;
      Double_t NCScale = 1;
      Double_t BeamNueScale = 1;
      Double_t nuMuScale = 1;
      Double_t nuEScale = 1;
      Double_t nuTauScale = 1;
      
      if(det==1){
	//for NearDet scale everything to MDC Challenge Set POT:
	totPOT       = NearFilesPOT;
	NCScale      = MDCChallengePOT/totPOT;
	BeamNueScale = MDCChallengePOT/totPOT;
	nuMuScale    = MDCChallengePOT/totPOT;
	nuEScale     = MDCChallengePOT/totPOT;
	nuTauScale   = MDCChallengePOT/totPOT;
      }
      else if(det==2){
	totPOT = NuMuFilesPOT + NueFilesPOT + NuTauFilesPOT;
	//normalise NC to MDC Challenge POT:
	NCScale = MDCChallengePOT/totPOT;
	//normalise beam nue to MDC Challenge POT:
	//there are no beam nue's in the nutau files!
	BeamNueScale = MDCChallengePOT/(NuMuFilesPOT + NueFilesPOT);
	//normalise numu to MDC Challenge POT: 
	nuMuScale = MDCChallengePOT/NuMuFilesPOT;
	//normalise nue to MDC Challenge POT:
	nuEScale = MDCChallengePOT/NueFilesPOT;
	//normalise nutau to MDC Challenge POT:
	nuTauScale = MDCChallengePOT/NuTauFilesPOT;
      }

      if(det==1){
	if(nuIntType==1){
	  man.Fill2d("nearHist",nuIntType,evt->ph.gev,
		     nuEScale);
	}
	else if(nuIntType==2){
	  man.Fill2d("nearHist",nuIntType,evt->ph.gev,
		     nuMuScale);
	}
	else if(nuIntType==3){
	  man.Fill2d("nearHist",nuIntType,evt->ph.gev,
		     nuTauScale);
	}
	else if(nuIntType==4){
	  man.Fill2d("nearHist",nuIntType,evt->ph.gev,
		     BeamNueScale);
	}
	else if(nuIntType==5){
	  man.Fill2d("nearHist",nuIntType,evt->ph.gev,
		     NCScale);
	}
	else if(nuIntType==0){ //assume that these events are Challenge
	  man.Fill2d("nearHist",nuIntType,evt->ph.gev,
		     MDCChallengePOT/ChallengeNearFilesPOT);  
	}
      }
      else if(det==2){
	for(int j=0;j<nThetaPoints;j++){
	  for(int k=0;k<nDeltaPoints;k++){
	    //histman names:
	    char name[256];
	    sprintf(name,"farHist%.0f_%.0f",1000*theta13[j],10000.*delta23[k]);
	    
	    Double_t nueAppearWeight = 1;
	    Double_t numuSurviveWeight = 1;
	    if(isMC){
	      //osc probs:
	      nueAppear->SetParameter(3,theta13[j]);
	      numuSurvive->SetParameter(3,theta13[j]);
	      nueAppear->SetParameter(4,delta23[k]);
	      numuSurvive->SetParameter(4,delta23[k]);
	      nueAppearWeight = nueAppear->Eval(mcth->p4neu[3]);
	      numuSurviveWeight = numuSurvive->Eval(mcth->p4neu[3]);
	    }
	    double UE32 = 0.5*(1 - sqrt(1-theta13[j]));
	    if(((4.*UE32*(1-UE32)) - theta13[j])>1e-7) {
	      cout << "calc not right: UE32 = " << UE32 
		   << " theta13 = " << theta13[j] << endl;
	    }
	    
// 	    cout << "nue appearance prob: " << endl;
// 	    cout << nueAppearWeight - 
// 	      Oscillate(12,14,mcth->p4neu[3],
// 			735.,delta23[k],0.65,UE32)
// 		 << endl;
// 	    cout << "numu dissappearance prob: " << endl;
// 	    cout << numuSurviveWeight-
// 	      Oscillate(14,14,mcth->p4neu[3],
// 			735.,delta23[k],0.65,UE32)
// 		 << endl;
// 	    cout << "nutau appearance prob: " << endl;
// 	    cout << 1.-numuSurviveWeight -
// 	      Oscillate(16,14,mcth->p4neu[3],
// 			      735.,delta23[k],0.65,UE32)
// 		 << endl;

	    if(nuIntType==1){
	      man.Fill2d(name,nuIntType,evt->ph.gev,
			 nueAppearWeight*nuEScale);
	      //man.Fill2d(name,nuIntType,evt->ph.gev,
	      //	 Oscillate(12,14,mcth->p4neu[3],
	      //		   735.,delta23[k],0.65,
	      //		   UE32)*nuEScale);
	    }
	    else if(nuIntType==2){
	      man.Fill2d(name,nuIntType,evt->ph.gev,
	      	 numuSurviveWeight*nuMuScale);
	      // man.Fill2d(name,nuIntType,evt->ph.gev,
	      // Oscillate(14,14,mcth->p4neu[3],
	      //	   735.,delta23[k],0.65,
	      //	   UE32)*nuMuScale);
	    }
	    else if(nuIntType==3){
	      man.Fill2d(name,nuIntType,evt->ph.gev,
	      	 (1.-numuSurviveWeight)*nuTauScale);
	      //man.Fill2d(name,nuIntType,evt->ph.gev,
	      // Oscillate(16,14,mcth->p4neu[3],
	      //	   735.,delta23[k],0.65,
	      //	   UE32)*nuTauScale);
	    }
	    else if(nuIntType==4){
	      man.Fill2d(name,nuIntType,evt->ph.gev,
	      	 BeamNueScale);
	      //man.Fill2d(name,nuIntType,evt->ph.gev,
	      // Oscillate(12,12,mcth->p4neu[3],
	      //	   735.,delta23[k],0.65,
	      //	   UE32)*BeamNueScale);
	    }
	    else if(nuIntType==5){
	      man.Fill2d(name,nuIntType,evt->ph.gev,
			 NCScale);
	    }
	    else if(nuIntType==0){ //assume that these events are Challenge
	      man.Fill2d(name,nuIntType,evt->ph.gev);
	    }
	  }
	}
      }
    }
  }
  delete [] nmc;
  return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

//......................................................................
void NueSensitivity::Analysis()
{
  
  HistMan man("sensitivity");
  
  MSG("NueSensitivity",Msg::kInfo) << "==============================" << endl;
  MSG("NueSensitivity",Msg::kInfo) 
    << "Job Summary: (considering deltam^{2}_{23} = 0.0025 eV^{2})" << endl;
  MSG("NueSensitivity",Msg::kInfo) << "------------------------------" << endl;
  MSG("NueSensitivity",Msg::kInfo) << "Near Detector: " << endl;
  MSG("NueSensitivity",Msg::kInfo) << " Total Events = " 
				   << man.Get<TH2F>("nearHist")->Integral()
				   << endl;
  MSG("NueSensitivity",Msg::kInfo) << "------------------------------" << endl;
  MSG("NueSensitivity",Msg::kInfo) << "Far Detector: " << endl;
  for(int i=0;i<nThetaPoints;i++){
    char name[256];
    sprintf(name,"farHist%.0f_%.0f",1000.*theta13[i],25.);
    MSG("NueSensitivity",Msg::kInfo) << "SinSq(2Theta13) = " << theta13[i] 
				     << " Total Events = " 
				     << man.Get<TH2F>(name)->Integral()
				     << endl;
  }

  MSG("NueSensitivity",Msg::kInfo) << "------------------------------" << endl;
  MSG("NueSensitivity",Msg::kInfo) << "Event Totals: (Near)  (Far)" << endl;
  MSG("NueSensitivity",Msg::kInfo) << "Unknown: " << nNearUnknownEvents << " "
				  << nFarUnknownEvents << endl;
  MSG("NueSensitivity",Msg::kInfo) << "Nue: " << nNearNueEvents << " "
				  << nFarNueEvents << endl;
  MSG("NueSensitivity",Msg::kInfo) << "NuMu: " << nNearNuMuEvents << " "
				  << nFarNuMuEvents << endl;
  MSG("NueSensitivity",Msg::kInfo) << "NuTau: " << nNearNuTauEvents << " "
				  << nFarNuTauEvents << endl;
  MSG("NueSensitivity",Msg::kInfo) << "BeamNue: " << nNearBeamNueEvents << " "
				  << nFarBeamNueEvents << endl;
  MSG("NueSensitivity",Msg::kInfo) << "NC: " << nNearNCEvents << " "
				  << nFarNCEvents << endl;


  MSG("NueSensitivity",Msg::kInfo) << "------------------------------" 
				  << endl;
  MSG("NueSensitivity",Msg::kInfo) << "POT Summary: (in units of 1e20)" 
				  << endl;
  MSG("NueSensitivity",Msg::kInfo) << "Near Total: " 
				  << NearFilesPOT/1e20 << endl;
  MSG("NueSensitivity",Msg::kInfo) << "Far Total: " 
				  << (NuMuFilesPOT + 
				      NueFilesPOT + 
				      NuTauFilesPOT)/1e20 
				  << endl;
  MSG("NueSensitivity",Msg::kInfo) << "consisting of: " << endl; 
  MSG("NueSensitivity",Msg::kInfo) << "numu POT: " 
				  << NuMuFilesPOT/1e20 << endl;
  MSG("NueSensitivity",Msg::kInfo) << "nue POT: " 
				  << NueFilesPOT/1e20 << endl;
  MSG("NueSensitivity",Msg::kInfo) << "nutau POT: " 
				  << NuTauFilesPOT/1e20 << endl;
  MSG("NueSensitivity",Msg::kInfo) << "==============================" 
				  << endl;


  //Make 2D sensitivity plot:
  TH2F *sensitivityHist1 = man.Get<TH2F>("sensitivityHist1");
  TH2F *sensitivityHist2 = man.Get<TH2F>("sensitivityHist2");

  for(int i=0;i<nDeltaPoints;i++){
    char name[256];
    sprintf(name,"farHist%.0f_%.0f",1000.*theta13[0],10000.*delta23[i]);
    TH2F *bkgFar = man.Get<TH2F>(name);
    TH1D *bkgHist = bkgFar->ProjectionY("bkgHist");
    
    for(int j=0;j<nThetaPoints;j++){
      sprintf(name,"farHist%.0f_%.0f",1000.*theta13[j],10000.*delta23[i]);
      TH2F *theFar = man.Get<TH2F>(name);
      TH1D *sigHist = theFar->ProjectionY("sigHist");

      Float_t chi2 = 0;
      for(int k=1;k<=200;k++){
	double sig = sigHist->GetBinContent(k);
	double bkg = bkgHist->GetBinContent(k);
	if(sig==0||bkg==0) chi2 += 2.*(bkg-sig);
	else chi2 += 2.*(bkg-sig)+2.*sig*TMath::Log(sig/bkg);
      }

      sensitivityHist1->Fill(theta13[j],delta23[i],chi2);
      sensitivityHist2->Fill(theta13[j],delta23[i],
			     2.*(bkgFar->Integral()-theFar->Integral()) +
			     2.*theFar->Integral() *
			     TMath::Log(theFar->Integral() / 
					bkgFar->Integral()));

      delete sigHist;
    }
    delete bkgHist;
  }

  //for calculating prob given difference in #events:
  TF1 *probgaus = new TF1("probgaus","gaus",-10,10);
  probgaus->SetParameters(1./TMath::Sqrt(2.*TMath::Pi()),0,1);
  
  //check that near det is present:
  TH2F *nearHist = man.Get<TH2F>("nearHist");
  Float_t normNear = nearHist->Integral();  
  if(normNear==0) return; //if no near det then don't try to do study
  
  //for estimating error on FD from ND:
  TF1 *gaus = new TF1("gaus","gaus",0,100);

  //calculate simple normalisation for Near to Far:
  Float_t normFar = man.Get<TH2F>("farHist0_25")->Integral();
  Float_t norm = normFar/normNear;

  //check to see if Near detector files are challenge set
  bool isChallenge = false;
  TH1D * normNearHistProj = nearHist->ProjectionX("normNearHistProj");
  if(normNearHistProj->GetBinContent(1) > 0) isChallenge = true;
  delete normNearHistProj;

  //remember some things for different theta_13's:
  Double_t *testStat_noSys  = new Double_t[nThetaPoints];
  Double_t *testProb_noSys  = new Double_t[nThetaPoints];
  Double_t *testStat_oscSys = new Double_t[nThetaPoints];
  Double_t *testProb_oscSys = new Double_t[nThetaPoints];
  Double_t *testStat_allSys = new Double_t[nThetaPoints];
  Double_t *testProb_allSys = new Double_t[nThetaPoints];
  for(int i=0;i<nThetaPoints;i++){
    testStat_noSys[i]  = 0;    testProb_noSys[i]  = 0;
    testStat_oscSys[i] = 0;    testProb_oscSys[i] = 0;
    testStat_allSys[i] = 0;    testProb_allSys[i] = 0;
  }

  //make TGraphs for a particular delta m^{2} (= 0.0025)
  for(int i=0;i<nThetaPoints;i++){
    char name[256];
    sprintf(name,"farHist%.0f_%.0f",1000.*theta13[i],25.);
    
    //predict number of far events from near measurement (using simple norm)
    Float_t nEventsFromNear = normNear*MDCNearToFar;
    
    TH2F *farHist = man.Get<TH2F>(name);
    Float_t nEvents = farHist->Integral();
    //if ND challenge present, use measured event rate to correct far rate
    if(isChallenge) nEvents *= MDCNearToFar/norm;
    Float_t statError = TMath::Sqrt(nEvents);
    
    //One-tail hypothesis test no Sys:
    testStat_noSys[i] = (nEvents-nEventsFromNear)/statError;
    testProb_noSys[i] = probgaus->Integral(testStat_noSys[i],10.);
    
    if(systematicHist_oscSys){
      Int_t bin = 1;
      bin = systematicHist_oscSys->GetXaxis()->FindBin(systematicHistNorm * 
						       nearHist->Integral() / 
						       MDCChallengePOT);
      
      TH1D *tempHist = systematicHist_oscSys->ProjectionY("tempHist",bin,bin);
      tempHist->Fit("gaus");
      Double_t mean=gaus->GetParameter(1)*MDCChallengePOT/systematicHistNorm;;
      Double_t sigma=gaus->GetParameter(2)*MDCChallengePOT/systematicHistNorm;;
      Double_t rms=tempHist->GetRMS()*MDCChallengePOT/systematicHistNorm;;
      //don't trust sigma, if rms is much smaller
      //if(sigma>2.*rms)
      sigma = rms;
      mean = nEventsFromNear;
      
      //One-tail hypothesis test osc Sys:
      testStat_oscSys[i] = (nEvents-mean)/sqrt(statError*statError + 
					       sigma*sigma);
      testProb_oscSys[i] = probgaus->Integral(testStat_oscSys[i],10.);
      delete tempHist;
    }
    
    if(systematicHist_allSys){
      Int_t bin = 1;
      bin = systematicHist_allSys->GetXaxis()->FindBin(systematicHistNorm * 
						       nearHist->Integral() / 
						       MDCChallengePOT);
      
      TH1D *tempHist = systematicHist_allSys->ProjectionY("tempHist",bin,bin);
      tempHist->Fit("gaus");
      Double_t mean=gaus->GetParameter(1)*MDCChallengePOT/systematicHistNorm;
      Double_t sigma=gaus->GetParameter(2)*MDCChallengePOT/systematicHistNorm;
      Double_t rms=tempHist->GetRMS()*MDCChallengePOT/systematicHistNorm;
      //don't trust sigma, if rms is much smaller
      if(sigma>2.*rms) sigma = rms;
      
      //One-tail hypothesis test all Sys:
      testStat_allSys[i] = (nEvents-mean)/sqrt(statError*statError + 
					       sigma*sigma);
      testProb_allSys[i] = probgaus->Integral(testStat_allSys[i],10.);
      delete tempHist;
    }

  }

  TGraph *sensitivityGraph1_noSys = new TGraph(40,theta13,testStat_noSys);
  sensitivityGraph1_noSys->SetName("sensitivityGraph1_noSys");
  sensitivityGraph1_noSys->SetTitle("One-Tailed Hypothesis Test Statistic vs Sin^{2}(2#theta_{13}) - No Systematics");
  man.Adopt("",sensitivityGraph1_noSys);

  TGraph *sensitivityGraph2_noSys = new TGraph(40,theta13,testProb_noSys);
  sensitivityGraph2_noSys->SetName("sensitivityGraph2_noSys");
  sensitivityGraph2_noSys->SetTitle("Probability Data Is Consistent with #theta_{13}=0 vs Sin^{2}(2#theta_{13}) - No Systematics");
  man.Adopt("",sensitivityGraph2_noSys);

  TGraph *sensitivityGraph1_oscSys = new TGraph(40,theta13,testStat_oscSys);
  sensitivityGraph1_oscSys->SetName("sensitivityGraph1_oscSys");
  sensitivityGraph1_oscSys->SetTitle("One-Tailed Hypothesis Test Statistic vs Sin^{2}(2#theta_{13}) - Oscillation Systematics");
  man.Adopt("",sensitivityGraph1_oscSys);

  TGraph *sensitivityGraph2_oscSys = new TGraph(40,theta13,testProb_oscSys);
  sensitivityGraph2_oscSys->SetName("sensitivityGraph2_oscSys");
  sensitivityGraph2_oscSys->SetTitle("Probability Data Is Consistent with #theta_{13}=0 vs Sin^{2}(2#theta_{13}) - Oscillation Systematics");
  man.Adopt("",sensitivityGraph2_oscSys);

  TGraph *sensitivityGraph1_allSys = new TGraph(40,theta13,testStat_allSys);
  sensitivityGraph1_allSys->SetName("sensitivityGraph1_allSys");
  sensitivityGraph1_allSys->SetTitle("One-Tailed Hypothesis Test Statistic vs Sin^{2}(2#theta_{13}) - Oscillation + Neugen Systematics");
  man.Adopt("",sensitivityGraph1_allSys);

  TGraph *sensitivityGraph2_allSys = new TGraph(40,theta13,testProb_allSys);
  sensitivityGraph2_allSys->SetName("sensitivityGraph2_allSys");
  sensitivityGraph2_allSys->SetTitle("Probability Data Is Consistent with #theta_{13}=0 vs Sin^{2}(2#theta_{13}) - Oscillation + Neugen Systematics");
  man.Adopt("",sensitivityGraph2_allSys);

  delete [] testStat_noSys;     delete [] testProb_noSys;
  delete [] testStat_oscSys;    delete [] testProb_oscSys;
  delete [] testStat_allSys;    delete [] testProb_allSys;
  
}

void NueSensitivity::EndJob()
{
  Analysis();
  HistMan man("sensitivity");
  man.WriteOut("SensitivityFile.root");
  systematicFile->Close();
}


////////////////////////////////////////////////////////////////////////
