////////////////////////////////////////////////////////////////////////
//
// FILL_IN: [Document your code!!]
//
// boehm@physics.harvard.edu
////////////////////////////////////////////////////////////////////////
#include "TH1F.h"
#include "TFile.h"
#include "Conventions/Detector.h"
#include "Conventions/SimFlag.h"
#include "NueAna/CompareMD.h"
#include "NueAna/NueRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "HistMan/HistMan.h"
#include <fstream>
#include "AnalysisNtuples/ANtpDefaultValue.h"

const int CALend=110;
const int SM1end=230;
const int SM2end=484;

JOBMODULE(CompareMD, "CompareMD",
          "Makes hists to All variables between detectors and truth classes");
CVSID("$Id:");
//......................................................................

CompareMD::CompareMD():
  counter(0),
  kNMCFiles(1)
{}

//......................................................................

CompareMD::~CompareMD()
{}

//......................................................................
void CompareMD::BeginJob()
{

  if(counter%1000==0){
    cout<<"On entry "<<counter<<endl;
  }
  counter++;

  TString dumName, dumtype;
  Float_t dumstart, dumend;
  Int_t dumbins, nvar;
  ifstream ins;
  ins.open("AllParam.txt");

  nvar = 0;

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


 
 static HistMan *hm = new HistMan("allcomp");

  for(int d=0;d<2;d++){
    char dm[3];
    if(d==0){
       sprintf(dm,"_d");
    }
    else if(d==1){
      sprintf(dm,"_m");
    }
    else{
      sprintf(dm,"_u");
    }
    if (d==0){//data
      char id[20];
      sprintf(id,"%s",dm);
      for(UInt_t l = 0; l < varName.size(); l++){
	TString param = varName[l] + id;
	hm->Book<TH1F>(param,param,nbins[l],beg[l],end[l]);
      }
    }
    else {//mc
      for(int b=0;b<5;b++){
	char bgc[6];
	if(b==0){
	  sprintf(bgc,"all");
	}
	else if(b==1){
	  sprintf(bgc,"numu");
	}
	else if(b==2){
	  sprintf(bgc,"nutau");
	}
	else if(b==3){
	  sprintf(bgc,"nc");
	}
	else if(b==4){
	  sprintf(bgc,"bnue");
	} 
	else{
	  sprintf(bgc,"u");
	}
     
	char id[20];
	sprintf(id,"%s_%s",dm,bgc);

        for(UInt_t l = 0; l < varName.size(); l++){
           TString param = varName[l] + id;
           hm->Book<TH1F>(param,param,nbins[l],beg[l],end[l]);
        }
      }
    }
  }
}
JobCResult CompareMD::Ana(const MomNavigator* mom)
{
  //get all NueRecords from mom 
  //may have more than one per go since mom reads in a snarl's worth of data
  //so, this is a little more complicated than just asking for a NueRecord
  static Float_t total_pot = 0;
  static Int_t runno = -1;
  static Int_t snarlno = -1;
  static Int_t printpot = 1;
  TObject *obj=0;
  TIter objiter = mom->FragmentIter();
   while((obj=objiter.Next())){
      NueRecord *nr = static_cast<NueRecord *>(obj);
      if(nr){
	 MSG("CompareMD",Msg::kDebug)<<"Found a NueRecord in MOM"<<endl;
      }
      else{
	 MSG("CompareMD",Msg::kError)<<"Didn't find a NueRecord in MOM"<<endl;
	 continue;
      }
      MSG("CompareMD",Msg::kDebug)<<"Found a NueRecord "<<nr<<endl;

      SimFlag::SimFlag_t s = nr->GetHeader().GetVldContext().GetSimFlag();

//  need at least 1 sucessfully reconstructed event
   
     if(counter%1000==0){
         cout<<"On entry "<<counter<<endl;
     }
      counter++;
	  
      
      if(nr->GetHeader().GetEventNo()<0){
	continue;
      }

      if(PassesBeamCuts(nr)){
	if (s==SimFlag::kData){
	  if (runno != nr->GetHeader().GetRun()||
	      runno == nr->GetHeader().GetRun()&&
	      snarlno != nr->GetHeader().GetSnarl()){
	    total_pot += nr->bmon.bI;
	    //cout<<"total_pot "<<total_pot<<endl;
	    runno = nr->GetHeader().GetRun();
	    snarlno = nr->GetHeader().GetSnarl();
	  }
	}
	if (PassesCuts(nr)){
	  if(s==SimFlag::kData){
	    TString id = MakeIdString(nr);
	    static HistMan *hm = new HistMan("allcomp");
	    FillFromList(nr,id,hm,1.);
	  }
	  else if (s==SimFlag::kMC){
	    if (printpot) {
	      cout<<"Total POTs "<<total_pot<<endl;
	      printpot=0;
	    }	   
	    TString id = MakeIdString(nr);
	    static HistMan *hm = new HistMan("allcomp");
	    FillFromList(nr,id,hm,total_pot/(kNMCFiles*550*24.));
	    TString idall = "_m_all";
	    FillFromList(nr,idall,hm,total_pot/(kNMCFiles*550*24.));
	  }
	}
      }
/*

      //fill histograms      
          FillShwFit(hm, param,nr);
          FillHitCalc(hm, param,nr);
          FillAng(hm, param,nr);
          FillFracVar(hm, param,nr);
          FillMST(hm, param,nr);
          FillSREvent(hm, param,nr);
          FillSRTrack(hm, param,nr);
          FillSRShw(hm, param,nr);
          FillBMon(hm,param,nr);
          FillTruth(hm,param,nr);
          FillAnalysis(hm,param,nr);
      }
*/


   }
   return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}

////////////////////////////////////////////////////////////////////////
void CompareMD::EndJob()
{}


bool CompareMD::PassesBeamCuts(NueRecord* nr)
{
  //bool passes = true;   
  if (nr->GetHeader().GetVldContext().GetSimFlag()!=SimFlag::kData) return true;
  if (nr->bmon.bI<0.1) return false;      //beam intensity (in 1e12)
  if (nr->bmon.hbw<0.8) return false;     //horizontal beam width
  if (nr->bmon.hbw>1.1) return false;
  if (nr->bmon.vbw<0.9) return false;     //vertical beam width
  if (nr->bmon.vbw>1.6) return false;
  if (nr->bmon.nuTarZ<500) return false;  //low energy position
  if (nr->bmon.nuTarZ>1000) return false;
  return true;
}


bool CompareMD::PassesCuts(NueRecord* nr)
{
    bool passes = true;   

    //dont even bother if there is a long track
    if(nr->srtrack.planes > 25) passes = false;
    if (nr->srevent.pulseHeight<2e4) return false;

    if  (nr->anainfo.isFullyContained != 1 && nr->anainfo.isFullyContained!= -2) passes = false;

    if (nr->srevent.phMeu>150) passes = false;
    
    if (nr->srevent.hotch==1) return false;
    
    if ((TMath::Max(nr->srtrack.pulseHeight,nr->srshower.pulseHeight)<5000)) 
	passes = false;

//    if(nr->anainfo.inFiducialVolume != 1)
//        passes = false;

    //fiducial volume cuts
//    if ((pow(nr->srevent.vtxX-1.4885,2)+pow(nr->srevent.vtxY-0.1397,2))>0.81) return false;
//    if (nr->srevent.vtxZ<1) return false;
//    if (nr->srevent.vtxZ>6) return false;
    


/*     if(nr->anainfo.isFullyContained == 1 || 
  Cuts += "NueRecord.srevent.phMeu > 1 && NueRecord.srevent.phMeu < 400 && ";

  string NC = "NueRecord.mctrue.interactionType == 0";
  string CC = "NueRecord.mctrue.interactionType == 1";


   Cuts += "NueRecord.anainfo.isFullyContained == 1 && ";
   Cuts += "(NueRecord.anainfo.isFullyContained == 1 || NueRecord.anainfo.isFullyContained == -2 ) && ";
     if(j == 2)  Cuts += "";

     TString isNumu = "&& abs(NueRecord.mctrue.nuFlavor) == 14";
     TString isBNue = "&& abs(NueRecord.mctrue.nuFlavor) == 12";
     TString isNue  = "&& abs(NueRecord.mctrue.nonOscNuFlavor) == 14";
*/

   //only look at events for which mst actually gets calculated
   if(nr->GetHeader().GetTrackLength()>25)
       passes = false;

    return passes;
     
}


TString CompareMD::MakeIdString(NueRecord *nr)
{
  Detector::Detector_t d = nr->GetHeader().GetVldContext().GetDetector();
  SimFlag::SimFlag_t s = nr->GetHeader().GetVldContext().GetSimFlag();
  TString det;
  TString dm;
  TString type;

  if(d==Detector::kFar){
    det = "f";
  }
  else if(d==Detector::kNear){
    det = "n";
  }
  else{
    det="u";
  }
  
  if(s==SimFlag::kData){
    dm = "_d";
    TString id = dm;
    return id;
  }
  else if (s==SimFlag::kMC){
    dm = "_m";
    
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
	type = "u";
      }
    }
    else{
      type = "nc";
    }
    TString id = dm+ "_" + type;
    return id;
  }      
  else {
    dm = "_u";
    return dm;
  }

}


void CompareMD::FillFromList(NueRecord* nr, TString id, HistMan* hm, Float_t weight)
{
    Int_t pos = 0;
    Int_t size = varName.size();
    if(varName.size() == 0) return;
    TString hname;
    Int_t lastpos = -1;
	
    for(int i = 0; i < 5; i++){ 
//	   cout<<"starting loop "<<i<<"   "<<pos<<endl;
     if(varName[pos] == "fHeader.fSnarl") {
         if ( !ANtpDefVal::IsDefault(nr->GetHeader().GetSnarl())) {
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->GetHeader().GetSnarl(),weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fHeader.fRun") {
         if(!ANtpDefVal::IsDefault(nr->GetHeader().GetRun())){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->GetHeader().GetRun(),weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fHeader.fSubRun") {
         if(!ANtpDefVal::IsDefault(nr->GetHeader().GetSubRun())){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->GetHeader().GetSubRun(),weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fHeader.fEvtNo") {
         if(!ANtpDefVal::IsDefault(nr->GetHeader().GetEventNo())){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->GetHeader().GetEventNo(),weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fHeader.fEvents") {
         if(!ANtpDefVal::IsDefault(nr->GetHeader().GetEvents())){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->GetHeader().GetEvents(),weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fHeader.fTrackLength") {
         if(!ANtpDefVal::IsDefault(nr->GetHeader().GetTrackLength())){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->GetHeader().GetTrackLength(),weight);
          } pos++; if(pos >= size) return;
     }

//  End of fHeader
//  Start of Shwfit
 
if(varName[pos] == "shwfit.par_a") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.par_a)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.par_a,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.par_b") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.par_b)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.par_b,weight);
          } pos++; if(pos >= size) return;
     }

     if(varName[pos] == "shwfit.par_e0") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.par_e0)){
          hname = varName[pos]+id;
//	  cout<<"Filling "<<nr->shwfit.par_e0<<endl;
          hm->Fill1d(hname,nr->shwfit.par_e0,weight);
          } pos++; if(pos >= size) return;
     }
 if(varName[pos] == "shwfit.chisq") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.chisq)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.chisq,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.shwmax") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.shwmax)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.shwmax,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.shwmaxplane") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.shwmaxplane)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.shwmaxplane,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.conv") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.conv)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.conv,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.chisq_ndf") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.chisq_ndf)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.chisq_ndf,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.e0_pe_ratio") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.e0_pe_ratio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.e0_pe_ratio,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.caldet_comp") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.caldet_comp)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.caldet_comp,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.max_pe_plane") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.max_pe_plane)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.max_pe_plane,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.shwmaxplane_diff") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.shwmaxplane_diff)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.shwmaxplane_diff,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.u_asym_peak") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.u_asym_peak)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.u_asym_peak,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.u_asym_vert") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.u_asym_vert)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.u_asym_vert,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.u_molrad_peak") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.u_molrad_peak)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.u_molrad_peak,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.u_molrad_vert") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.u_molrad_vert)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.u_molrad_vert,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.u_mean") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.u_mean)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.u_mean,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.u_rms") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.u_rms)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.u_rms,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.u_skew") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.u_skew)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.u_skew,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.u_kurt") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.u_kurt)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.u_kurt,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.v_asym_peak") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.v_asym_peak)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.v_asym_peak,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.v_asym_vert") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.v_asym_vert)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.v_asym_vert,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.v_molrad_peak") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.v_molrad_peak)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.v_molrad_peak,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.v_molrad_vert") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.v_molrad_vert)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.v_molrad_vert,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.v_mean") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.v_mean)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.v_mean,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.v_rms") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.v_rms)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.v_rms,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.v_skew") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.v_skew)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.v_skew,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.v_kurt") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.v_kurt)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.v_kurt,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.uv_asym_peak") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.uv_asym_peak)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.uv_asym_peak,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.uv_asym_vert") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.uv_asym_vert)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.uv_asym_vert,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.uv_molrad_peak") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.uv_molrad_peak)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.uv_molrad_peak,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.uv_molrad_vert") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.uv_molrad_vert)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.uv_molrad_vert,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.uv_mean") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.uv_mean)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.uv_mean,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.uv_rms") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.uv_rms)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.uv_rms,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.uv_skew") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.uv_skew)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.uv_skew,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.uv_kurt") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.uv_kurt)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.uv_kurt,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "shwfit.uv_ratio") {
       if(!ANtpDefVal::IsDefault(nr->shwfit.uv_ratio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->shwfit.uv_ratio,weight);
          } pos++; if(pos >= size) return;
     }

//end of shwfitcalc

//start of hitcalc
if(varName[pos] == "hitcalc.fHitTotalEnergy") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitTotalEnergy)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitTotalEnergy,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitTransEnergy") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitTransEnergy)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitTransEnergy,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitLongEnergy") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitLongEnergy)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitLongEnergy,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitTransCMEnergy") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitTransCMEnergy)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitTransCMEnergy,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitTransEnergyRatio") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitTransEnergyRatio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitTransEnergyRatio,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitLongEnergyRatio") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitLongEnergyRatio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitLongEnergyRatio,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitTransLongEnergyRatio") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitTransLongEnergyRatio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitTransLongEnergyRatio,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitTransCMEnergyRatio") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitTransCMEnergyRatio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitTransCMEnergyRatio,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitFarMomBalance") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitFarMomBalance)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitFarMomBalance,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitPeakMomBalance") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitPeakMomBalance)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitPeakMomBalance,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitFarAngle") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitFarAngle)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitFarAngle,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "hitcalc.fHitPeakAngle") {
       if(!ANtpDefVal::IsDefault(nr->hitcalc.fHitPeakAngle)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->hitcalc.fHitPeakAngle,weight);
          } pos++; if(pos >= size) return;
     }

//end of hitcalc
//start of angcluster
if(varName[pos] == "angcluster.fACluRmsShwAxis") {
       if(!ANtpDefVal::IsDefault(nr->angcluster.fACluRmsShwAxis)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angcluster.fACluRmsShwAxis,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angcluster.fACluRmsZAxis") {
       if(!ANtpDefVal::IsDefault(nr->angcluster.fACluRmsZAxis)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angcluster.fACluRmsZAxis,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angcluster.fACluShwDirX") {
       if(!ANtpDefVal::IsDefault(nr->angcluster.fACluShwDirX)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angcluster.fACluShwDirX,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angcluster.fACluShwDirY") {
       if(!ANtpDefVal::IsDefault(nr->angcluster.fACluShwDirY)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angcluster.fACluShwDirY,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angcluster.fACluShwDirZ") {
       if(!ANtpDefVal::IsDefault(nr->angcluster.fACluShwDirZ)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angcluster.fACluShwDirZ,weight);
          } pos++; if(pos >= size) return;
     }
     if(varName[pos] == "angclusterfit.fACluFitParA") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitParA)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitParA,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitParB") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitParB)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitParB,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitParLongE0") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitParLongE0)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitParLongE0,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitShwMax") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitShwMax)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitShwMax,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitE0EnergyRatio") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitE0EnergyRatio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitE0EnergyRatio,weight);
          } pos++; if(pos >= size) return;
     }
     if(varName[pos] == "angclusterfit.fACluFitParL1") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitParL1)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitParL1,weight);
          } pos++; if(pos >= size) return;
     }
     if(varName[pos] == "angclusterfit.fACluFitParL2") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitParL2)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitParL2,weight);
          } pos++; if(pos >= size) return;
     }
     if(varName[pos] == "angclusterfit.fACluFitParC12") {
          if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitParC12)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitParC12,weight);
          } pos++; if(pos >= size) return;
     }
     if(varName[pos] == "angclusterfit.fACluFitParTransE0") {
          if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitParTransE0)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitParTransE0,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitTransChiSq") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitTransChiSq)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitTransChiSq,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitTransNDF") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitTransNDF)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitTransNDF,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitLongNDF") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitLongNDF)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitLongNDF,weight);
          } pos++; if(pos >= size) return;
     }    
     if(varName[pos] == "angclusterfit.fACluFitMolRadPeak") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitMolRadPeak)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitMolRadPeak,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitMolRadVert") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitMolRadVert)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitMolRadVert,weight);
          } pos++; if(pos >= size) return;     
     }if(varName[pos] == "angclusterfit.fACluFitMean") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitMean)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitMean,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitRMS") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitRMS)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitRMS,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitSkew") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitSkew)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitSkew,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitKurt") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitKurt)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitKurt,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitLongChiSq") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitLongChiSq)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitLongChiSq,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitLongConv") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitLongConv)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitLongConv,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "angclusterfit.fACluFitTransConv") {
       if(!ANtpDefVal::IsDefault(nr->angclusterfit.fACluFitTransConv)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->angclusterfit.fACluFitTransConv,weight);
          } pos++; if(pos >= size) return;
     }

//end of ang
//start of mst
//     if(nr->mstvars.enn1 >3 && nr->mstvars.onn1>3){
           
     if(varName[pos] == "mstvars.enmsg") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.enmsg)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.enmsg,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.ew1") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.ew1)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.ew1,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.enn1") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.enn1)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.enn1,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.esm1") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.esm1)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.esm1,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.ewtot") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.ewtot)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.ewtot,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.enntot") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.enntot)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.enntot,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.esmtot") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.esmtot)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.esmtot,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.e4w") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.e4w)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.e4w,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.e4sm") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.e4sm)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.e4sm,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.e4nn") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.e4nn)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.e4nn,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.eb1") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.eb1)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.eb1,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.eb25") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.eb25)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.eb25,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.enbranch") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.enbranch)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.enbranch,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.onmsg") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.onmsg)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.onmsg,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.ow1") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.ow1)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.ow1,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.onn1") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.onn1)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.onn1,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.osm1") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.osm1)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.osm1,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.owtot") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.owtot)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.owtot,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.onntot") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.onntot)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.onntot,weight);
          } pos++; if(pos >= size) return;
     }
     if(varName[pos] == "mstvars.osmtot") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.osmtot)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.osmtot,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.o4w") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.o4w)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.o4w,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.o4sm") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.o4sm)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.o4sm,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.o4nn") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.o4nn)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.o4nn,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.ob1") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.ob1)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.ob1,weight);
          } pos++; if(pos >= size) return;
     }
     if(varName[pos] == "mstvars.ob25") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.ob25)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.ob25,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.onbranch") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.onbranch)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.onbranch,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.eallw1") {
        for(int i=0;i<nr->mstvars.enn1;i++){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.eallw1[i],weight);           
        }
        pos++;          
       if(pos >= size) return;
     }if(varName[pos] == "mstvars.oallw1") {
        for(int i=0;i<nr->mstvars.onn1;i++){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.oallw1[i],weight); 
        }
        pos++;       if(pos >= size) return;
     }if(varName[pos] == "mstvars.eallm1") {
        for(int i=0;i<nr->mstvars.enn1;i++){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.eallm1[i],weight);
        }
        pos++;          if(pos >= size) return;
     }if(varName[pos] == "mstvars.oallm1") {
        for(int i=0;i<nr->mstvars.onn1;i++){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.oallm1[i],weight); 
        }
        pos++;          if(pos >= size) return;
     }if(varName[pos] == "mstvars.eeprob") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.eeprob)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.eeprob,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.oeprob") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.oeprob)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.oeprob,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.ealpha") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.ealpha)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.ealpha,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.oalpha") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.oalpha)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.oalpha,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.ebeta") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.ebeta)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.ebeta,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "mstvars.obeta") {
       if(!ANtpDefVal::IsDefault(nr->mstvars.obeta)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->mstvars.obeta,weight);
          } pos++; if(pos >= size) return;
     }
//     }
//end of mst vars
//start of fracvars

     if(varName[pos] == "fracvars.fract_1_plane") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_1_plane)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_1_plane,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_2_planes") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_2_planes)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_2_planes,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_3_planes") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_3_planes)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_3_planes,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_4_planes") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_4_planes)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_4_planes,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_5_planes") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_5_planes)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_5_planes,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_6_planes") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_6_planes)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_6_planes,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_2_counters") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_2_counters)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_2_counters,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_4_counters") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_4_counters)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_4_counters,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_6_counters") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_6_counters)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_6_counters,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_8_counters") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_8_counters)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_8_counters,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_10_counters") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_10_counters)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_10_counters,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_12_counters") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_12_counters)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_12_counters,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.fract_road") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.fract_road)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.fract_road,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.shw_max") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.shw_max)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.shw_max,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.shw_nstp") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.shw_nstp)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.shw_nstp,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.shw_npl") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.shw_npl)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.shw_npl,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "fracvars.pid") {
       if(!ANtpDefVal::IsDefault(nr->fracvars.pid)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->fracvars.pid,weight);
          } pos++; if(pos >= size) return;
     }


//end of fracvar
//start of anainfo

     if(varName[pos] == "anainfo.inFiducialVolume") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.inFiducialVolume)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.inFiducialVolume,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.isFullyContained") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.isFullyContained)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.isFullyContained,weight);
          } pos++; if(pos >= size) return;
/*     }if(varName[pos] == "anainfo.passesCuts") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.passesCuts)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.passesCuts,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.pass") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.pass)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.pass,weight);
          } pos++; if(pos >= size) return;
//     }if(varName[pos] == "anainfo.isNC") {
//       if(!ANtpDefVal::IsDefault(nr->anainfo.isNC)){
//          hname = varName[pos]+id;
//          hm->Fill1d(hname,nr->anainfo.isNC,weight);
//          } pos++; if(pos >= size) return;
//     }if(varName[pos] == "anainfo.isCC") {
//       if(!ANtpDefVal::IsDefault(nr->anainfo.isCC)){
//          hname = varName[pos]+id;
//          hm->Fill1d(hname,nr->anainfo.isCC,weight);
//          } pos++; if(pos >= size) return;
//     }if(varName[pos] == "anainfo.separationParameterCut") {
//       if(!ANtpDefVal::IsDefault(nr->anainfo.separationParameterCut)){
//          hname = varName[pos]+id;
//          hm->Fill1d(hname,nr->anainfo.separationParameterCut,weight);
//          } pos++; if(pos >= size) return;
//     }if(varName[pos] == "anainfo.separationParameter") {
//       if(!ANtpDefVal::IsDefault(nr->anainfo.separationParameter)){
//          hname = varName[pos]+id;
//          hm->Fill1d(hname,nr->anainfo.separationParameter,weight);
//          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoNuEnergy") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoNuEnergy)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoNuEnergy,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoMuEnergy") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoMuEnergy)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoMuEnergy,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoShowerEnergy") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoShowerEnergy)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoShowerEnergy,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoQENuEnergy") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoQENuEnergy)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoQENuEnergy,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoQEQ2") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoQEQ2)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoQEQ2,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoHadronicY") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoHadronicY)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoHadronicY,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoNuDCos") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoNuDCos)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoNuDCos,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoMuDCosZVtx") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoMuDCosZVtx)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoMuDCosZVtx,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoVtxX") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoVtxX)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoVtxX,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoVtxY") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoVtxY)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoVtxY,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoVtxZ") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoVtxZ)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoVtxZ,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoEventLength") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoEventLength)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoEventLength,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoTrackLength") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoTrackLength)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoTrackLength,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoTrackMomentum") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoTrackMomentum)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoTrackMomentum,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoTrackRange") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoTrackRange)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoTrackRange,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "anainfo.recoSigmaQoverP") {
       if(!ANtpDefVal::IsDefault(nr->anainfo.recoSigmaQoverP)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->anainfo.recoSigmaQoverP,weight);
          } pos++; if(pos >= size) return;
  */   }

//end of anainfo
//start of srevent
//     cout<<varName[pos]<<"  "<<pos<<endl;

     if(varName[pos] == "srevent.event") {
       if(!ANtpDefVal::IsDefault(nr->srevent.index)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.index,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.pulseHeight") {
       if(!ANtpDefVal::IsDefault(nr->srevent.pulseHeight)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.pulseHeight,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.energyGeV") {
       if(!ANtpDefVal::IsDefault(nr->srevent.energyGeV)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.energyGeV,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.begPlane") {
       if(!ANtpDefVal::IsDefault(nr->srevent.begPlane)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.begPlane,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.endPlane") {
       if(!ANtpDefVal::IsDefault(nr->srevent.endPlane)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.endPlane,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.planes") {
       if(!ANtpDefVal::IsDefault(nr->srevent.planes)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.planes,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.totalStrips") {
       if(!ANtpDefVal::IsDefault(nr->srevent.totalStrips)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.totalStrips,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.passStrips") {
       if(!ANtpDefVal::IsDefault(nr->srevent.passStrips)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.passStrips,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.showers") {
       if(!ANtpDefVal::IsDefault(nr->srevent.showers)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.showers,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.tracks") {
       if(!ANtpDefVal::IsDefault(nr->srevent.tracks)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.tracks,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.vtxX") {
       if(!ANtpDefVal::IsDefault(nr->srevent.vtxX)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.vtxX,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.vtxY") {
       if(!ANtpDefVal::IsDefault(nr->srevent.vtxY)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.vtxY,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.vtxZ") {
       if(!ANtpDefVal::IsDefault(nr->srevent.vtxZ)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.vtxZ,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.timeLength") {
       if(!ANtpDefVal::IsDefault(nr->srevent.timeLength)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.timeLength,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.phMeu") {
       if(!ANtpDefVal::IsDefault(nr->srevent.phMeu)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srevent.phMeu,weight);
        } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.phNueGeV") {
        if(!ANtpDefVal::IsDefault(nr->srevent.phNueGeV)){
           hname = varName[pos]+id;
           hm->Fill1d(hname,nr->srevent.phNueGeV,weight);
        } pos++; if(pos >= size) return;
     }if(varName[pos] == "srevent.triggerPass") {
        if(!ANtpDefVal::IsDefault(nr->srevent.triggerPass)){
           hname = varName[pos]+id;
           hm->Fill1d(hname,nr->srevent.triggerPass,weight);
        } pos++; if(pos >= size) return;
     }
     

//end of srevent
//start of srshw

     if(varName[pos] == "srshower.planes") {
       if(!ANtpDefVal::IsDefault(nr->srshower.planes)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.planes,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.totalStrips") {
       if(!ANtpDefVal::IsDefault(nr->srshower.totalStrips)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.totalStrips,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.begPlane") {
       if(!ANtpDefVal::IsDefault(nr->srshower.begPlane)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.begPlane,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.endPlane") {
       if(!ANtpDefVal::IsDefault(nr->srshower.endPlane)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.endPlane,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.vtxX") {
       if(!ANtpDefVal::IsDefault(nr->srshower.vtxX)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.vtxX,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.vtxY") {
       if(!ANtpDefVal::IsDefault(nr->srshower.vtxY)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.vtxY,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.vtxZ") {
       if(!ANtpDefVal::IsDefault(nr->srshower.vtxZ)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.vtxZ,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.dcosX") {
       if(!ANtpDefVal::IsDefault(nr->srshower.dcosX)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.dcosX,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.dcosY") {
       if(!ANtpDefVal::IsDefault(nr->srshower.dcosY)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.dcosY,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.dcosZ") {
       if(!ANtpDefVal::IsDefault(nr->srshower.dcosZ)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.dcosZ,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.pulseHeight") {
       if(!ANtpDefVal::IsDefault(nr->srshower.pulseHeight)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.pulseHeight,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.stripRatio") {
       if(!ANtpDefVal::IsDefault(nr->srshower.stripRatio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.stripRatio,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.planeRatio") {
       if(!ANtpDefVal::IsDefault(nr->srshower.planeRatio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.planeRatio,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srshower.pulseHeightRatio") {
       if(!ANtpDefVal::IsDefault(nr->srshower.pulseHeightRatio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srshower.pulseHeightRatio,weight);
          } pos++; if(pos >= size) return;
     }

//end of srshower

//start of srtrack
if(varName[pos] == "srtrack.planes") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.planes)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.planes,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.totalStrips") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.totalStrips)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.totalStrips,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.pulseHeight") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.pulseHeight)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.pulseHeight,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.fitMomentum") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.fitMomentum)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.fitMomentum,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.rangeMomentum") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.rangeMomentum)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.rangeMomentum,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.sigmaQoverP") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.sigmaQoverP)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.sigmaQoverP,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.begPlane") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.begPlane)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.begPlane,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.endPlane") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.endPlane)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.endPlane,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.length") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.length)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.length,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.vtxX") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.vtxX)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.vtxX,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.vtxY") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.vtxY)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.vtxY,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.vtxZ") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.vtxZ)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.vtxZ,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.dcosXVtx") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.dcosXVtx)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.dcosXVtx,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.dcosYVtx") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.dcosYVtx)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.dcosYVtx,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.dcosZVtx") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.dcosZVtx)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.dcosZVtx,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.endX") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.endX)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.endX,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.endY") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.endY)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.endY,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.endZ") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.endZ)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.endZ,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.dcosXEnd") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.dcosXEnd)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.dcosXEnd,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.dcosYEnd") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.dcosYEnd)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.dcosYEnd,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.dcosZEnd") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.dcosZEnd)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.dcosZEnd,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.reducedChi2") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.reducedChi2)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.reducedChi2,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.passedFit") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.passedFit)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.passedFit,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.traceVtx") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.traceVtx)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.traceVtx,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.traceVtxZ") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.traceVtxZ)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.traceVtxZ,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.traceEnd") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.traceEnd)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.traceEnd,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.traceEndZ") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.traceEndZ)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.traceEndZ,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.trklikePlanes") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.trklikePlanes)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.trklikePlanes,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.trklikeRatio") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.trklikeRatio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.trklikeRatio,weight);
          } pos++; if(pos >= size) return;
     }if(varName[pos] == "srtrack.pulseHeightRatio") {
       if(!ANtpDefVal::IsDefault(nr->srtrack.pulseHeightRatio)){
          hname = varName[pos]+id;
          hm->Fill1d(hname,nr->srtrack.pulseHeightRatio,weight);
          } pos++; if(pos >= size) return;
     }


//end of srtrack


     if(pos == lastpos) pos++;  //it couldnt find it in the list skip it
          lastpos = pos;  
   }

   return;     
}     
           
           
        
const Registry& CompareMD::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================

  static Registry r; // Default configuration for module

  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  // Set values in configuration
  r.UnLockValues();
  r.Set("NMCFiles","1");
  r.LockValues();

  return r;
}

//......................................................................

void CompareMD::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================

  int imps;
  if(r.Get("NMCFiles",imps)) { kNMCFiles=imps;}

}


         
