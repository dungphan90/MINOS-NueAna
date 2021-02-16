#include "NueAna/Extrapolation/NueSystematic.h"
#include "NueAna/Extrapolation/Background.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueAna/NueStandard.h"
#include "MCReweight/ReweightHelpers.h"
#include "MCReweight/MCReweight.h"
#include "Registry/Registry.h"
#include <stdlib.h>
#include <string>
#include <iostream>
#include "TFile.h"
#include "TH1D.h"
#include "TMath.h"

SKZPWeightCalculator* NueSystematic::skzp = 0;

NueSystematic::NueSystematic(std::string name) :
  //  skzpcfg("PiMinus_CedarDaikon"),
  skzpcfg("DetXs"),
  fTempDouble(-9999)
{
  sprintf(fName,"%s",name.c_str());
  this->Init();
}
 
NueSystematic::~NueSystematic()
{
  fSysList.clear();
  //delete fNWC->GetReweightConfig();
  //delete fNWC;
}

void NueSystematic::Init()
{
  fTheta12 = 0.59365; 
  fTheta23 = 0.785398; 
  fTheta13 = 0.19885;
  
  fDeltaMSq12 = 8.0e-5;  
  fDeltaMSq23 = 2.32e-3;

  fDeltaCP = 0; 
  fMassHierarchy = 1;

  if(skzp == 0){
     skzp = new SKZPWeightCalculator(skzpcfg, true);
  }

}

void NueSystematic::SetSKZPParams(std::string cfg)
{
  skzpcfg = cfg;
}

std::string NueSystematic::GetSKZPParams()
{
  return skzpcfg;
}

void NueSystematic::SetOscParams(Double_t theta12,Double_t theta23,Double_t theta13,
				 Double_t deltam12,Double_t deltam23,Double_t deltaCP,
				 Int_t massH)
{
  fTheta12 = theta12; fTheta23 = theta23; fTheta13 = theta13;
  fDeltaMSq12 = deltam12; fDeltaMSq23 = deltam23;
  fDeltaCP = deltaCP; fMassHierarchy = massH;

  Double_t par[9] = {0};
  par[OscPar::kL] = 735.0;
  par[OscPar::kTh23] = fTheta23;
  par[OscPar::kTh12] = fTheta12;
  par[OscPar::kTh13] = fTheta13; // TMath::ASin(TMath::Sqrt(ss2th13))/2.;
  par[OscPar::kDeltaM23] = massH*fDeltaMSq23;
  par[OscPar::kDeltaM12] = fDeltaMSq12;
  par[OscPar::kDensity] = 2.75; //standard rock density
  par[OscPar::kDelta] = deltaCP;
  par[OscPar::kNuAntiNu] = 1;
  fOscCalc.SetOscParam(par);
}

void NueSystematic::SetOscParams(double *par)
{

  fOscCalc.SetOscParam(par);
}


void NueSystematic::GetOscParams(Double_t &theta12,Double_t &theta23,Double_t &theta13,
				 Double_t &deltam12,Double_t &deltam23,Double_t &deltaCP,
				 Int_t &massH)
{
  fTheta12 = theta12; fTheta23 = theta23; fTheta13 = theta13;
  fDeltaMSq12 = deltam12; fDeltaMSq23 = deltam23;
  fDeltaCP = deltaCP; fMassHierarchy = massH;
  
  Double_t par[10];
  fOscCalc.GetOscParam(par);

  par[OscPar::kL] = 735.0;
  theta23 = par[OscPar::kTh23];
  theta12 = par[OscPar::kTh12];
  theta13 = par[OscPar::kTh13];
  deltam23 = par[OscPar::kDeltaM23];
  deltam12 = par[OscPar::kDeltaM12];
  deltaCP = par[OscPar::kDelta];

  massH = 1;
  if(deltam23 < 0){ massH = -1; deltam23 = -deltam23;}
}

Double_t NueSystematic::GetSysValue(Systematic::Systematic_t sys)
{
  std::map<Systematic::Systematic_t,Double_t>::iterator beg = fSysList.begin();
  std::map<Systematic::Systematic_t,Double_t>::iterator end = fSysList.end();
  while (beg!=end) {
    if(beg->first==sys) return beg->second;
    beg++;
  }
  return Systematic::GetDefaultValue(sys);
}

void NueSystematic::SetSysValue(Systematic::Systematic_t sys, Double_t val)
{
  std::map<Systematic::Systematic_t,Double_t>::iterator beg = fSysList.begin();
  std::map<Systematic::Systematic_t,Double_t>::iterator end = fSysList.end();

  while (beg!=end) {
    if(beg->first==sys){ beg->second = val;  return; }
    beg++;
  }
}

Double_t NueSystematic::UpdateRecord(NueRecord *rec,Selection::Selection_t sel)
{  
  Double_t totWeight = 1;
  Double_t enShift = 0;
  Double_t pidShift = 0;
  Double_t trkLength = 0;
  Double_t trkLike = 0;
  double shwEnShift = 0;

  std::map<Systematic::Systematic_t,Double_t>::iterator beg = fSysList.begin();
  std::map<Systematic::Systematic_t,Double_t>::iterator end = fSysList.end();
  while (beg!=end) {
    Systematic::Systematic_t sys = beg->first;
    Double_t val = beg->second;
    switch (sys) {
    case Systematic::kNorm     : totWeight *= this->DoNormCalc(rec,val);
    case Systematic::kEMCalib  : enShift   += this->DoCalibShift(rec,sys,val);  break;
    case Systematic::kHadCalib : enShift   += this->DoCalibShift(rec,sys,val);  break;
    case Systematic::kRelCalib : enShift   += this->DoCalibShift(rec,sys,val);  break;
    case Systematic::kMA_QE    : totWeight *= this->DoNeugenCalc(rec,sys,val);  break;
    case Systematic::kMA_RES   : totWeight *= this->DoNeugenCalc(rec,sys,val);  break;
    case Systematic::kKNO      : totWeight *= this->DoNeugenCalc(rec,sys,val);  break;
    case Systematic::kTrkLike  : trkLike   += val;                              break;
    case Systematic::kTrkPlane : trkLength += val;                              break;
    case Systematic::kPIDShift : pidShift  += val;                              break;
    case Systematic::kSKZP     : totWeight *= this->DoSKZPCalc(rec,val);        break;
    case Systematic::kOscProb  : totWeight *= this->DoOscCalc(rec,val);         break;
    case Systematic::kShwDev   : totWeight *= this->DoShwDevCalc(rec,val,sel);  break;
    case Systematic::kTauProd  : totWeight *= this->DoTauProd(rec,val);         break;
    case Systematic::kPIDSkew  : totWeight *= this->DoPIDSkew(rec,val,sel);     break;
    case Systematic::kNCScale  : totWeight *= this->DoNCScale(rec,val, sel);           break;
    case Systematic::kCCShwE   : shwEnShift += this->DoCCShwEnergyScale(rec,val,sel);  break;
    default:                                                                    break;
    }
    beg++;
//    std::cout<<Systematic::AsString(sys)<<"  "<<val<<"  "<<totWeight<<std::endl;
  }

  rec->srevent.phNueGeV         *= 1+enShift;
  rec->srshower.phCCGeV         *= 1+shwEnShift;

  rec->subshowervars.pid += pidShift;
  rec->ann.pid_30inp            += pidShift;
  rec->ann.pid_6inp             += pidShift;
  rec->ann.pid_11inp            += pidShift;
  rec->ann.pid_11inp_daikon04   += pidShift;
  rec->mcnnv.mcnn_var1      += pidShift;

//  rec->mdadiscrim.fMdaPIDnue += pidShift;
//  rec->fracvars.pid          += pidShift;
//  rec->fracvars.pid1         += pidShift;
  rec->srtrack.endPlane      += Int_t(trkLength);
  rec->srtrack.trklikePlanes += Int_t(trkLike);

  return totWeight*rec->fluxweights.totbeamweight;
}

Double_t NueSystematic::UpdateRecord(NueRecord *rec,Selection::Selection_t sel,
      Background::Background_t bg)
{
  Double_t totWeight = 1;
  Double_t enShift = 0;
  Double_t pidShift = 0;
  Double_t trkLength = 0;
  Double_t trkLike = 0;
  double shwEnShift = 0;

  //Need to Turn off oscillatons
  Double_t OscSysVal = this->GetSysValue(Systematic::kOscProb);
  this->SetSysValue(Systematic::kOscProb, 0);

  std::map<Systematic::Systematic_t,Double_t>::iterator beg = fSysList.begin();
  std::map<Systematic::Systematic_t,Double_t>::iterator end = fSysList.end();
  while (beg!=end) {
    Systematic::Systematic_t sys = beg->first;
    Double_t val = beg->second;
    switch (sys) {
    case Systematic::kNorm     : totWeight *= this->DoNormCalc(rec,val);
    case Systematic::kEMCalib  : enShift   += this->DoCalibShift(rec,sys,val);  break;
    case Systematic::kHadCalib : enShift   += this->DoCalibShift(rec,sys,val);  break;
    case Systematic::kRelCalib : enShift   += this->DoCalibShift(rec,sys,val);  break;
    case Systematic::kMA_QE    : totWeight *= this->DoNeugenCalc(rec,sys,val);
break;
    case Systematic::kMA_RES   : totWeight *= this->DoNeugenCalc(rec,sys,val);
break;
    case Systematic::kKNO      : totWeight *= this->DoNeugenCalc(rec,sys,val); break;
    case Systematic::kTrkLike  : trkLike   += val;                              break;
    case Systematic::kTrkPlane : trkLength += val;                              break;
    case Systematic::kPIDShift : pidShift  += val;                              break;
    case Systematic::kSKZP     : totWeight *= this->DoSKZPCalc(rec,val);        break;
    case Systematic::kOscProb  : totWeight *= this->DoOscCalc(rec,val);         break;
    case Systematic::kShwDev   : totWeight *= this->DoShwDevCalc(rec,val,sel);  break;
    case Systematic::kTauProd  : totWeight *= this->DoTauProd(rec,val);         break;
    case Systematic::kPIDSkew  : totWeight *= this->DoPIDSkew(rec,val,sel);     break;
    case Systematic::kNCScale  : totWeight *= this->DoNCScale(rec,val, sel);           break;
    case Systematic::kCCShwE   : shwEnShift += this->DoCCShwEnergyScale(rec,val,sel);  break;
    default:                                                                    break;
    }
    beg++;
  }

  this->SetSysValue(Systematic::kOscProb, OscSysVal);

  totWeight *= this->GetAppearanceWeight(rec, bg); 

  rec->srevent.phNueGeV      *= 1+enShift;
  rec->srshower.phCCGeV      *= 1+shwEnShift;
  rec->subshowervars.pid  += pidShift;
  rec->ann.pid_30inp             += pidShift;
  rec->ann.pid_6inp             += pidShift;
  rec->ann.pid_11inp            += pidShift;
  rec->ann.pid_11inp_daikon04   += pidShift;
  rec->mcnnv.mcnn_var1      += pidShift;

//  rec->mdadiscrim.fMdaPIDnue += pidShift;
//  rec->fracvars.pid          += pidShift;
//  rec->fracvars.pid1         += pidShift;
  rec->srtrack.endPlane      += Int_t(trkLength);
  rec->srtrack.trklikePlanes += Int_t(trkLike);

  return totWeight*rec->fluxweights.totbeamweight;
}


Double_t NueSystematic::GetAppearanceWeight(NueRecord *rec, Background::Background_t bg)
{
  if(bg == Background::kNueCC || bg == Background::kNuTauCC)  rec->mctrue.nonOscNuFlavor = 14;
  if(bg == Background::kNueCC) rec->mctrue.nuFlavor = 12;
  if(bg == Background::kNuTauCC) rec->mctrue.nuFlavor = 16;
  if(bg == Background::kNuTauCC) rec->mctrue.resonanceCode = 1001;
  if(rec->mctrue.nonOscNuFlavor < 0) rec->mctrue.nuFlavor *= -1;
  std::map<Systematic::Systematic_t,Double_t>::iterator beg = fSysList.begin();
  std::map<Systematic::Systematic_t,Double_t>::iterator end = fSysList.end();

  double totWeight = 1.0;
                                                                                                                    
  while (beg!=end) {
    Systematic::Systematic_t sys = beg->first;
    Double_t val = beg->second;
    switch (sys) {
    case Systematic::kOscProb  : totWeight *= this->DoOscCalc(rec,val);         break;
    case Systematic::kTauProd  : totWeight *= this->DoTauProd(rec,val);         break;
    default:                                                                    break;
    }
    beg++;
  }
                                                                                                                    
  return totWeight;
}
                                                                                                       
Double_t NueSystematic::DoSKZPCalc(NueRecord *record,Double_t val)
{
  // val < -999 -> unweight to base MC (use -1000)
  // val == 0 -> Do nothing - > standard SKZP
  // val > 0 -> distort spectrum by val*1sigma change
  double baseweight = record->fluxweights.totskzpweight;
  if(val < -999){
     if(baseweight <= 0) return 0;
     else return 1.0/baseweight;
  }
  else{
    if(val != 0 ) {
      if(record->mctrue.nuFlavor<-9998) return 1;
      float energy = record->mctrue.nuEnergy;
      int inu = record->mctrue.nonOscNuFlavor;
                                                                                
      Detector::Detector_t det = record->GetHeader().GetVldContext().GetDetector();
      //the 2 here is the IBEAM value of L010185
      // so for numu/antinumu we have the whole envelope
      // for nue we have the HadPrd + 1.77% envelope from Bob and Masaki
      // error returned here is the weight to change the spectrum up by 1 sigma
                                                                                
      double err = 1.0;
                                                                                
      if(TMath::Abs(inu) == 14){
         err = skzp->GetFluxError(det,2,inu,energy,SKZPWeightCalculator::kTotalErrorAfterTune);
         err -= 1.0;
      }
      if(TMath::Abs(inu) == 12){
         err = skzp->GetFluxError(det,2,inu,energy,SKZPWeightCalculator::kHadProdAfterTune);
         err = TMath::Sqrt((err - 1)*(err - 1) + 0.0177*0.0177);
      }
                                                                                
      double weight = 1.0 + err*val;
                                                                                
      return weight;
    }
  }

  return 1;
}

Double_t NueSystematic::DoOscCalc(NueRecord *record,Double_t val)
{
  Double_t osc_prob = 1;
  if(record->mctrue.nuFlavor<-9998) return osc_prob;
  Double_t L = 735.0;
  if(record->GetHeader().GetVldContext().GetDetector()==Detector::kNear){
     L = 1.0;  return 1.0;
  }

  fOscCalc.SetOscParam(OscPar::kNuAntiNu, 1);
  if(record->mctrue.nonOscNuFlavor < 0)   fOscCalc.SetOscParam(OscPar::kNuAntiNu, -1);                                                                                                                    
  double ue32 = TMath::Sin(fTheta13);
  ue32 *= ue32;
 

  if( int(val+0.5) ==1 ) {
    osc_prob = fOscCalc.Oscillate(record->mctrue.nuFlavor,
					record->mctrue.nonOscNuFlavor,
					record->mctrue.nuEnergy);
  }
  else if( int(val+0.5) == 2 )
    osc_prob = fOscCalc.Oscillate(record->mctrue.nuFlavor,
					      record->mctrue.nonOscNuFlavor,
					      record->mctrue.nuEnergy);

  return osc_prob;
}

Double_t NueSystematic::DoShwDevCalc(NueRecord *record,Double_t val, Selection::Selection_t sel)
{
  // Note for later usage:
  //  there was also a change in the nutau xsec carrot->Daikon, really this 
  //   should be disentangled, but going to ignore it for now
  //     E      r (c/d)      
  //    2.25       0
  //    2.75       0
  //    3.25       0
  //    3.75    1.05014
  //    4.25    1.0897
  //    4.75    1.05697
  //    5.25    1.07126
  //    5.75    1.09773
  //    6.25    1.13736
  //    6.75    1.19135
  //    7.25    1.23808
  //    7.75    1.27554
  //    8.25    1.30237
  //    8.75    1.32485
  //    9.25    1.34172
  //    9.75    1.35351
  //   10.25    1.36185
  //   10.75    1.36765
  //   11.25    1.37102
  //   11.75    1.37229
  if(sel == Selection::kCC) return 1;

  if(int(val+0.5)==1) {
    //weight based on ratio of N(F)D_daikon/N(F)D_carrot 

    const int NSEL = 10;
    static vector<TH1D*> ND_NC;
    static vector<TH1D*> ND_NUMU;
    static vector<TH1D*> ND_BNUE;

    static vector<TH1D*> FD_NC;
    static vector<TH1D*> FD_NUMU;
    static vector<TH1D*> FD_BNUE;
    static vector<TH1D*> FD_NUE;
    static vector<TH1D*> FD_NUTAU;

    if(ND_NC.size() == 0){

       for(int i = 0; i < NSEL; i++){
           TH1D* dum = 0;
           ND_NC.push_back(dum);
           ND_NUMU.push_back(dum);
           ND_BNUE.push_back(dum);
           FD_NC.push_back(dum);
           FD_NUMU.push_back(dum);
           FD_BNUE.push_back(dum);
           FD_NUE.push_back(dum);
           FD_NUTAU.push_back(dum);
       }

       const int NDET = 2;
       const int NUC = 5;
       const int NCUT = 4;
       string det[2]  = {"Near", "Far"};
       string nuC[5]  = {"numu", "nc", "bnue", "nue", "nutau"};
//       string cut[4]  = {"2", "6", "8", "7"};
       string name[NCUT] = {"Fid", "Presel", "ANN", "SSPID"};

       TFile Input("NueAna/data/CarrotDaikonWeights.root", "READ");
       TH1D* temp;
       TH1D* temp2;
     
       for(int i = 0; i < NDET; i++){
         for(int j = 0; j < NUC; j++){
          if(i < 1 && j >= 3) continue;
          for(int k = 0; k < NCUT; k++){
           TString id = det[i]+nuC[j]+name[k];
           Input.GetObject(id, temp);
           if(temp != 0){
             temp2 = (TH1D*) temp->Clone("ratio");
             temp2->SetDirectory(0);
 
             int iSel = (int) Selection::StringToEnum(name[k].c_str());
                       
             if(i+1 ==Detector::kFar ){
               if(j==1)  FD_NC[iSel] = temp2;
               if(j==0)  FD_NUMU[iSel] = temp2;
               if(j==4)  FD_NUTAU[iSel] = temp2;
               if(j==2)  FD_BNUE[iSel] = temp2;
               if(j==3)  FD_NUE[iSel] = temp2;
             }
             if(i+1 ==Detector::kNear ){
               if(j==1)  ND_NC[iSel] = temp2;
               if(j==0)  ND_NUMU[iSel] = temp2;
               if(j==2)  ND_BNUE[iSel] = temp2;
             }
          }
          else
            cout<<"Unable to Load: "<<id<<endl;
         }
       }
     }
    }

    Background::Background_t bg = 
      Background::TranslateFromMC(record->mctrue.interactionType,
				  record->mctrue.nuFlavor,
				  record->mctrue.nonOscNuFlavor);

    Double_t recoEnergy = record->srevent.phNueGeV;

    TH1D* temp = 0;
    int iSel = (int) sel;
    
    if(record->GetHeader().GetVldContext().GetDetector()==Detector::kFar){
      if(bg==Background::kNC)  temp = FD_NC[iSel];
      if(bg==Background::kNuMuCC) temp = FD_NUMU[iSel];
      if(bg==Background::kNuTauCC) temp = FD_NUTAU[iSel];
      if(bg==Background::kBNueCC) temp = FD_BNUE[iSel];
      if(bg==Background::kNueCC) temp = FD_NUE[iSel];
    }
    if(record->GetHeader().GetVldContext().GetDetector()==Detector::kNear){
      if(bg==Background::kNC) temp = ND_NC[iSel];
      if(bg==Background::kNuMuCC) temp = ND_NUMU[iSel];
      if(bg==Background::kBNueCC) temp = ND_BNUE[iSel];
    }

    if(temp == 0) cout<<"Massive loading error "<<bg<<"  "<<record->GetHeader().GetVldContext().GetDetector()<<endl;
    for(int i=0;i<temp->GetNbinsX();i++){
      if(recoEnergy>=temp->GetBinLowEdge(i) && recoEnergy< temp->GetBinLowEdge(i+1)) {
         return temp->GetBinContent(i);
      }
    }

    return 1;
  }
  if(int(val+0.5)==2){
    //weight using MODBYRS4
    Double_t weight = this->DoNeugenCalc(record,Systematic::kShwDev,4);
    if(record->xsecweights.xsecweight>0) weight/=record->xsecweights.xsecweight;
    return weight;
  }
  if(int(val+0.5)==3){
    //weight based on piZero energy reweighting (carrot->daikon)
    Double_t piZeroEnergy[8] = {0.0,0.75,1.5,2.25,3.0,3.75,4.5,100};
    Double_t weight[8] = {1.19675,1.06378,2.00025,1.91159,
			  1.61631,1.2359,1.84444,1.00};
    Double_t pi0E = record->shi.epi0;
    for(int i=0;i<7;i++){
      if(pi0E>=piZeroEnergy[i] && pi0E<piZeroEnergy[i+1])
	return weight[i];
    }
  }
  return 1;
}


Double_t NueSystematic::DoNeugenCalc(NueRecord *record,
				     Systematic::Systematic_t sys,Double_t val)
{

  MCReweight *mcr = &MCReweight::Instance();
//  mcr->ResetAllReweightConfigs();  

  static double oldVal = -999;
  bool reset = false;
  Registry rwtconfig;

  if(val != oldVal){
    rwtconfig.UnLockValues();
    rwtconfig.UnLockKeys();
    //rwtconfig.Set("neugen:config_name","MODBYRS");
    //rwtconfig.Set("neugen:config_no",3);
    rwtconfig.Set("neugen:use_scale_factors",1);
    if(sys==Systematic::kMA_QE) rwtconfig.Set("neugen:ma_qe",val);
    else if(sys==Systematic::kMA_RES) rwtconfig.Set("neugen:ma_res",val);
    else if(sys==Systematic::kKNO) rwtconfig.Set("neugen:scale_kno_all",val);
    else if(sys==Systematic::kShwDev) {
      rwtconfig.Set("neugen:config_no",int(val+0.5));
      rwtconfig.Set("neugen:config_name","MODBYRS");
    }
    rwtconfig.LockValues();
    rwtconfig.LockKeys();
    reset = true; 
    oldVal = val;
  }

  MCEventInfo ei;
  ei.UseStoredXSec(true);
  ei.nuE           = record->mctrue.nuEnergy;
  ei.nuPx          = record->mctrue.nuDCosX*ei.nuE;
  ei.nuPy          = record->mctrue.nuDCosY*ei.nuE;
  ei.nuPz          = record->mctrue.nuDCosZ*ei.nuE;
  ei.tarE          = record->mctrue.targetEnergy;
  ei.tarPx         = record->mctrue.targetPX;
  ei.tarPx         = record->mctrue.targetPY;
  ei.tarPx         = record->mctrue.targetPZ;  
  ei.y             = record->mctrue.hadronicY;
  ei.x             = record->mctrue.bjorkenX;
  ei.q2            = record->mctrue.q2;
  ei.w2            = record->mctrue.w2;  
  ei.iaction       = record->mctrue.interactionType;
  ei.inu           = record->mctrue.nuFlavor;
  ei.iresonance    = record->mctrue.resonanceCode;
  ei.initial_state = record->mctrue.initialState;
  ei.nucleus       = ReweightHelpers::
    FindNucleusNumber(int(record->mctrue.atomicNumber),
		      int(record->mctrue.atomicWeight));
  ei.had_fs        = abs(record->mctrue.hadronicFinalState);  
  ei.stdXSec       = record->xsecweights.xsecweight;

  //std::cout << ei.nuE << " " << ei.tarE << " " << ei.y << " " << ei.q2 << " " 
  //    << ei. inu << " " << ei.initial_state << " " << ei.nucleus 
  //    << " " << ei.had_fs << std::endl;
  if(ei.inu<-9998 || 
     (ei.iresonance==1003 && TMath::Abs(ei.had_fs)<200) ) return 1;
  NuParent *nuparent = 0;  

  double weight = 1.0;

  if(reset) weight = mcr->ComputeWeight(&ei,nuparent,&rwtconfig);
  else      weight = mcr->ComputeWeight(&ei,nuparent);

  return weight;
}

Double_t NueSystematic::DoCalibShift(NueRecord *record,
				     Systematic::Systematic_t sys,Double_t val) 
{
/*
  This is Chris's old code, but I don't quite know why he chose this params and
setup
    so I'm going a different way....  (JAB) 10/12/2007
                                                                                
  Double_t epi_ratio = 1.26; //use a rough e/pi ratio to adjust reco energy
  Double_t true_shw_energy = record->mctrue.trueVisibleE;
  if(true_shw_energy == 0 ) return 0;
  Double_t true_em_energy  = record->mctrue.emShowerFraction*true_shw_energy;
  Double_t respcor_true_energy = (true_shw_energy + (epi_ratio-1)*true_em_energy);
                                                                                
  Double_t true_elec_energy = 0;
  if(TMath::Abs(record->mctrue.nuFlavor)==12)
    true_elec_energy = TMath::Abs(record->mctrue.leptonMomentum);
                                                                                
  //Double_t respcor_em_frac = (epi_ratio*true_em_energy) / respcor_true_energy;
   Double_t respcor_em_frac = (epi_ratio*true_elec_energy) / respcor_true_energy;                                                                                
  Double_t cor_em_frac = respcor_em_frac;
*/
                                                                                
  Double_t cor_em_frac = 0.0;
                                                                                
  if(record->mctrue.fNueClass != 0){
    double EMen = record->shi.emenergy - record->shi.epi0;
    cor_em_frac = EMen/(record->mctrue.nuEnergy);
  }
                                                                                
  if(sys==Systematic::kEMCalib)
    return cor_em_frac*val;
  else if(sys==Systematic::kHadCalib)
    return (1-cor_em_frac)*val;
  else if(sys==Systematic::kRelCalib) {
    if(record->GetHeader().GetVldContext().GetDetector()!=2) return 0;
    return val;
  }
  return 0;
}

Double_t NueSystematic::DoTauProd(NueRecord *record,Double_t val)
{
  //Correspondance with Hugh G. on 1/23/2008
  //use 50% for QEl, RES, but 10% for DIS
  double scale = 0.5;

  if(record->mctrue.interactionType==1 &&
     TMath::Abs(record->mctrue.nuFlavor)==16){

     if(record->mctrue.resonanceCode == 1001) scale = 0.5;
     if(record->mctrue.resonanceCode == 1002) scale = 0.5;
     if(record->mctrue.resonanceCode == 1003) scale = 0.1;
     return 1 + scale*val;
  }
  return 1;
}

Double_t NueSystematic::DoPIDSkew(NueRecord *record,Double_t val,
				  Selection::Selection_t sel)
{
  double pid = NueStandard::GetPIDValue(record, sel);

  if(pid > -1000)
  {
    return (1 + pid*val);
  }

  return 1;
}

Double_t NueSystematic::DoNormCalc(NueRecord *record,Double_t val)
{
  if(record->GetHeader().GetVldContext().GetDetector()==2) return val;
  return 1;
}

Double_t NueSystematic::DoNCScale(NueRecord *record, 
                                  Double_t val,
                                  Selection::Selection_t sel)
{
  if(sel != Selection::kCC) return 1;
  if(record->mctrue.fNueClass == 0) return val;
  return 1;
}

Double_t NueSystematic::DoCCShwEnergyScale(NueRecord * /* record */,
                                           Double_t val,
                                           Selection::Selection_t sel)
{
  if(sel != Selection::kCC) return 0;
  return val;
}

void NueSystematic::MakeBranches(TTree *tree) 
{  
  TBranch *branch = 0;
  if((branch = tree->GetBranch("SysName"))) branch->SetAddress(fName);
  else tree->Branch("SysName",fName,"SysName/C");

  std::map<Systematic::Systematic_t,Double_t>::iterator beg = fSysList.begin();
  std::map<Systematic::Systematic_t,Double_t>::iterator end = fSysList.end();
  Int_t max_sys_index = 0;
  while(strcmp(Systematic::AsString(Systematic::ESystematic(max_sys_index)),
	       "?Unknown?")!=0) {
    //std::cout << max_sys_index << " " 
    //      << Systematic::AsString(Systematic::ESystematic(max_sys_index))
    //      << std::endl;
    if (beg!=end && beg->first==max_sys_index) {
      //first check if branch already exists:
      if((branch = tree->GetBranch(Systematic::AsString(beg->first)))) {
	branch->SetAddress(&(beg->second));
      }
      else {
	char leafname[256]; sprintf(leafname,"%s/D",Systematic::AsString(beg->first));
	tree->Branch(Systematic::AsString(beg->first),&(beg->second),leafname);
      }
      beg++;
    }
    else {
      if((branch = 
	  tree->GetBranch(Systematic::AsString(Systematic::ESystematic(max_sys_index))))) {
	branch->SetAddress(&fTempDouble);
      }
      else {
	char leafname[256];
	sprintf(leafname,"%s/D",
		Systematic::AsString(Systematic::ESystematic(max_sys_index))); 
	tree->Branch(Systematic::AsString(Systematic::ESystematic(max_sys_index)),
		     &fTempDouble,leafname);
      }
    }
    max_sys_index++;
  }

  if((branch = tree->GetBranch("Theta12"))) branch->SetAddress(&fTheta12);
  else tree->Branch("Theta12",&fTheta12,"Theta12/D");
  if((branch = tree->GetBranch("Theta23"))) branch->SetAddress(&fTheta23);
  else tree->Branch("Theta23",&fTheta23,"Theta23/D");
  if((branch = tree->GetBranch("Theta13"))) branch->SetAddress(&fTheta13);
  else tree->Branch("Theta13",&fTheta13,"Theta13/D");
  if((branch = tree->GetBranch("DeltaMSq23"))) branch->SetAddress(&fDeltaMSq23);
  else tree->Branch("DeltaMSq23",&fDeltaMSq23,"DeltaMSq23/D");
  if((branch = tree->GetBranch("DeltaMSq12"))) branch->SetAddress(&fDeltaMSq12);
  else tree->Branch("DeltaMSq12",&fDeltaMSq12,"DeltaMSq12/D");
  if((branch = tree->GetBranch("DeltaCP"))) branch->SetAddress(&fDeltaCP);
  else tree->Branch("DeltaCP",&fDeltaCP,"DeltaCP/D");
  if((branch = tree->GetBranch("MassHierarchy"))) branch->SetAddress(&fMassHierarchy);
  else tree->Branch("MassHierarchy",&fMassHierarchy,"MassHierarchy/I");
}
