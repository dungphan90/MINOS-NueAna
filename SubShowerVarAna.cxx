#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/SubShowerVarAna.h"
#include "NueAna/SubShowerNN.h"
#include "TPad.h"
#include "TLatex.h"
#include "TNtuple.h"
#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>

using namespace std;


const Float_t STP_WIDTH = 0.0412; //should be the same for the near and far detectors


SubShowerVarAna::SubShowerVarAna(SubShowerVar &sv):
  fSubShowerVar(sv)
{
}

SubShowerVarAna::~SubShowerVarAna()
{}

void SubShowerVarAna::Analyze(int evtn,  RecRecordImp<RecCandHeader> *srobj)
{
  if(srobj==0){
    return;
  }
  
  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }
  
  NtpSRShower *shower = SntpHelpers::GetPrimaryShower(evtn,srobj);
  Analyze(srobj,shower);
  
}

void SubShowerVarAna::Analyze(RecRecordImp<RecCandHeader> *srobj, NtpSRShower* ntpShower=0)
{    
  
  fSubShowerVar.Reset();
  
  if(srobj==0){
    return;
  }
  
  if (ntpShower==0) {
    return;
  }
  
  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }
  
  fSubShowerVar.Zero();
  
  fSubShowerVar.ncluster = ntpShower->ncluster;
  fSubShowerVar.nclusterU = ntpShower->nUcluster;
  fSubShowerVar.nclusterV = ntpShower->nVcluster;
  
  Double_t totphU = 0.;
  Double_t totphV = 0.;
  Double_t E0U = 0.;
  Double_t E0V = 0.;
  Double_t maxU = 0.;
  Double_t maxV = 0.;
  Double_t phystotphU = 0.;
  Double_t phystotphV = 0.;
  
  for(int k=0; k<fSubShowerVar.ncluster; k++){
    Int_t index = ntpShower->clu[k];
    NtpSRCluster *ntpCluster = SntpHelpers::GetCluster(index,srobj);
    
    if(ntpCluster->planeview==2){
      
      totphU += ntpCluster->ph.gev;
      fSubShowerVar.PHAvgIDU += (ntpCluster->id)*ntpCluster->ph.gev;
      fSubShowerVar.PHAvgProbEMU += (ntpCluster->probem)*ntpCluster->ph.gev;
      fSubShowerVar.PHAvgDevU += (ntpCluster->avgdev)*ntpCluster->ph.gev;
      
      if(ntpCluster->ph.gev>0.5) fSubShowerVar.nclusterU_he += 1;
      
      if(ntpCluster->id==0||ntpCluster->id==1||ntpCluster->id==3) {
	fSubShowerVar.nPhysClusterU+=1;
	fSubShowerVar.PHFracRMSU += (ntpCluster->ph.gev*ntpCluster->ph.gev);
	phystotphU += ntpCluster->ph.gev;
      }
      if(ntpCluster->ph.gev>maxU){
	maxU = ntpCluster->ph.gev;
	fSubShowerVar.nstp0U = ntpCluster->nstrip;
	E0U = maxU;
      }
      
      if(ntpCluster->id==0) {
	fSubShowerVar.nEMU+=1;
	if(ntpCluster->ph.gev>0.5) fSubShowerVar.nEMU_he+=1;
	fSubShowerVar.eEMU+=ntpCluster->ph.gev;
      }
      else if(ntpCluster->id==1) {
	fSubShowerVar.nHadU+=1;
	if(ntpCluster->ph.gev>0.5) fSubShowerVar.nHadU_he+=1;
	fSubShowerVar.eHadU+=ntpCluster->ph.gev;
      }
      else if(ntpCluster->id==3) {
	fSubShowerVar.nTrkU+=1;
	if(ntpCluster->ph.gev>0.5) fSubShowerVar.nTrkU_he+=1;
	fSubShowerVar.eTrkU+=ntpCluster->ph.gev;
      }
      else if(ntpCluster->id==2) {
	fSubShowerVar.nXTkU+=1;
	if(ntpCluster->ph.gev>0.5) fSubShowerVar.nXTkU_he+=1;
	fSubShowerVar.eXTkU+=ntpCluster->ph.gev;
      }
      else if(ntpCluster->id==4 || ntpCluster->id==5) {
	fSubShowerVar.nHalU+=1;
	if(ntpCluster->ph.gev>0.5) fSubShowerVar.nHalU_he+=1;
	fSubShowerVar.eHalU+=ntpCluster->ph.gev;
      }
      
    }
    else if(ntpCluster->planeview==3){
      
      totphV += ntpCluster->ph.gev;
      fSubShowerVar.PHAvgIDV += (ntpCluster->id)*ntpCluster->ph.gev;
      fSubShowerVar.PHAvgProbEMV += (ntpCluster->probem)*ntpCluster->ph.gev;
      fSubShowerVar.PHAvgDevV += (ntpCluster->avgdev)*ntpCluster->ph.gev;
      
      if(ntpCluster->ph.gev>0.5) fSubShowerVar.nclusterV_he += 1;
      
      if(ntpCluster->id==0||ntpCluster->id==1||ntpCluster->id==3) {
	fSubShowerVar.nPhysClusterV+=1;
	fSubShowerVar.PHFracRMSV += (ntpCluster->ph.gev*ntpCluster->ph.gev);
	phystotphV += ntpCluster->ph.gev;
      }
      if(ntpCluster->ph.gev>maxV){
	maxV = ntpCluster->ph.gev;
	fSubShowerVar.nstp0V = ntpCluster->nstrip;
	E0V = maxV;
      }
      
      if(ntpCluster->id==0) {
	fSubShowerVar.nEMV+=1;
	if(ntpCluster->ph.gev>0.5) fSubShowerVar.nEMV_he+=1;
	fSubShowerVar.eEMV+=ntpCluster->ph.gev;
      }
      else if(ntpCluster->id==1) {
	fSubShowerVar.nHadV+=1;
	if(ntpCluster->ph.gev>0.5) fSubShowerVar.nHadV_he+=1;
	fSubShowerVar.eHadV+=ntpCluster->ph.gev;
      }
      else if(ntpCluster->id==3) {
	fSubShowerVar.nTrkV+=1;
	if(ntpCluster->ph.gev>0.5) fSubShowerVar.nTrkV_he+=1;
	fSubShowerVar.eTrkV+=ntpCluster->ph.gev;
      }
      else if(ntpCluster->id==2) {
	fSubShowerVar.nXTkV+=1;
	if(ntpCluster->ph.gev>0.5) fSubShowerVar.nXTkV_he+=1;
	fSubShowerVar.eXTkV+=ntpCluster->ph.gev;
      }
      else if(ntpCluster->id==4 || ntpCluster->id==5) {
	fSubShowerVar.nHalV+=1;
	if(ntpCluster->ph.gev>0.5) fSubShowerVar.nHalV_he+=1;
	fSubShowerVar.eHalV+=ntpCluster->ph.gev;
      }    
    }
  }
  if(totphU>0){
    fSubShowerVar.PHAvgIDU /= totphU;
    fSubShowerVar.PHAvgProbEMU /= totphU;
    fSubShowerVar.PHAvgDevU /= totphU;
  }
  else{
    fSubShowerVar.PHAvgIDU = ANtpDefVal::kDouble;
    fSubShowerVar.PHAvgProbEMU = ANtpDefVal::kDouble;
    fSubShowerVar.PHAvgDevU = ANtpDefVal::kDouble;
  }
  if(totphV>0){
    fSubShowerVar.PHAvgIDV /= totphV;
    fSubShowerVar.PHAvgProbEMV /= totphV;
    fSubShowerVar.PHAvgDevV /= totphV;
  }
  else{
    fSubShowerVar.PHAvgIDV = ANtpDefVal::kDouble;
    fSubShowerVar.PHAvgProbEMV = ANtpDefVal::kDouble;
    fSubShowerVar.PHAvgDevV = ANtpDefVal::kDouble;
  }
  
  Double_t maxU1 = 0.;
  Double_t maxV1 = 0.;
  for(int k=0; k<fSubShowerVar.ncluster; k++){
    Int_t index = ntpShower->clu[k];
    NtpSRCluster *ntpCluster = SntpHelpers::GetCluster(index,srobj);
    if(ntpCluster->planeview==2 && E0U>0){
      if(ntpCluster->ph.gev>maxU1&&ntpCluster->ph.gev<E0U){
	fSubShowerVar.E2to1U = ntpCluster->ph.gev/E0U;
	maxU1 = ntpCluster->ph.gev;
      }
    }
    else if(ntpCluster->planeview==3 && E0V>0){
      if(ntpCluster->ph.gev>maxV1&&ntpCluster->ph.gev<E0V){
	fSubShowerVar.E2to1V = ntpCluster->ph.gev/E0V;
	maxV1 = ntpCluster->ph.gev;
      }
    }
  }
  
  if(ntpShower->nUcluster>0&&phystotphU>0){
    fSubShowerVar.PHFracRMSU /= (phystotphU*phystotphU*fSubShowerVar.nclusterU);
    fSubShowerVar.PHFracRMSU -= TMath::Power(1./fSubShowerVar.nclusterU,2);
    if(fSubShowerVar.PHFracRMSU>0) 
      fSubShowerVar.PHFracRMSU  = TMath::Sqrt(fSubShowerVar.PHFracRMSU);
    else fSubShowerVar.PHFracRMSU = 0;
  }
  else{
    fSubShowerVar.PHFracRMSU = ANtpDefVal::kDouble;
  }
  
  if(ntpShower->nVcluster>0&&phystotphV>0){
    fSubShowerVar.PHFracRMSV /= (phystotphV*phystotphV*fSubShowerVar.nclusterV);
    fSubShowerVar.PHFracRMSV -= TMath::Power(1./fSubShowerVar.nclusterV,2);
    if(fSubShowerVar.PHFracRMSV>0) 
      fSubShowerVar.PHFracRMSV  = TMath::Sqrt(fSubShowerVar.PHFracRMSV);
    else fSubShowerVar.PHFracRMSV = 0;
  }
  else{
    fSubShowerVar.PHFracRMSV = ANtpDefVal::kDouble;
  }
  
  //calculate SS pid...
  if(phystotphU>0&&phystotphV>0){
    SubShowerNN snn;
    Double_t params[4];
    params[0] = fSubShowerVar.PHFracRMSU + fSubShowerVar.PHFracRMSV;
    params[1] = (fSubShowerVar.PHAvgProbEMU*fSubShowerVar.PHAvgProbEMU + 
		 fSubShowerVar.PHAvgProbEMV*fSubShowerVar.PHAvgProbEMV)/2.;
    params[2] = (fSubShowerVar.PHAvgIDU + fSubShowerVar.PHAvgIDV)/2.;
    params[3] = fSubShowerVar.PHAvgDevU + fSubShowerVar.PHAvgDevV;
    fSubShowerVar.pid = snn.Value(0,params[0],params[1],params[2],params[3]);
  }
  else { 
   fSubShowerVar.pid = ANtpDefVal::kDouble;
  }

  //effective number:
  if(maxU>0) {
    fSubShowerVar.wEMU  = fSubShowerVar.eEMU/maxU;
    fSubShowerVar.wHadU = fSubShowerVar.eHadU/maxU;
    fSubShowerVar.wTrkU = fSubShowerVar.eTrkU/maxU;
    fSubShowerVar.wXTkU = fSubShowerVar.eXTkU/maxU;
    fSubShowerVar.wHalU = fSubShowerVar.eHalU/maxU;
  }
  if(maxV>0) {
    fSubShowerVar.wEMV  = fSubShowerVar.eEMV/maxV;
    fSubShowerVar.wHadV = fSubShowerVar.eHadV/maxV;
    fSubShowerVar.wTrkV = fSubShowerVar.eTrkV/maxV;
    fSubShowerVar.wXTkV = fSubShowerVar.eXTkV/maxV;
    fSubShowerVar.wHalV = fSubShowerVar.eHalV/maxV;
  }
}

			  
