#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "CandNtupleSR/NtpSRShieldStrip.h"
#include "CandNtupleSR/NtpSRShieldSummary.h"
#include "CandShield/CandShieldSR.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/ShieldRejVarAna.h"
#include "TClonesArray.h"
#include "TPad.h"
#include "TLatex.h"
#include "TNtuple.h"
#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"

using namespace std;


ShieldRejVarAna::ShieldRejVarAna(ShieldRejVar &sv):
  fShieldRejVar(sv)
{
}

ShieldRejVarAna::~ShieldRejVarAna()
{}

void ShieldRejVarAna::Analyze(int evtn,  RecRecordImp<RecCandHeader> *srobj)
{
  if(srobj==0){
    return;
  }
  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }

  fShieldRejVar.Reset();

  Detector::Detector_t det = srobj->GetHeader().GetVldContext().GetDetector();
  if(det != Detector::kFar) return;

  if(srobj->GetHeader().GetVldContext().GetSimFlag()==SimFlag::kMC) return;


  const NtpStRecord* record = dynamic_cast<const NtpStRecord *>(srobj);

  ShieldGeom *shge = new ShieldGeom((srobj->GetHeader()).GetVldContext());

  NtpSREvent *event = SntpHelpers::GetEvent(evtn,srobj);

  Double_t vtxTime=0;
  //Shield 
  Int_t vetoHitPlane=0;
  Int_t vetoHitStrip0=0;
  Double_t vetoTime[2]={0};
  Int_t vetoHitSec=-1;
  Int_t contFound=0;
  
  Float_t vtxZ = event->vtx.z;
  Float_t vtxT = event->vtx.t;
                                                                                
  if(ReleaseType::IsCedar(release)){
    NtpStRecord* st = dynamic_cast<NtpStRecord *>(srobj);
    NtpVtxFinder vtxf(evtn, st);
    if(vtxf.FindVertex() > 0){
       vtxZ = vtxf.VtxZ();
       vtxT = vtxf.VtxT();
    }
  }

  //Checking if hits in shield are found
  //------------------------------------
  //  NtpSRShieldSummary *shieldSummary = SntpHelpers::shieldSummary(srobj);
  const long veto_stp = record->vetostp->GetEntries(); 
  vtxTime = 1.e9*vtxT;	
  contFound=0;

  Int_t secabove0 = shge->ClosestTwoSections(vtxZ,0);
  Int_t secabove1 = shge->ClosestTwoSections(vtxZ,1);

  for(int idig = 0; idig<veto_stp; ++idig){
    NtpSRShieldStrip* ntpShieldStrip = SntpHelpers::GetShieldStrip(idig,srobj);
    if(ntpShieldStrip == NULL) {
      continue;
    }    
    vetoHitPlane = ntpShieldStrip->pln;
    vetoHitStrip0 = ntpShieldStrip->plank;
    vetoTime[0]=1.e9*(ntpShieldStrip->time[0]);
    vetoTime[1]=1.e9*(ntpShieldStrip->time[1]);
    vetoHitSec = shge->WhatSection(vetoHitPlane);

    if(((vetoHitSec==secabove0&&secabove0>0)||(vetoHitSec==secabove1&&secabove1>0))&& (fabs(vetoTime[0]-vtxTime)<100||fabs(vetoTime[1]-vtxTime)<100)){//same section && in time
      contFound+=1;
    }
  }
  fShieldRejVar.ShieldHit = contFound;

  if(shge){
   delete shge;
  }
}
  
