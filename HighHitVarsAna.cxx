#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"
#include "HighHitVars.h"
#include "HighHitVarsAna.h"
#include "MessageService/MsgService.h"
#include "TPad.h"
#include "TLatex.h"
#include "TNtuple.h"
#include "TFile.h"
#include "TMath.h"
#include <cassert>
#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;


HighHitVarsAna::HighHitVarsAna(HighHitVars &hhv):
  fHighHitVars(hhv)
{
}


HighHitVarsAna::~HighHitVarsAna()
{}

void HighHitVarsAna::Analyze(int evtn,  RecRecordImp<RecCandHeader> *srobj)
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

void HighHitVarsAna::Analyze(RecRecordImp<RecCandHeader> *srobj, NtpSRShower* ntpShower=0)
{  

 fHighHitVars.Reset();

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
  

  float A1, A2, A3, A4, A5, A6;
  float A0, B0, C0, D0, E0, F0;

  A0=0;
  B0=0;
  C0=0;
  D0=0;
  E0=0;
  F0=0;

  int showerstrips;
  showerstrips=ntpShower->nstrip;

  int j, k, l, m, n;
  j = k = l = m = n = 0;

  for(int i=0;i<showerstrips;i++){
    Int_t index = ntpShower->stp[i]; 
     NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
      if(!strip){
	 continue;
      }

     ///first hit
      A1 = evtstp0mip[index] + evtstp1mip[index];

      if(A1>A0){
	fHighHitVars.high_hit_1_ph = strip->ph0.sigcor+strip->ph1.sigcor;
        fHighHitVars.high_hit_1_mip = A1;
        fHighHitVars.high_hit_1_strip = strip->strip;
        fHighHitVars.high_hit_1_plane = strip->plane;
        fHighHitVars.high_hit_1_planeview = strip->planeview;
        fHighHitVars.high_hit_1_zpos = strip->z;
        fHighHitVars.high_hit_1_tpos = strip->tpos;
        j=i;
      }

      if(A1<=A0){
        continue;
      }

    A0 = fHighHitVars.high_hit_1_mip;

  }

  ///1st hit loop finished


  Float_t hit1;
  hit1=fHighHitVars.high_hit_1_mip;

  for(int i=0;i<showerstrips;i++){
      Int_t index = ntpShower->stp[i];
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
      if(!strip){
	 continue;
      }

  if(showerstrips>1){
  
     A2 = evtstp0mip[index] + evtstp1mip[index];

     if(A2<=hit1&&A2>B0&&i!=j){
       	fHighHitVars.high_hit_2_ph = strip->ph0.sigcor+strip->ph1.sigcor;
        fHighHitVars.high_hit_2_mip = A2;
        fHighHitVars.high_hit_2_strip = strip->strip;
        fHighHitVars.high_hit_2_plane = strip->plane;
        fHighHitVars.high_hit_2_planeview = strip->planeview;
        fHighHitVars.high_hit_2_zpos = strip->z;
        fHighHitVars.high_hit_2_tpos = strip->tpos;
        k=i;
     }
     if(A2<B0){
       continue;
     }

    B0 = fHighHitVars.high_hit_2_mip;
  }

  if(showerstrips<=1){
       	fHighHitVars.high_hit_2_ph = ANtpDefVal::kFloat;
       	fHighHitVars.high_hit_2_mip = ANtpDefVal::kFloat;
        fHighHitVars.high_hit_2_strip = ANtpDefVal::kInt;
        fHighHitVars.high_hit_2_plane = ANtpDefVal::kInt;
        fHighHitVars.high_hit_2_planeview = 7;
        fHighHitVars.high_hit_2_zpos = ANtpDefVal::kFloat;
        fHighHitVars.high_hit_2_tpos = ANtpDefVal::kFloat;
  }
  }
   ///2nd hit loop finished 


  Float_t hit2;
  hit2=fHighHitVars.high_hit_2_mip;


  for(int i=0;i<showerstrips;i++){
      Int_t index = ntpShower->stp[i];
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
      if(!strip){
	 continue;
      }


  if(showerstrips>2){

     A3 = evtstp0mip[index] + evtstp1mip[index];

     if(A3<=hit2&&A3>C0&&i!=j&&i!=k){
       	fHighHitVars.high_hit_3_ph = strip->ph0.sigcor+strip->ph1.sigcor;
        fHighHitVars.high_hit_3_mip = A3;
        fHighHitVars.high_hit_3_strip = strip->strip;
        fHighHitVars.high_hit_3_plane = strip->plane;
        fHighHitVars.high_hit_3_planeview = strip->planeview;
        fHighHitVars.high_hit_3_zpos = strip->z;
        fHighHitVars.high_hit_3_tpos = strip->tpos;
        l=i;
     }
     if(A3<C0){
       continue;
     }

    C0 = fHighHitVars.high_hit_3_mip;
  }

  if(showerstrips<=2){
       	fHighHitVars.high_hit_3_ph = ANtpDefVal::kFloat;
       	fHighHitVars.high_hit_3_mip = ANtpDefVal::kFloat;
        fHighHitVars.high_hit_3_strip = ANtpDefVal::kInt;
        fHighHitVars.high_hit_3_plane = ANtpDefVal::kInt;
        fHighHitVars.high_hit_3_planeview = 7;
        fHighHitVars.high_hit_3_zpos = ANtpDefVal::kFloat;
        fHighHitVars.high_hit_3_tpos = ANtpDefVal::kFloat;
  }
  }
   ///3rd hit loop finished 


  Float_t hit3;
  hit3=fHighHitVars.high_hit_3_mip;

  for(int i=0;i<showerstrips;i++){
      Int_t index = ntpShower->stp[i];
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
      if(!strip){
	 continue;
      }


  if(showerstrips>3){

     A4 = evtstp0mip[index] + evtstp1mip[index];

     if(A4<=hit3&&A4>D0&&i!=j&&i!=k&&i!=l){
       	fHighHitVars.high_hit_4_ph = strip->ph0.sigcor+strip->ph1.sigcor;
        fHighHitVars.high_hit_4_mip = A4;
        fHighHitVars.high_hit_4_strip = strip->strip;
        fHighHitVars.high_hit_4_plane = strip->plane;
        fHighHitVars.high_hit_4_planeview = strip->planeview;
        fHighHitVars.high_hit_4_zpos = strip->z;
        fHighHitVars.high_hit_4_tpos = strip->tpos;
        m=i;
     }
     if(A4<D0){
       continue;
     }

    D0 = fHighHitVars.high_hit_4_mip;
  }

  if(showerstrips<=3){
       	fHighHitVars.high_hit_4_ph = ANtpDefVal::kFloat;
       	fHighHitVars.high_hit_4_mip = ANtpDefVal::kFloat;
        fHighHitVars.high_hit_4_strip = ANtpDefVal::kInt;
        fHighHitVars.high_hit_4_plane = ANtpDefVal::kInt;
        fHighHitVars.high_hit_4_planeview = 7;
        fHighHitVars.high_hit_4_zpos = ANtpDefVal::kFloat;
        fHighHitVars.high_hit_4_tpos = ANtpDefVal::kFloat;
  }
  }
   ///4th hit loop finished 

  Float_t hit4;
  hit4=fHighHitVars.high_hit_4_mip;


  for(int i=0;i<showerstrips;i++){
      Int_t index = ntpShower->stp[i];
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
      if(!strip){
	 continue;
      }

  if(showerstrips>4){

     A5 = evtstp0mip[index] + evtstp1mip[index];

     if(A5<=hit4&&A5>E0&&i!=j&&i!=k&&i!=l&&i!=m){
       	fHighHitVars.high_hit_5_ph = strip->ph0.sigcor+strip->ph1.sigcor;
        fHighHitVars.high_hit_5_mip = A5;
        fHighHitVars.high_hit_5_strip = strip->strip;
        fHighHitVars.high_hit_5_plane = strip->plane;
        fHighHitVars.high_hit_5_planeview = strip->planeview;
        fHighHitVars.high_hit_5_zpos = strip->z;
        fHighHitVars.high_hit_5_tpos = strip->tpos;
        n=i;
     }
     if(A5<E0){
       continue;
     }

    E0 = fHighHitVars.high_hit_5_mip;
  }

  if(showerstrips<=4){
       	fHighHitVars.high_hit_5_ph = ANtpDefVal::kFloat;
       	fHighHitVars.high_hit_5_mip = ANtpDefVal::kFloat;
        fHighHitVars.high_hit_5_strip = ANtpDefVal::kInt;
        fHighHitVars.high_hit_5_plane = ANtpDefVal::kInt;
        fHighHitVars.high_hit_5_planeview = 7;
        fHighHitVars.high_hit_5_zpos = ANtpDefVal::kFloat;
        fHighHitVars.high_hit_5_tpos = ANtpDefVal::kFloat;
  }
  }
   ///5th hit loop finished 

  Float_t hit5;
  hit5=fHighHitVars.high_hit_5_mip;

  for(int i=0;i<showerstrips;i++){
      Int_t index = ntpShower->stp[i];
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
      if(!strip){
	 continue;
      }

  if(showerstrips>5){

     A6 = evtstp0mip[index] + evtstp1mip[index];

     if(A6<=hit5&&A6>F0&&i!=j&&i!=k&&i!=l&&i!=m&&i!=n){
       	fHighHitVars.high_hit_6_ph = strip->ph0.sigcor+strip->ph1.sigcor;
        fHighHitVars.high_hit_6_mip = A6;
        fHighHitVars.high_hit_6_strip = strip->strip;
        fHighHitVars.high_hit_6_plane = strip->plane;
        fHighHitVars.high_hit_6_planeview = strip->planeview;
        fHighHitVars.high_hit_6_zpos = strip->z;
        fHighHitVars.high_hit_6_tpos = strip->tpos;
     }
     if(A6<F0){
       continue;
     }

    F0 = fHighHitVars.high_hit_6_mip;
  }

  if(showerstrips<=5){
       	fHighHitVars.high_hit_6_ph = ANtpDefVal::kFloat;
       	fHighHitVars.high_hit_6_mip = ANtpDefVal::kFloat;
        fHighHitVars.high_hit_6_strip = ANtpDefVal::kInt;
        fHighHitVars.high_hit_6_plane = ANtpDefVal::kInt;
        fHighHitVars.high_hit_6_planeview = 7;
        fHighHitVars.high_hit_6_zpos = ANtpDefVal::kFloat;
        fHighHitVars.high_hit_6_tpos = ANtpDefVal::kFloat;
  }
  }
   ///6th hit loop finished 

  //Variables added at Jenny's request - find all the strips above a certain threshold and count them or total their ph:                                                  

 float total_ph_1mip = 0;
 float total_ph_2mip = 0;
 float total_ph_3mip = 0;
 float total_ph_4mip = 0;
 float total_ph_5mip = 0;
 float total_ph_6mip = 0;

 float total_mip_1mip = 0;
 float total_mip_2mip = 0;
 float total_mip_3mip = 0;
 float total_mip_4mip = 0;
 float total_mip_5mip = 0;
 float total_mip_6mip = 0;

 int total_strips_1mip = 0;
 int total_strips_2mip = 0;
 int total_strips_3mip = 0;
 int total_strips_4mip = 0;
 int total_strips_5mip = 0;
 int total_strips_6mip = 0;

  float tempMIP=0;

  for(int i=0;i<showerstrips;i++){
    Int_t index = ntpShower->stp[i]; 
     NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
      if(!strip){
	 continue;
      }
 
    tempMIP = evtstp0mip[index] + evtstp1mip[index];

    if(tempMIP>1.0){
      total_ph_1mip=total_ph_1mip+strip->ph0.sigcor+strip->ph1.sigcor;
      total_mip_1mip=total_mip_1mip+tempMIP;
      total_strips_1mip=total_strips_1mip+1;
    }

    if(tempMIP>2.0){
      total_ph_2mip=total_ph_2mip+strip->ph0.sigcor+strip->ph1.sigcor;
      total_mip_2mip=total_mip_2mip+tempMIP;
      total_strips_2mip=total_strips_2mip+1;
    }
    if(tempMIP>3.0){
      total_ph_3mip=total_ph_3mip+strip->ph0.sigcor+strip->ph1.sigcor;
      total_mip_3mip=total_mip_3mip+tempMIP;
      total_strips_3mip=total_strips_3mip+1;
    }

    if(tempMIP>4.0){
      total_ph_4mip=total_ph_4mip+strip->ph0.sigcor+strip->ph1.sigcor;
      total_mip_4mip=total_mip_4mip+tempMIP;
      total_strips_4mip=total_strips_4mip+1;
    }

    if(tempMIP>5.0){
      total_ph_5mip=total_ph_5mip+strip->ph0.sigcor+strip->ph1.sigcor;
      total_mip_5mip=total_mip_5mip+tempMIP;
      total_strips_5mip=total_strips_5mip+1;
    }
    if(tempMIP>6.0){
      total_ph_6mip=total_ph_6mip+strip->ph0.sigcor+strip->ph1.sigcor;
      total_mip_6mip=total_mip_6mip+tempMIP;
      total_strips_6mip=total_strips_6mip+1;
    }
  
  }


  if(total_strips_1mip>0){
    fHighHitVars.hitsabove1MIP_total_ph = total_ph_1mip;
    fHighHitVars.hitsabove1MIP_total_mip = total_mip_1mip;
    fHighHitVars.hitsabove1MIP_total_strips = total_strips_1mip;
  }
  if(total_strips_1mip<=0){
    fHighHitVars.hitsabove1MIP_total_ph = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove1MIP_total_mip = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove1MIP_total_strips = ANtpDefVal::kInt;
  }

  if(total_strips_2mip>0){
    fHighHitVars.hitsabove2MIP_total_ph = total_ph_2mip;
    fHighHitVars.hitsabove2MIP_total_mip = total_mip_2mip;
    fHighHitVars.hitsabove2MIP_total_strips = total_strips_2mip;
  }
  if(total_strips_2mip<=0){
    fHighHitVars.hitsabove2MIP_total_ph = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove2MIP_total_mip = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove2MIP_total_strips = ANtpDefVal::kInt;
  }

  if(total_strips_3mip>0){
    fHighHitVars.hitsabove3MIP_total_ph = total_ph_3mip;
    fHighHitVars.hitsabove3MIP_total_mip = total_mip_3mip;
    fHighHitVars.hitsabove3MIP_total_strips = total_strips_3mip;
  }
  if(total_strips_3mip<=0){
    fHighHitVars.hitsabove3MIP_total_ph = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove3MIP_total_mip = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove3MIP_total_strips = ANtpDefVal::kInt;
  }

  if(total_strips_4mip>0){
    fHighHitVars.hitsabove4MIP_total_ph = total_ph_4mip;
    fHighHitVars.hitsabove4MIP_total_mip = total_mip_4mip;
    fHighHitVars.hitsabove4MIP_total_strips = total_strips_4mip;
  }
  if(total_strips_4mip<=0){
    fHighHitVars.hitsabove4MIP_total_ph = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove4MIP_total_mip = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove4MIP_total_strips = ANtpDefVal::kInt;
  }


  if(total_strips_5mip>0){
    fHighHitVars.hitsabove5MIP_total_ph = total_ph_5mip;
    fHighHitVars.hitsabove5MIP_total_mip = total_mip_5mip;
    fHighHitVars.hitsabove5MIP_total_strips = total_strips_5mip;
  }
  if(total_strips_5mip<=0){
    fHighHitVars.hitsabove5MIP_total_ph = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove5MIP_total_mip = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove5MIP_total_strips = ANtpDefVal::kInt;
  }


  if(total_strips_6mip>0){
    fHighHitVars.hitsabove6MIP_total_ph = total_ph_6mip;
    fHighHitVars.hitsabove6MIP_total_mip = total_mip_6mip;
    fHighHitVars.hitsabove6MIP_total_strips = total_strips_6mip;
  }
  if(total_strips_6mip<=0){
    fHighHitVars.hitsabove6MIP_total_ph = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove6MIP_total_mip = ANtpDefVal::kFloat;
    fHighHitVars.hitsabove6MIP_total_strips = ANtpDefVal::kInt;
  }



  //filling out the remaining vars
  fHighHitVars.showervtx_plane=ntpShower->vtx.plane;
  fHighHitVars.showervtx_u_pos=ntpShower->vtx.u;
  fHighHitVars.showervtx_v_pos=ntpShower->vtx.v;
  fHighHitVars.showerbeg_plane=ntpShower->plane.beg;

 


}

