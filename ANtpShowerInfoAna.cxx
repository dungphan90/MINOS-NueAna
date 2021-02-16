#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "MessageService/MsgService.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/ANtpShowerInfoAna.h"
#include "AnalysisNtuples/Module/ANtpRecoNtpManipulator.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "AtNuOutput/MOISolution.h"
#include "DataUtil/EnergyCorrections.h"
#include <cmath>

//CVSID("$Id: ANtpShowerInfoAna.cxx,v 1.36 2009/07/03 14:45:34 vahle Exp $");
using namespace EnergyCorrections;

ANtpShowerInfoAna::ANtpShowerInfoAna(ANtpShowerInfoNue &ansi):
   fANtpShowerInfo(ansi)
{}

ANtpShowerInfoAna::~ANtpShowerInfoAna()
{
}

//void ANtpShowerInfoAna::Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord * /*mcobj*/, NtpTHRecord * /*thobj*/)
void ANtpShowerInfoAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{    
  if(srobj==0){
    return;
  }
  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }

  NtpSREvent *event = 0;
  NtpSRShower *shower = 0;
  NtpSRTrack *track = 0;
  Detector::Detector_t  det;

  //and now an ugly bit of code to deal with either NtpStRecord or NtpSRRecord
  bool foundst=false;
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(srobj);
    fInfoFiller = new ANtpInfoObjectFiller();

  SimFlag::SimFlag_t s = SimFlag::kUnknown;

  if(ReleaseType::IsDogwood(release))
     SetCorrectionVersion(EnergyCorrections::kDogwood);
  if(ReleaseType::IsCedar(release))
     SetCorrectionVersion(EnergyCorrections::kCedar);
  if(ReleaseType::IsBirch(release))
     SetCorrectionVersion(EnergyCorrections::kBirch);

  if(st!=0){
    foundst=true;    
    //instansiate a NtpHelper object to help you get the info you want
    ANtpRecoNtpManipulator ntpManipulator(st);  
    fInfoFiller->SetStripArray(ntpManipulator.GetStripArray());

    event = SntpHelpers::GetEvent(evtn, st);
    track = SntpHelpers::GetPrimaryTrack(evtn, st);
    shower = SntpHelpers::GetPrimaryShower(evtn, st);

    det = st->GetHeader().GetVldContext().GetDetector();
    s = st->GetHeader().GetVldContext().GetSimFlag();
  }


  if(!foundst){  //Old Code, only good for R1.12 or earlier
      NtpSRRecord *sr=dynamic_cast<NtpSRRecord *>(srobj);    
      //instansiate a NtpHelper object to help you get the info you want
      ANtpRecoNtpManipulator ntpManipulator(sr);
      fInfoFiller->SetStripArray(ntpManipulator.GetStripArray());
	
      //set up which flags you want to use to determine the primary shower or track
      //a value of 0 for a flag means it will not be used
      ntpManipulator.SetPrimaryShowerCriteria(0,1); // nplanes, total pulse height
      //get the primary shower for the event - if no track is present it
      //returns 0     
    ANtpEventManipulator * ntpEventManipulator =
                           ntpManipulator.GetEventManipulator();
                                   
    ntpEventManipulator->SetEventInSnarl(evtn);
    event=ntpEventManipulator->GetEvent();
    shower = ntpEventManipulator->GetPrimaryShower();
  }  //end of old code


  if(shower) {
    //          Fill information from the base ANtpShowerInfo  
    fInfoFiller->FillShowerInformation(shower, &fANtpShowerInfo);
    //          Fill information specific to the ANtpShowerInfoNue class
    FillNueShowerInformation(shower, event, &fANtpShowerInfo,st);
    VldContext vc = st->GetHeader().GetVldContext();

    fANtpShowerInfo.phCCGeV = RecoShwEnergyNew(shower, 0, vc);
    fANtpShowerInfo.phNCGeV = shower->shwph.linNCgev;

    FillGapInformation(event, shower, track, st);

  }

  if(fInfoFiller){
    delete fInfoFiller;
    fInfoFiller=0;
  }

  
}

//----------------------------------------------------------------------
void ANtpShowerInfoAna::FillNueShowerInformation(NtpSRShower *ntpShower, NtpSREvent *ntpEvent, ANtpShowerInfoNue *showerInfoNue, RecRecordImp<RecCandHeader> *srobj)
{

    if(ntpEvent->nstrip){
        showerInfoNue->stripRatio=static_cast<Float_t>(ntpShower->nstrip)/static_cast<Float_t>(ntpEvent->nstrip); 
    }
    if(ntpEvent->plane.n){
        showerInfoNue->planeRatio=static_cast<Float_t>(ntpShower->plane.n)/static_cast<Float_t>(ntpEvent->plane.n);
    }
    if(ntpEvent->ph.mip){
        showerInfoNue->pulseHeightRatio=ntpShower->ph.mip/ntpEvent->ph.mip; 
    }
 
    showerInfoNue->phMeu = ntpShower->ph.mip;
    showerInfoNue->phMip = ntpShower->ph.mip;    
    showerInfoNue->phNueGeV = showerInfoNue->phMip/MeuPerGeV;    
    
    //borrowing moment of intertia calculation from AtNuOutput package
    vector<double> UStripsU, UStripsZ, UStripsE;
    vector<double> VStripsU, VStripsZ, VStripsE;
    for(int i=0;i<ntpShower->nstrip;i++){
      NtpSRStrip *strip = SntpHelpers::GetStrip(ntpShower->stp[i],srobj);
      float energy = ntpShower->stpph1mip[i] + ntpShower->stpph0mip[i];

      if(!strip) continue;
      if(strip->planeview==2){
	UStripsU.push_back(strip->tpos);
	UStripsZ.push_back(strip->z);
	UStripsE.push_back(energy);
      }
      else if(strip->planeview==3){
	VStripsU.push_back(strip->tpos);
	VStripsZ.push_back(strip->z);
	VStripsE.push_back(energy);
      }
    }

    MOISolution MOIUZ(UStripsU,UStripsZ,UStripsE);
    showerInfoNue->EValUZ0 = MOIUZ.EigenValues[0];
    showerInfoNue->EValUZ1 = MOIUZ.EigenValues[1];
    showerInfoNue->EVecUZ0[0] = MOIUZ.EigenVectors[0][0];
    showerInfoNue->EVecUZ1[0] = MOIUZ.EigenVectors[1][0];
    showerInfoNue->EVecUZ0[1] = MOIUZ.EigenVectors[0][1];
    showerInfoNue->EVecUZ1[1] = MOIUZ.EigenVectors[1][1];

    MOISolution MOIVZ(VStripsU,VStripsZ,VStripsE);
    showerInfoNue->EValVZ0 = MOIVZ.EigenValues[0];
    showerInfoNue->EValVZ1 = MOIVZ.EigenValues[1];
    showerInfoNue->EVecVZ0[0] = MOIVZ.EigenVectors[0][0];
    showerInfoNue->EVecVZ1[0] = MOIVZ.EigenVectors[1][0];
    showerInfoNue->EVecVZ0[1] = MOIVZ.EigenVectors[0][1];
    showerInfoNue->EVecVZ1[1] = MOIVZ.EigenVectors[1][1];

    return;
}



Float_t ANtpShowerInfoAna::RecoShwEnergy(NtpSRShower * ntpShower, Int_t opt, Int_t det){
  //use SR reco
  Float_t shower_ene=0;
  Float_t shwEtemp = GetShwEnergy(ntpShower, opt);
  CandShowerHandle::EShowerType type = GetShwHandleType(opt);

  if(det==1){
    if(opt==-1 || opt == 4) shower_ene = shwEtemp;
    else{
        shower_ene = CorrectShowerEnergyNear(shwEtemp, type);
    }
    return shower_ene;
  }
  if(det==2){
    if(opt==-1 || opt == 4) shower_ene = shwEtemp;
    else{
        shower_ene = CorrectShowerEnergyFar(shwEtemp, type);
    }
    return shower_ene;
  }
  return 0;
}

Float_t ANtpShowerInfoAna::GetShwEnergy(NtpSRShower* ntpShower, Int_t opt)
{
  Float_t shower_ene = 0.0;
  if(opt==-1) shower_ene = ntpShower->ph.gev;
  if(opt==0)  shower_ene = ntpShower->shwph.linCCgev;
  if(opt==1)  shower_ene = ntpShower->shwph.wtCCgev;
  if(opt==2)  shower_ene = ntpShower->shwph.linNCgev;
  if(opt==3)  shower_ene = ntpShower->shwph.wtNCgev;
  if(opt==4)  shower_ene = ntpShower->shwph.EMgev;

  return shower_ene;
}

CandShowerHandle::EShowerType ANtpShowerInfoAna::GetShwHandleType(Int_t opt)
{
  CandShowerHandle::EShowerType type = CandShowerHandle::kCC;
  if(opt==0)  type = CandShowerHandle::kCC;
  if(opt==1)  type = CandShowerHandle::kWtCC;
  if(opt==2)  type = CandShowerHandle::kNC;
  if(opt==3)  type = CandShowerHandle::kWtCC;
 
  return type;
}


Float_t ANtpShowerInfoAna::RecoShwEnergy(Float_t linearCCGeV, SimFlag::SimFlag_t s, const Detector::Detector_t& det, int mode){

  Float_t result = 0;
  // Return Jim's calculation
  Bool_t ok = !ANtpDefVal::IsDefault(linearCCGeV);
  bool isdata = (s == SimFlag::kData);
  if(!ok) return 0;

  // test for shwph.linCCGeV<=0.0
  // if so, assume that this is an R1.16 object
  // and use ph.gev
  result = linearCCGeV;
  if(result> 0.0 )
  {
     result = CorrectShowerEnergy(result,det,CandShowerHandle::kCC,mode,isdata); 
  }
  return result;

}


Float_t ANtpShowerInfoAna::RecoShwEnergyNew(NtpSRShower * ntpShower,
   Int_t opt, VldContext cx) 
{
  EnergyCorrections::WhichCorrection_t corrver = EnergyCorrections::kDefault; 
  Float_t shower_ene=0;

  if(opt==-1) shower_ene = ntpShower->ph.gev;
  if(opt==0)  shower_ene = FullyCorrectShowerEnergy(ntpShower->shwph.linCCgev,CandShowerHandle::kCC,cx,release,corrver);
  if(opt==1)  shower_ene = FullyCorrectShowerEnergy(ntpShower->shwph.wtCCgev,CandShowerHandle::kWtCC,cx,release,corrver);
  if(opt==2)  shower_ene = FullyCorrectShowerEnergy(ntpShower->shwph.linNCgev,CandShowerHandle::kNC,cx,release,corrver);
  if(opt==3)  shower_ene = FullyCorrectShowerEnergy(ntpShower->shwph.wtNCgev,CandShowerHandle::kWtNC,cx,release,corrver);
  if(opt==4)  shower_ene = ntpShower->shwph.EMgev;
  return shower_ene;

  return 0;
}

void ANtpShowerInfoAna::FillGapInformation(NtpSREvent* /*evt*/, NtpSRShower *shw, NtpSRTrack* trk, NtpStRecord * st)
{
   float planes[500];
   int dir[500];

   for(int i = 0; i < 500; i++){
     planes[i] = dir[i] = 0;
   }

   Detector::Detector_t det = st->GetHeader().GetVldContext().GetDetector();

   for(int i = 0; i < 500; i++){
     if(i < 250 || det == Detector::kNear){
       if(i%2 == 0) dir[i] = PlaneView::kU;
       else         dir[i] = PlaneView::kV;
     }else{
       if(i%2 == 0) dir[i] = PlaneView::kV;
       else         dir[i] = PlaneView::kU;
     }
   }                                                                                
   int shwendplane = shw->plane.end;

   if(trk == 0) return;
   if(trk->plane.end <= shwendplane){
     fANtpShowerInfo.longestTrackGapU = 0;
     fANtpShowerInfo.longestTrackGapV = 0;
     fANtpShowerInfo.gapPlanesU = 0;
     fANtpShowerInfo.gapPlanesV = 0;
     fANtpShowerInfo.longestTrackGap = 0;
     fANtpShowerInfo.gapPlanes = 0;
     return;
   }
       
   for(int i=0; i<trk->nstrip; i++){
      Int_t index = SntpHelpers::GetStripIndex(i,trk);
      NtpSRStrip *ntpStrip = SntpHelpers::GetStrip(index,st);
      if(ntpStrip==0) continue;
      if(ntpStrip->plane < shwendplane) continue;
      float charge =  ntpStrip->ph0.sigcor + ntpStrip->ph1.sigcor;
      planes[ntpStrip->plane] = charge;
   }
 
   int ucount = 0;
   int vcount = 0;  
   int count = 0;

   int longestU = 0;
   int longestV = 0;
   int longest = 0;

   int totalU = 0;
   int totalV = 0;   
 
   bool start = true;

   for(int i = shwendplane - 1; i < trk->plane.end + 1; i++)
   {
      if(planes[i] == 0 && start)  continue;  //Keep reading till we find a hit

      if(planes[i] == 0){
        if(dir[i] == PlaneView::kU && !start)   ucount++;
        if(dir[i] == PlaneView::kV && !start)   vcount++;
        if(!start) count++;
        continue;
      }
      
      if(planes[i] > 0){
        start = false;

        if(dir[i] == PlaneView::kU){
           if(ucount > longestU) longestU = ucount;
           totalU += ucount;
           if(count > longest) longest = count;
           ucount = count = 0;
        }
        if(dir[i] == PlaneView::kV){
           if(ucount > longestV) longestV = vcount;
           totalV += vcount;
           if(count > longest) longest = count;
           vcount = count = 0;
        }
      }
   } 

   fANtpShowerInfo.longestTrackGapU = longestU;
   fANtpShowerInfo.longestTrackGapV = longestV;
   fANtpShowerInfo.gapPlanesU = totalU;
   fANtpShowerInfo.gapPlanesV = totalV;
   fANtpShowerInfo.longestTrackGap = longest;
   fANtpShowerInfo.gapPlanes = totalU + totalV;
}
