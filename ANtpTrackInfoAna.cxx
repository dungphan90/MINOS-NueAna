#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "MessageService/MsgService.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/ANtpTrackInfoAna.h"
#include "AnalysisNtuples/Module/ANtpRecoNtpManipulator.h"
#include "DataUtil/EnergyCorrections.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "VertexFinder/NtpVtxFinder/VertexHelper.h"                                       
//CVSID("$Id: ANtpTrackInfoAna.cxx,v 1.29 2009/07/03 14:45:34 vahle Exp $");

using namespace EnergyCorrections;

ANtpTrackInfoAna::ANtpTrackInfoAna(ANtpTrackInfoNue &anti):
   fANtpTrackInfo(anti)
{}

ANtpTrackInfoAna::~ANtpTrackInfoAna()
{

}

//void ANtpTrackInfoAna::Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord * /*mcobj*/, NtpTHRecord * /*thobj*/)
void ANtpTrackInfoAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
  if(srobj==0){
    return;
  }
  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }

  NtpSREvent *event =0;
  NtpSRTrack *track = 0;

  const RecCandHeader *ntpHeader = &(srobj->GetHeader());
  VldContext vc = ntpHeader->GetVldContext();
  fDetectorType = ntpHeader->GetVldContext().GetDetector();
  //SimFlag::SimFlag_t s = ntpHeader->GetVldContext().GetSimFlag();

  //and now an ugly bit of code to deal with either NtpStRecord or NtpSRRecord
  bool foundst=false;
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(srobj);
  fInfoFiller = new ANtpInfoObjectFiller();

  ANtpRecoNtpManipulator *ntpManipulator = 0;
  NtpSRRecord *sr = 0;
                                                                                
  if(st != 0){
      foundst = true;
      ntpManipulator = new ANtpRecoNtpManipulator(st);
  }
  else{
      sr = dynamic_cast<NtpSRRecord *>(srobj);
      ntpManipulator = new ANtpRecoNtpManipulator(sr);
  }

  if(ReleaseType::IsDogwood(release))
     SetCorrectionVersion(EnergyCorrections::kDogwood);
  if(ReleaseType::IsCedar(release))
     SetCorrectionVersion(EnergyCorrections::kCedar);
  if(ReleaseType::IsBirch(release))
     SetCorrectionVersion(EnergyCorrections::kBirch);
                                                                                
  fInfoFiller->SetStripArray(ntpManipulator->GetStripArray());
  //set up which flags you want to use to determine the primary shower or track
  //a value of 0 for a flag means it will not be used
  ntpManipulator->SetPrimaryTrackCriteria(0,1,0); // nplanes, length, total pulse height
  //get the primary track for the event - if no track is present it
  //returns 0
  ANtpEventManipulator * ntpEventManipulator =
                           ntpManipulator->GetEventManipulator();
                                                                                          
  ntpEventManipulator->SetEventInSnarl(evtn);
  event = ntpEventManipulator->GetEvent();
  track = SntpHelpers::GetPrimaryTrack(evtn, st);
  
  if(track){ 
//          Fill information from the base ANtpTrackInfo  
      fInfoFiller->FillTrackInformation(track, &fANtpTrackInfo);
//          Fill information specific to the ANtpTrackInfoNue class

     int ntrklike = track->plane.ntrklike;
     if(ReleaseType::IsCedar(release)){
       ntrklike = VertexHelper::CalculateNTrkLike(track, st);
     }
     fANtpTrackInfo.trklikePlanes = ntrklike;
                                                                                
     if(event->plane.n){
        fANtpTrackInfo.trklikeRatio= static_cast<Float_t>(ntrklike)/static_cast<Float_t>(event->plane.n);
     }

      FillNueTrackInformation(track, event, &fANtpTrackInfo);
      DetermineSigInOut(track, srobj);
  
      bool IsCont = (IsFidAll(track) || track->momentum.qp==0);
      Int_t Method  = 0;
      if(IsCont) Method = 2; //stoppers use range
      else Method = 1; //Punch Through -> use curvature
      if(track->momentum.qp == 0) Method = 2;
      fANtpTrackInfo.muonEnergyMethod = Method;

      fANtpTrackInfo.phCCGeV = RecoMuEnergyNew(vc);
  }

  if(fInfoFiller){
    delete fInfoFiller;
    fInfoFiller=0;
  }

  if(ntpManipulator){
    delete ntpManipulator;
    ntpManipulator=0;
  }

}

//----------------------------------------------------------------------
void ANtpTrackInfoAna::FillNueTrackInformation(NtpSRTrack *ntpTrack, NtpSREvent *ntpEvent, ANtpTrackInfoNue *trackInfoNue)
{
    if(ntpEvent->ph.mip){
        trackInfoNue->pulseHeightRatio=ntpTrack->ph.mip/ntpEvent->ph.mip; 
    }

    trackInfoNue->phMeu = ntpTrack->ph.mip;
    trackInfoNue->phMip = ntpTrack->ph.mip;    
    trackInfoNue->phNueGeV = trackInfoNue->phMip/MeuPerGeV;

    trackInfoNue->deltaUVVtx = TMath::Abs(ntpTrack->plane.begu - ntpTrack->plane.begv);
    

    return;
}

void ANtpTrackInfoAna::DetermineSigInOut(NtpSRTrack *ntpTrack, RecRecordImp<RecCandHeader> *srobj){
  // compute the amount of signal in the partially instrumented region
  // and the amount in the fully instrumented region
  //
  // this region is defined as:
  // v planes: (strip<=4 || strip>=67)
  // partial u: (strip==0 || strip=63)
  // full u: (strip<=26 || strip>=88)
                                                                                
  Float_t sigfull,sigpart;
  sigfull=sigpart=0;
                                                                                
  // loop over all strips in the event
  // and sum the signals in the partial and full regions
  if(fDetectorType == Detector::kNear){
    for(int i=0; i<ntpTrack->nstrip; i++){
      Int_t index = SntpHelpers::GetStripIndex(i,ntpTrack);
      NtpSRStrip *ntpStrip = SntpHelpers::GetStrip(index,srobj);
      if(ntpStrip==0) continue;
      Int_t pr = NueConvention::InPartialRegion(ntpStrip->plane, ntpStrip->strip);
      float charge =  ntpTrack->stpph1mip[i] + ntpTrack->stpph0mip[i];
      if(pr==1)      {  sigpart += charge; }
      else if(pr==-1){  sigfull += charge; }
    }
  }
                                                                                
  if(fDetectorType == Detector::kFar)
    sigfull = ntpTrack->ph.mip;

  fANtpTrackInfo.trackSignalFull = sigfull;
  fANtpTrackInfo.trackSignalPartial = sigpart;
}


Bool_t ANtpTrackInfoAna::IsFidAll(NtpSRTrack *ntpTrack){
                                                        
  if(ntpTrack == 0) return false;
 
  float x = ntpTrack->end.x;
  float y = ntpTrack->end.y;
  float z = ntpTrack->end.z;

 // new definition - coil hole cuts removed for cedar
  if(fDetectorType ==Detector::kNear) {//near det
    if (!(z<15   &&
          x<2.7  && x>-1.65 &&
          y<1.65 && y>-1.65 &&
          y>(-x)-1.65 && y< x+1.65   &&
          y<(-x)+3.55 && y>x-3.55)) 
      {return false;}
  }
  else if(pow(x,2)+pow(y,2)>14 || ntpTrack->end.plane>475 
           || ntpTrack->end.plane<5) return false;

  return true;
}


Float_t ANtpTrackInfoAna::RecoMuEnergy(SimFlag::SimFlag_t s, const Detector::Detector_t det)
{
  const float mumass=0.10566; // GeV

  bool isdata = (s == SimFlag::kData);
  int opt = fANtpTrackInfo.muonEnergyMethod;

  if(!ANtpDefVal::IsDefault(fANtpTrackInfo.muonEnergyMethod)){
     float mr= fANtpTrackInfo.rangeMomentum;
     float mc= fANtpTrackInfo.fitMomentum;

     mr=CorrectMomentumFromRange(mr,isdata,det);
     mc=CorrectSignedMomentumFromCurvature(mc,isdata,det);

    if(opt==0){
      //return the most appropriate measure of momentum
      // assign opt based on our choice
    }
    else if(opt==1) { //return curvature measurement
        if(fANtpTrackInfo.fitMomentum > 10000.0) return 10000.0;
        else return sqrt(mc*mc + mumass*mumass);
    }
    else if(opt==2) //return range measurement
      return sqrt(mr*mr + mumass*mumass);
    else return 0;
  }
  return 0.;
}

Float_t ANtpTrackInfoAna::RecoMuEnergyNew(VldContext cx, EnergyCorrections::WhichCorrection_t corrver) {

  // using the new version of energy corrections developed for Cedar/Daikon
  corrver = EnergyCorrections::kDefault;

  const float mumass=0.10566;
  int method = fANtpTrackInfo.muonEnergyMethod;

  float mr= fANtpTrackInfo.rangeMomentum;
  float mc= fANtpTrackInfo.fitMomentum;

  if(mr>0) mr=FullyCorrectMomentumFromRange(mr,cx,release,corrver);
  mc=FullyCorrectSignedMomentumFromCurvature(mc,cx,release,corrver);

  if(method == 2) {
     return sqrt(mr*mr+ mumass*mumass);
  }
  if(method == 1){
      if(fANtpTrackInfo.fitMomentum > 10000.0) return 10000.0;
      else return sqrt(mc*mc+ mumass*mumass);
  }
  
  return 0.;
}

