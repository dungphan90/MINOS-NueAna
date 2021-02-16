/**
 *
 * $Id: ANtpAnalysisInfoAna.cxx,v 1.2 2008/11/19 18:22:51 rhatcher Exp $
 *
 * \class ANtpAnalysisInfoAna
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: J. Boehm
 *
 * Created on: Mon Apr 18 17:17:36 2005
 *
 */

#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "MessageService/MsgService.h"
#include "StandardNtuple/NtpStRecord.h"
#include "NueAna/SntpHelpers.h"
#include "NueAna/ANtpAnalysisInfoAna.h"
#include "AnalysisNtuples/ANtpEventInfo.h"
#include "AnalysisNtuples/ANtpTrackInfo.h"
#include "AnalysisNtuples/ANtpShowerInfo.h"
#include "AnalysisNtuples/Module/ANtpRecoNtpManipulator.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "Mad/MadNsID.h"
#include "Mad/MadDpID.h"

MadNsID ANtpAnalysisInfoAna::nsid;
MadDpID ANtpAnalysisInfoAna::dpid;


CVSID("$ID:");

ANtpAnalysisInfoAna::ANtpAnalysisInfoAna(ANtpAnalysisInfoNue &anai):
   fANtpAnalysisInfo(anai)
{
    fDetectorType =  Detector::kUnknown;
}

ANtpAnalysisInfoAna::~ANtpAnalysisInfoAna()
{}

//void ANtpAnalysisInfoAna::Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord * /*mcobj*/, NtpTHRecord * /*thobj*/)
void ANtpAnalysisInfoAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
  MSG("ANtpAnalysisInfoAna",Msg::kDebug)<<"Entering ANtpAnalysisInfoAna::Analyze"<<endl;
	  
	
  if(srobj==0){
      return;
  }  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
        ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
     return;
  }
//	Reset(srobj->GetHeader().GetSnarl(),evtn);

  NtpSREvent *event = 0;
  event = SntpHelpers::GetEvent(evtn,srobj);
  if(!event){
      MSG("ANtpAnalysisInfoAna",Msg::kError)<<"Couldn't get event "<<evtn
         <<" from Snarl "<<srobj->GetHeader().GetSnarl()<<endl;
     return;
  }

 const RecCandHeader *ntpHeader = &(srobj->GetHeader());
 fDetectorType = ntpHeader->GetVldContext().GetDetector();

 NtpSRTrack *track = 0;
 NtpSRShower *shower = 0; 
 NtpSREventSummary *eventSummary;

  //and now an ugly bit of code to deal with either NtpStRecord or NtpSRRecord
  bool foundst=false;
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(srobj);
  ANtpInfoObjectFiller* filla =  new ANtpInfoObjectFiller();
  ANtpRecoNtpManipulator *ntpManipulator = 0;

  //ANtpRecoNtpManipulator ntpManipulator;
  NtpSRRecord *sr = 0; 
 
  if(st != 0)
  {
      foundst = true;
      ntpManipulator =  new ANtpRecoNtpManipulator(st);
      eventSummary = &(st->evthdr);
//.SetRecord(st);
  }
  else{
      sr=dynamic_cast<NtpSRRecord *>(srobj);
      ntpManipulator =  new ANtpRecoNtpManipulator(sr);
      eventSummary = &(sr->evthdr);
//.SetRecord(sr);
  }
 
  //instansiate a NtpHelper object to help you get the info you want
  filla->SetStripArray(ntpManipulator->GetStripArray());
    
  //set up which flags you want to use to determine the primary shower or track
  //a value of 0 for a flag means it will not be used
  ntpManipulator->SetPrimaryShowerCriteria(0,1); // nplanes, total pulse height
  ntpManipulator->SetPrimaryTrackCriteria(0,1,0); // nplanes, length, total pulse height           
  //get the primary shower for the event - if no track is present it
  //returns 0
  ANtpEventManipulator * ntpEventManipulator = 
                            ntpManipulator->GetEventManipulator();
 
  ntpEventManipulator->SetEventInSnarl(evtn);   
  event=ntpEventManipulator->GetEvent();
  shower = ntpEventManipulator->GetPrimaryShower();
  track = ntpEventManipulator->GetPrimaryTrack(); 
  
  FillNueAnalysisInformation(event, track, shower, &fANtpAnalysisInfo, filla);


  Double_t niki_cc_pid;

  
  if(fDetectorType==Detector::kFar ){
       dpid.SetPHCorrection(1.018);
  } 

  BeamType::BeamType_t current_beam = BeamType::kLE;
  if(dpid.ChoosePDFs(fDetectorType,current_beam))
    fANtpAnalysisInfo.dpCCPID = dpid.CalcPID(track, event, eventSummary, fDetectorType, 0);

 //  MadNsID nsid;
  if(nsid.ChooseWeightFile(fDetectorType,current_beam)){
       if(!nsid.GetPID(event,track,shower,st,fDetectorType,fANtpAnalysisInfo.nsCCPID))
           niki_cc_pid=-999;
   }
   else niki_cc_pid=-999;


  if(filla){
     delete filla;
     filla=0;
  }
  if(ntpManipulator){
     delete ntpManipulator;
     ntpManipulator = 0;
  }

  MSG("ANtpAnalysisInfoAna",Msg::kDebug)<<"Leaving ANtpAnalysisInfoAna::Analyze"<<endl;    	  
}
//----------------------------------------------------------------------
void ANtpAnalysisInfoAna::FillNueAnalysisInformation(NtpSREvent *ntpEvent, NtpSRTrack *ntpTrack, NtpSRShower *ntpShower, ANtpAnalysisInfoNue *analysisInfoNue, ANtpInfoObjectFiller * filla)
{
//    ANtpInfoObjectFiller filla;
//    ANtpHeaderInfo *antpHeader       = new ANtpHeaderInfo();
//    ANtpTruthInfoBeam *antpTruth     = new ANtpTruthInfoBeam();
//    ANtpAnalysisInfo *antpAnalysis   = new ANtpAnalysisInfo();
    ANtpEventInfo *antpEvent         = new ANtpEventInfo();
    ANtpTrackInfo *antpTrack         = new ANtpTrackInfo();
    ANtpShowerInfo *antpShower       = new ANtpShowerInfo();

    filla->FillEventInformation(ntpEvent, antpEvent);
    if(ntpTrack != 0)  filla->FillTrackInformation(ntpTrack,antpTrack);
    if(ntpShower != 0) filla->FillShowerInformation(ntpShower,antpShower);      


    analysisInfoNue->isFullyContained  = -10;
    analysisInfoNue->passesCuts        = 1;
    analysisInfoNue->recoEventLength   = TMath::Abs(antpEvent->endPlane -
		       antpEvent->begPlane);

    if(ntpTrack != 0){
     analysisInfoNue->recoTrackLength   = TMath::Abs(antpTrack->endPlane -
		       antpTrack->begPlane);
     analysisInfoNue->recoTrackMomentum = antpTrack->fitMomentum;
     analysisInfoNue->recoTrackRange    = antpTrack->rangeMomentum;
     analysisInfoNue->recoSigmaQoverP   = antpTrack->sigmaQoverP;
      
     analysisInfoNue->recoMuEnergy      = RecoMuEnergy(ntpTrack);

     if(ntpShower != 0){
       analysisInfoNue->recoShowerEnergy  = RecoShwEnergy(ntpShower);
       analysisInfoNue->recoNuEnergy      = ( analysisInfoNue->recoMuEnergy +
			analysisInfoNue->recoShowerEnergy );
     }
     analysisInfoNue->recoQENuEnergy     = RecoQENuEnergy(ntpTrack);
     analysisInfoNue->recoQEQ2           = RecoQEQ2(ntpTrack);

     if (analysisInfoNue->recoNuEnergy>0) {
        analysisInfoNue->recoHadronicY    = ( analysisInfoNue->recoShowerEnergy /
		       analysisInfoNue->recoNuEnergy );
     }
	
       analysisInfoNue->recoMuDCosZVtx    = RecoMuDCosZ(ntpTrack);
       analysisInfoNue->recoNuDCos        = RecoMuDCosNeu(ntpTrack);
    }

    analysisInfoNue->recoVtxX          = antpEvent->vtxX;
    analysisInfoNue->recoVtxY          = antpEvent->vtxY;
    analysisInfoNue->recoVtxZ          = antpEvent->vtxZ;

    analysisInfoNue->isFullyContained  = 
         IsFidAll(antpEvent->vtxX, antpEvent->vtxY, antpEvent->vtxZ, ntpEvent);
//       analysisInfoNue->isFullyContained  = 0;  
    //fiducial criteria on vtx

    analysisInfoNue->inFiducialVolume = IsFidVtxEvt(ntpEvent);
//       analysisInfoNue->inFiducialVolume = 0;
    

    if(antpEvent){
      delete antpEvent;
      antpEvent=0;
    }

    if(antpTrack){
      delete antpTrack;
      antpTrack=0;
    }
    if(antpShower){
      delete antpShower;
      antpShower=0;
    }
  
    return;
}

//*******************************************************
Float_t ANtpAnalysisInfoAna::RecoMuEnergy(NtpSRTrack *ntpTrack){

     if(IsFidAll(ntpTrack->vtx.x, ntpTrack->vtx.y, ntpTrack->vtx.z)
                    || ntpTrack->momentum.qp==0) {
         return sqrt(ntpTrack->momentum.range*ntpTrack->momentum.range
                  + 0.10555*0.10555);
     }
     else {
       if(ntpTrack->momentum.qp < 1e-5) return 10000.0;
       else 
          return sqrt(1./(ntpTrack->momentum.qp*ntpTrack->momentum.qp)
                                + 0.10555*0.10555);
     }

     return 0;
}

//*******************************************************
Float_t ANtpAnalysisInfoAna::RecoShwEnergy(NtpSRShower *ntpShower){
	  //use SR reco
  Float_t theGeV=ntpShower->ph.gev;
  return theGeV;
}

//*******************************************************
Float_t ANtpAnalysisInfoAna::RecoQENuEnergy(NtpSRTrack *ntpTrack){
	
	Float_t nucleonMass = 0.93956563; //mass of neutron by default
	if(GetChargeSign(ntpTrack)==1) nucleonMass = 0.93827231; //proton mass for nubar
	Float_t muonEnergy = RecoMuEnergy(ntpTrack);
	Float_t muonMass = 0.10555;
        if(muonEnergy < muonMass){
          MSG("ANtpAnalysisInfoAna",Msg::kError)
                 << "muon Energy < muon mass | Large reco failure"
                 << endl;
          return ANtpDefVal::kFloat;
        }

        if(TMath::Abs(muonEnergy) > 1e10){
           MSG("ANtpAnalysisInfoAna",Msg::kError)
                  << "muon Energy too big "
		  << muonEnergy<< " stopping"
                  << endl;
           return ANtpDefVal::kFloat;
        }
		
	
	Float_t muonMomentum = sqrt(muonEnergy*muonEnergy - muonMass*muonMass);
	Float_t costhbl = RecoMuDCosNeu(ntpTrack);

	Float_t Eneu = ANtpDefVal::kFloat;
        if(nucleonMass - muonEnergy + muonMomentum*costhbl > 1e-8) 
	Eneu = (nucleonMass*muonEnergy - muonMass*muonMass/2.)
		/(nucleonMass - muonEnergy + muonMomentum*costhbl);
	
	return Eneu;
}
                                                         
//********************************************************       
Float_t ANtpAnalysisInfoAna::RecoQEQ2(NtpSRTrack *ntpTrack){

     Float_t Eneu = RecoQENuEnergy(ntpTrack);
     if(ANtpDefVal::IsDefault(Eneu)) return ANtpDefVal::kFloat;
     Float_t muonEnergy = RecoMuEnergy(ntpTrack);
     Float_t muonMass = 0.10555;
     Float_t muonMomentum = sqrt(muonEnergy*muonEnergy - muonMass*muonMass);
     Float_t costhbl = RecoMuDCosNeu(ntpTrack);
     return -2.*Eneu*(muonEnergy-muonMomentum*costhbl)+muonMass*muonMass;
  
}

//*******************************************************
Int_t ANtpAnalysisInfoAna::GetChargeSign(NtpSRTrack *ntpTrack){
  if(RecoMuQP(ntpTrack)>0) return 1;
  return -1;
}

//*******************************************************
Float_t ANtpAnalysisInfoAna::RecoMuQP(NtpSRTrack *ntpTrack){
  return ntpTrack->momentum.qp;
}

//*******************************************************
Float_t ANtpAnalysisInfoAna::RecoMuDCosNeu(NtpSRTrack *ntpTrack){ //ds_mu/ds_neu
  Float_t bl_z = TMath::Cos(TMath::Pi()*3./180.); //3degree beam
  Float_t bl_y = sqrt(1. - bl_z*bl_z);
  Float_t costhbl = ntpTrack->vtx.dcosz*bl_z + ntpTrack->vtx.dcosy*bl_y;
	
  return  costhbl;
}

//*******************************************************
Float_t ANtpAnalysisInfoAna::RecoMuDCosZ(NtpSRTrack *ntpTrack){ //dz/ds
   return ntpTrack->vtx.dcosz;
}

//*******************************************************
Bool_t ANtpAnalysisInfoAna::IsFidVtx(NtpSRTrack *ntpTrack){
  if(fDetectorType == Detector::kFar){
	  
     if(ntpTrack->vtx.z<0.5 || ntpTrack->vtx.z>29.4 ||   //ends
	  (ntpTrack->vtx.z<16.5&&ntpTrack->vtx.z>14.5) ||  //between SMs
	  sqrt((ntpTrack->vtx.x*ntpTrack->vtx.x)           //radial cut
               +(ntpTrack->vtx.y*ntpTrack->vtx.y))>3.5 ||
             sqrt((ntpTrack->vtx.x*ntpTrack->vtx.x)           //radial cut
                 +(ntpTrack->vtx.y*ntpTrack->vtx.y))<0.4) return false;
   
  }
  else if(fDetectorType == Detector::kNear){

      if(ntpTrack->vtx.z<0.5 || ntpTrack->vtx.z>6.5 ||
           sqrt(((ntpTrack->vtx.x-1.3)*(ntpTrack->vtx.x-1.3)) +
              (ntpTrack->vtx.y*ntpTrack->vtx.y))>1) return false;
       }                                                                                                                    
  return true;
}

//*******************************************************
Int_t ANtpAnalysisInfoAna::IsFidVtxEvt(NtpSREvent *ntpEvent){

  Bool_t contained = true;
	
  Int_t test = 0;
  MSG("ANtpAnalysisInfoAna",Msg::kDebug) << "DetectorType " << fDetectorType << endl;
  if(fDetectorType == Detector::kNear)
       test = InsideNearFiducial(ntpEvent->vtx.x, ntpEvent->vtx.y,
		                  ntpEvent->vtx.z);
  if(fDetectorType == Detector::kFar)
       test = InsideFarFiducial(ntpEvent->vtx.x, ntpEvent->vtx.y,
		                  ntpEvent->vtx.z);
  if(test <= 0) contained = false;
  MSG("ANtpAnalysisInfoAna",Msg::kDebug) << " IsFidVtxEvt " << test << endl;
  return test;
}

Int_t ANtpAnalysisInfoAna::IsFidVtxEvt(NtpSREvent *ntpEvent, Int_t detType){

  Bool_t contained = true;
	
  Int_t test = 0;
  MSG("ANtpAnalysisInfoAna",Msg::kDebug) << "DetectorType " << detType << endl;
  if(detType == Detector::kNear)
       test = InsideNearFiducial(ntpEvent->vtx.x, ntpEvent->vtx.y,
		                  ntpEvent->vtx.z);
  if(detType == Detector::kFar)
       test = InsideFarFiducial(ntpEvent->vtx.x, ntpEvent->vtx.y,
		                  ntpEvent->vtx.z);
  if(test <= 0) contained = false;
  
  MSG("ANtpAnalysisInfoAna",Msg::kDebug) << " IsFidVtxEvt " << test << endl;
  return test;
}

void ANtpAnalysisInfoAna::Set3DHit(DeqFloat_t &x
                             , DeqFloat_t &y
                             , DeqFloat_t &z
                             , DeqFloat_t &e){
    	
     if(x.size()!=0 && y.size()!=0 && z.size()!=0 && e.size()!=0){
 
           fX=x;
           fY=y;
           fZ=z;
           fE=e;
     }
}


Int_t ANtpAnalysisInfoAna::IsFidAll(Float_t vtxX, Float_t vtxY, Float_t vtxZ, NtpSREvent *event)
{
  
  Bool_t contained = true;
  Int_t test = 0;
  if(fX.size()==0 || fY.size()==0 || fZ.size()==0 || fE.size()==0){
//       MAXMSG("ANtpAnalysisInfoAna",Msg::kWarning, 1)
//            << "3D Hits not set for event,"
//            << "cannot evaluate isFullyContained."
//            << endl;
        return ANtpDefVal::kInt;
    }   	

  for(UInt_t i = 0; i < fX.size() && contained; i++){
     test = 0;
     Float_t x = fX[i] + vtxX;
     Float_t y = fY[i] + vtxY;
     Float_t z = fZ[i] + vtxZ;
     if(fDetectorType == Detector::kNear)
	     test = InsideNearFiducial(x, y, z);
     if(fDetectorType == Detector::kFar)
             test = InsideFarFiducial(x, y, z);
     if(test <= 0) contained = false;
  }

   if(event) 
   if(contained && event->plane.n > 66) //this is approximately 4 meteres which is the 3D hit cutoff
   {
     if(fDetectorType == Detector::kNear)
                test = InsideNearFiducial(event->end.x, event->end.y, event->end.z);
     if(fDetectorType == Detector::kFar)
                test = InsideFarFiducial(event->end.x, event->end.y, event->end.z);
     if(test <= 0) contained = false;
   }	  
  
  if(contained) test = 1;
  return test;
}
     
Int_t ANtpAnalysisInfoAna::InsideFarFiducial(Float_t x, Float_t y, Float_t z)
{
  Float_t SuperModule1Beg =  0.35;
  Float_t SuperModule2Beg = 16.20;
  Float_t SuperModule1End = 14.57;
  Float_t SuperModule2End = 29.62;

  Float_t radialInner = 0.40;
  Float_t radialOuter = 3.87;
  Bool_t zContained = false;
  Bool_t xyContained = false;

  Float_t r = TMath::Sqrt(x*x + y*y);
    
  if( (z >= SuperModule1Beg && z <=SuperModule1End) || 
      (z >= SuperModule2Beg && z <=SuperModule2End) )
     zContained = true;

  if( r >= radialInner && r <= radialOuter)
     xyContained = true;
  
  Int_t retVal = 0;
  if(zContained && xyContained) retVal = 1;
  if(!zContained) retVal = -1;
  if(!xyContained) retVal -= 2;
 
  return retVal;  //  1 contained, -1 out of bounds z 
                  //  -2 oob xy, -3 oob both
}

Int_t ANtpAnalysisInfoAna::InsideNearFiducial(Float_t x, Float_t y, Float_t z)
{
  Float_t SuperModule1Beg = 0.40;
  Float_t SuperModule1End = 6.50;
                                                                                
  Float_t radialInner = 0;
  Float_t radialOuter = 1;
  Float_t xCenter = 1.4885;
  Float_t yCenter = 0.1397;
  Bool_t zContained = false;
  Bool_t xyContained = false;

  Float_t r = TMath::Sqrt((x-xCenter)*(x-xCenter) + (y-yCenter)*(y-yCenter));
  if( z >= SuperModule1Beg && z <=SuperModule1End)
     zContained = true;
                                                                                
  if( r >= radialInner && r <= radialOuter)
     xyContained = true;

  Int_t retVal = 0;
  if(zContained && xyContained) retVal = 1;
  if(!zContained) retVal = -1;
  if(!xyContained) retVal -= 2;		  
  
  return retVal;
}

