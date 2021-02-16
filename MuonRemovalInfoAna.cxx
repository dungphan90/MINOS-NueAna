#include <cassert>

#include "TClonesArray.h"
#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSREventSummary.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "TruthHelperNtuple/NtpTHEvent.h"
#include "MuonRemoval/NtpMRRecord.h"
#include "MuonRemoval/NtpMREvent.h"
#include "MuonRemoval/NtpMRTruth.h"
#include "MessageService/MsgService.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/MuonRemovalInfoAna.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "RecoBase/CandShowerHandle.h"
#include "DataUtil/EnergyCorrections.h"
#include "TSystem.h"

#include "NueAna/Module/SetKNNModule.h"  
#include "PhysicsNtuple/Handle.h"
#include "PhysicsNtuple/Factory.h"

MadDpID MuonRemovalInfoAna::dpid;
MadNsID MuonRemovalInfoAna::nsid;
MadAbID MuonRemovalInfoAna::abid;
                                                                                
bool MuonRemovalInfoAna::readabidfile=false;

using namespace EnergyCorrections;

CVSID("$Id: MuonRemovalInfoAna.cxx,v 1.31 2011/07/15 14:09:30 annah1 Exp $");

MuonRemovalInfoAna::MuonRemovalInfoAna(MuonRemovalInfo &mri):
  fMuonRemovalInfo(mri)
{}

MuonRemovalInfoAna::~MuonRemovalInfoAna()
{}

void MuonRemovalInfoAna::Analyze(int event,RecRecordImp<RecCandHeader> *mrobj)
{
  if(mrobj==0){
    return;
  }
  NtpMRRecord *mr = dynamic_cast<NtpMRRecord *>(mrobj);
  Analyze(event,mr,0);
}

void MuonRemovalInfoAna::Analyze(int event,RecRecordImp<RecCandHeader> *mrobj,
				 RecRecordImp<RecCandHeader> *oldstobj)
{
  if(mrobj==0){
    return;
  }
  NtpMRRecord *mr = dynamic_cast<NtpMRRecord *>(mrobj);
  NtpStRecord *oldst = dynamic_cast<NtpStRecord *>(oldstobj);
  Analyze(event,mr,oldst);
}

void MuonRemovalInfoAna::Analyze(int evtn, NtpMRRecord *mrobj, 
				 NtpStRecord *oldstobj)
{
  fMuonRemovalInfo.Reset();
  if(mrobj==0){
    return;
  }

  if(ReleaseType::IsDogwood(release))
     SetCorrectionVersion(EnergyCorrections::kDogwood);
  if(ReleaseType::IsCedar(release))
     SetCorrectionVersion(EnergyCorrections::kCedar);
  if(ReleaseType::IsBirch(release))
     SetCorrectionVersion(EnergyCorrections::kBirch);
  
  //find the best matching rmmu entry for this event
  Int_t best_rmmu = -1;
  Float_t best_com = 0;
  for(int i=0;i<mrobj->mrhdr.nmrevt;i++){
    NtpMREvent *ev = SntpHelpers::GetMREvent(i,mrobj);
    if(ev && ev->best_event==evtn && ev->best_complete>best_com) {
      best_com = ev->best_complete;
      best_rmmu = i;
    }
  }
  if(best_rmmu<0) return;
  
  NtpMREvent *ev = SntpHelpers::GetMREvent(best_rmmu,mrobj);
  NtpMRTruth *tru = SntpHelpers::GetMRTruth(best_rmmu,mrobj);

  if(!ev) return;

  VldContext vc = mrobj->GetHeader().GetVldContext();
  double rangemom = ev->prng;
  double curvemom = ev->pcrv;
  EnergyCorrections::WhichCorrection_t corrver = EnergyCorrections::kDefault;
  

  if(rangemom>0) rangemom=FullyCorrectMomentumFromRange(rangemom,vc,release,corrver);
  curvemom =    FullyCorrectSignedMomentumFromCurvature(curvemom,vc,release,corrver);

  fMuonRemovalInfo.ndigit = ev->ndigit;
  fMuonRemovalInfo.nstrip = ev->nstrip;
  fMuonRemovalInfo.orig_event = ev->orig_event;
  fMuonRemovalInfo.best_purity = ev->best_purity;
  fMuonRemovalInfo.best_complete = ev->best_complete;
  fMuonRemovalInfo.elec_complete = ev->elec_complete;
  fMuonRemovalInfo.best_purity_phw = ev->best_purity_phw;
  fMuonRemovalInfo.best_complete_phw = ev->best_complete_phw;
  fMuonRemovalInfo.elec_complete_phw = ev->elec_complete_phw;
  fMuonRemovalInfo.vtxx = ev->vtxx;
  fMuonRemovalInfo.vtxy = ev->vtxy;
  fMuonRemovalInfo.vtxz = ev->vtxz;
  fMuonRemovalInfo.vtxp = ev->vtxp;
  fMuonRemovalInfo.npln = ev->npln;
  fMuonRemovalInfo.prng = rangemom;
  fMuonRemovalInfo.pcrv = curvemom;
  fMuonRemovalInfo.pvdx = ev->pvdx;
  fMuonRemovalInfo.pvdy = ev->pvdy;
  fMuonRemovalInfo.pvdz = ev->pvdz;
  fMuonRemovalInfo.fitp = ev->fitp;
  fMuonRemovalInfo.endc = ev->endc; 
  fMuonRemovalInfo.pass = ev->pass;
  fMuonRemovalInfo.pmux = ev->pmux;
  fMuonRemovalInfo.pmuy = ev->pmuy;
  fMuonRemovalInfo.pmuz = ev->pmuz;
  fMuonRemovalInfo.mxpl = ev->mxpl;
  fMuonRemovalInfo.vtxdistance = ev->vtxdistance;
  fMuonRemovalInfo.endx = ev->endx;
  fMuonRemovalInfo.endy = ev->endy;
  fMuonRemovalInfo.endz = ev->endz;
  fMuonRemovalInfo.enddistance = ev->enddistance;
  fMuonRemovalInfo.endp = ev->endp;
  fMuonRemovalInfo.zenith = ev->zenith;
  fMuonRemovalInfo.azimuth = ev->azimuth;
  fMuonRemovalInfo.mrmpmux = ev->mrmpmux;
  fMuonRemovalInfo.mrmpmuy = ev->mrmpmuy;
  fMuonRemovalInfo.mrmpmuz = ev->mrmpmuz;
  fMuonRemovalInfo.mrmQ2   = ev->mrmQ2;
  fMuonRemovalInfo.mrmEshw = ev->mrmEshw;
 
  if(tru) {
    fMuonRemovalInfo.nMuonDig = tru->nMuonDig;
    fMuonRemovalInfo.nMuonDigRetained = tru->nMuonDigRetained;
    fMuonRemovalInfo.nShwDig = tru->nShwDig;
    fMuonRemovalInfo.nShwDigRetained = tru->nShwDigRetained;
    fMuonRemovalInfo.nShwDigAtVtx = tru->nShwDigAtVtx;
    fMuonRemovalInfo.nShwDigRetainedAtVtx = tru->nShwDigRetainedAtVtx;
    fMuonRemovalInfo.nShwPE = tru->nShwPE;
    fMuonRemovalInfo.nShwPERetained = tru->nShwPERetained;
    fMuonRemovalInfo.nShwPEAtVtx = tru->nShwPEAtVtx;
    fMuonRemovalInfo.nShwPERetainedAtVtx = tru->nShwPERetainedAtVtx;
    fMuonRemovalInfo.nRetained = tru->nRetained;
    fMuonRemovalInfo.nRetainedMuon = tru->nRetainedMuon;
    fMuonRemovalInfo.nRetainedShw = tru->nRetainedShw;
    fMuonRemovalInfo.nRetainedBoth = tru->nRetainedBoth;
    fMuonRemovalInfo.nPERetained = tru->nPERetained;
    fMuonRemovalInfo.nPERetainedMuon = tru->nPERetainedMuon;
    fMuonRemovalInfo.nPERetainedShw = tru->nPERetainedShw;
    fMuonRemovalInfo.nPERetainedBoth = tru->nPERetainedBoth;
    fMuonRemovalInfo.nRejected = tru->nRejected;
    fMuonRemovalInfo.nRejectedMuon = tru->nRejectedMuon;
    fMuonRemovalInfo.nRejectedShw = tru->nRejectedShw;
    fMuonRemovalInfo.nRejectedBoth = tru->nRejectedBoth;
    fMuonRemovalInfo.nRejShw = tru->nRejShw;
    fMuonRemovalInfo.nRejShwMaxTrk = tru->nRejShwMaxTrk;
    fMuonRemovalInfo.nRejShwFakeTrk = tru->nRejShwFakeTrk;
    fMuonRemovalInfo.nRejShwMix = tru->nRejShwMix;
  }

  //check if an NtpOldRecord was passed
  if(oldstobj==0){
    return;
  }

  NtpSREventSummary *evtSum = &(oldstobj->evthdr);
  NtpSREvent *oldevt = SntpHelpers::GetEvent(ev->orig_event,oldstobj);
  if(oldevt == 0){
      MAXMSG("MuonRemovalAna",Msg::kError,10)<<"Original event ("<<ev->orig_event
          <<") isn't here ("<<evtn<<")anymore, that is odd\n"
          <<" Lost event had comp/pur "
          <<ev->best_purity_phw<<"  "<<ev->best_complete_phw<<endl;
    return;
  }

  NtpSRTrack *oldtrk = 0;
  NtpSRShower *oldshw = 0;
  if(oldevt->ntrack>0)  oldtrk = SntpHelpers::GetTrack(oldevt->trk[0],oldstobj);
  if(oldevt->nshower>0) oldshw = SntpHelpers::GetPrimaryShower(ev->orig_event,oldstobj);
  
  Detector::Detector_t detType = 
    oldstobj->GetHeader().GetVldContext().GetDetector();
  if(detType==Detector::kFar) dpid.SetPHCorrection(1.018);
  Detector::Detector_t fDetectorType = detType;
 

  string reco_version = ReleaseType::AsString(release); 
  string mc_version = ""; 
  if(mrobj->GetHeader().GetVldContext().GetSimFlag()==SimFlag::kMC) 
    mc_version = ReleaseType::AsString(release); 
  if(ReleaseType::IsCarrot(release)) reco_version = "birch"; 
  else if(ReleaseType::IsCedar(release)) reco_version = "cedar"; 
  else if(ReleaseType::IsDogwood(release)) reco_version = "cedar";  //<-- ! no dogwood file! 
  else reco_version = "birch"; 

  BeamType::BeamType_t current_beam = beam;

  if(dpid.ChoosePDFs(detType,BeamType::kLE,
		     reco_version,mc_version))
    fMuonRemovalInfo.orig_cc_pid = 
      dpid.CalcPID(oldtrk,oldevt,evtSum,detType,0);
 
  if(nsid.ChooseWeightFile(detType,current_beam)){
    if(!nsid.GetPID(oldevt,oldtrk,oldshw,oldstobj,detType,fMuonRemovalInfo.orig_nsCCPID))
       fMuonRemovalInfo.orig_nsCCPID=ANtpDefVal::kFloat;
  }

  if(!readabidfile){
    string base=getenv("SRT_PRIVATE_CONTEXT");
    if(base!=""&&base!="."){
      // check if directory exists in SRT_PRIVATE_CONTEXT
      std::string path = base + "/Mad/data";
      void *dir_ptr = gSystem -> OpenDirectory(path.c_str());
      if(!dir_ptr){
        base=getenv("SRT_PUBLIC_CONTEXT");
      }
    }
    else{  base=getenv("SRT_PUBLIC_CONTEXT"); }
    if(base=="") {
      MSG("AnalysisInfoAna",Msg::kFatal)<<"No SRT_PUBLIC_CONTEXT set "
                                        <<"Do not know where to look "
                                        <<"for AB pdf files "<<std::endl;
      assert(false);
    }
    base+="/Mad/data";
    string fname=base;
    string rmc="";
    if(ReleaseType::IsCedar(release)&&mc_version=="daikon"){
      rmc="cedar_daikon";
    }
    else{
      MSG("AnalysisInfoAna",Msg::kWarning)<<"Dont know reco/mc versions "
                                          <<"defaulting to cedar_daikon"<<endl;
      rmc="cedar_daikon";
    }
    string sdet="";
    if(fDetectorType==Detector::kNear){
      sdet="near";
    }
    else if(fDetectorType==Detector::kFar){
      sdet="far";
    }
    else{
      MSG("AnalysisInfoAna",Msg::kWarning)<<"Dont know detector type "
                                          <<"defaulting to far"<<endl;
      sdet="far";
    }
    string sbeam="";
    switch (current_beam){
    case BeamType::kL010z000i:
      sbeam="le0";
      break;
    case BeamType::kL010z170i:
      sbeam="le170";
      break;
    case BeamType::kL010z185i:
      sbeam="le";
      break;
    case BeamType::kL010z200i:
      sbeam="le200";
      break;
    case BeamType::kL100z200i:
      sbeam="pme";
      break;
    case BeamType::kL150z200i:
      sbeam="pme";
      break;
    case BeamType::kL250z200i:
      sbeam="phe";
      break;
    case BeamType::kLE:
      sbeam="le";
      break;
    default:
      MSG("AnalysisInfoAna",Msg::kWarning)<<"Don't know beam type "
                                         <<" defaulting to LE"<<endl;
      sbeam="le";
      break;
    }
                                                                                
    fname+="/ab_pdf_"+sdet+"_"+sbeam+"_"+rmc+".root";
                                                                                
    abid.ReadPDFs(fname.c_str());
    readabidfile=true;
  }
                                                                                
  fMuonRemovalInfo.orig_abCCPID=abid.CalcPID(oldevt,oldstobj);
  fMuonRemovalInfo.orig_roCCPID=LoadROPID(ev->orig_event);
  
  

  ///////////////////////////////////  
  // This block of code is to fix a bug in the original 
  // MuonRemoval/SelectEvent/TrackEndContained function
  // It will recalculate the muon momentum based on the 
  // track end containment. As of 17/11/06 these cuts are
  // the same as those in MuonRemoval
  if(oldtrk) {    

    fMuonRemovalInfo.qp = oldtrk->momentum.qp;
    fMuonRemovalInfo.SigmaQP = oldtrk->momentum.eqp;

    if(detType==Detector::kFar) {
      if(oldtrk->end.plane<=475 &&
	 TMath::Power(oldtrk->end.x,2) + 
	 TMath::Power(oldtrk->end.y,2) <= 12.25){
	fMuonRemovalInfo.endc = true;
	fMuonRemovalInfo.pmux = fMuonRemovalInfo.pvdx*fMuonRemovalInfo.prng;
	fMuonRemovalInfo.pmuy = fMuonRemovalInfo.pvdy*fMuonRemovalInfo.prng;
	fMuonRemovalInfo.pmuz = fMuonRemovalInfo.pvdz*fMuonRemovalInfo.prng;
      }
      else {
	fMuonRemovalInfo.endc = false;
	if(fMuonRemovalInfo.fitp==1){
	  fMuonRemovalInfo.pmux = fMuonRemovalInfo.pvdx*fMuonRemovalInfo.pcrv;
	  fMuonRemovalInfo.pmuy = fMuonRemovalInfo.pvdy*fMuonRemovalInfo.pcrv;
	  fMuonRemovalInfo.pmuz = fMuonRemovalInfo.pvdz*fMuonRemovalInfo.pcrv;
	}
	else {
	  fMuonRemovalInfo.pmux = fMuonRemovalInfo.pvdx*fMuonRemovalInfo.prng;
	  fMuonRemovalInfo.pmuy = fMuonRemovalInfo.pvdy*fMuonRemovalInfo.prng;
	  fMuonRemovalInfo.pmuz = fMuonRemovalInfo.pvdz*fMuonRemovalInfo.prng;
	}
      }
    }
    else if(detType==Detector::kNear) {
      //pitt fiducial volume for ND
      Float_t trk_x  = oldtrk->end.x;
      Float_t trk_y  = oldtrk->end.y;
      Float_t trk_u  = oldtrk->end.u;
      Float_t trk_v  = oldtrk->end.v;
      Float_t trk_z  = oldtrk->end.z;
      Float_t trk_r2 = trk_x*trk_x + trk_y*trk_y;
      Bool_t endc = false;
      if(trk_z<7.0) {
	if( trk_u>0.3 && trk_u<1.8 && trk_v>-1.8 && trk_v<-0.3 && 
	    trk_x<2.4 && trk_r2>0.8*0.8 ) endc = true;
	else endc = false;
      }
      else{
	static const Float_t coil_cut=0.8*0.8;
	static const Float_t x0=0.8;
	static const Float_t y0=0.0;
	static const Float_t a=1.7;
	static const Float_t b=1.4;
	const Float_t xsc = (trk_x-x0)/a; // rescale ellipse to unit circle
	const Float_t ysc = (trk_y-y0)/b;
	if( (sqrt(xsc*xsc + ysc*ysc)<1.0) && 
	    (trk_r2>coil_cut) && (trk_z<15.6) ) endc = true;
	else endc = false;
      }
      
      if(endc){
	fMuonRemovalInfo.endc = true;
	fMuonRemovalInfo.pmux = fMuonRemovalInfo.pvdx*fMuonRemovalInfo.prng;
	fMuonRemovalInfo.pmuy = fMuonRemovalInfo.pvdy*fMuonRemovalInfo.prng;
	fMuonRemovalInfo.pmuz = fMuonRemovalInfo.pvdz*fMuonRemovalInfo.prng;
      }
      else {
	fMuonRemovalInfo.endc = false;
	if(fMuonRemovalInfo.fitp==1){
	  fMuonRemovalInfo.pmux = fMuonRemovalInfo.pvdx*fMuonRemovalInfo.pcrv;
	  fMuonRemovalInfo.pmuy = fMuonRemovalInfo.pvdy*fMuonRemovalInfo.pcrv;
	  fMuonRemovalInfo.pmuz = fMuonRemovalInfo.pvdz*fMuonRemovalInfo.pcrv;
	}
	else {
	  fMuonRemovalInfo.pmux = fMuonRemovalInfo.pvdx*fMuonRemovalInfo.prng;
	  fMuonRemovalInfo.pmuy = fMuonRemovalInfo.pvdy*fMuonRemovalInfo.prng;
	  fMuonRemovalInfo.pmuz = fMuonRemovalInfo.pvdz*fMuonRemovalInfo.prng;
	}
      }
    }
  }
  ///////////////////////////
  
  if(oldshw){
    if(!ANtpDefVal::IsDefault(oldshw->shwph.linCCgev)){
        fMuonRemovalInfo.shwe = FullyCorrectShowerEnergy(oldshw->shwph.linCCgev,
                                   CandShowerHandle::kCC,vc,release,EnergyCorrections::kDefault);
    fMuonRemovalInfo.origShwPlanes = oldshw->plane.n;
    fMuonRemovalInfo.origShwBegPlane = oldshw->plane.beg;
    fMuonRemovalInfo.origShwEndPlane = oldshw->plane.end;
    fMuonRemovalInfo.origShwStrips = oldshw->nstpcnt;
    fMuonRemovalInfo.origShwVtxPlane = oldshw->vtx.plane;
    fMuonRemovalInfo.origShwVtxX = oldshw->vtx.x;
    fMuonRemovalInfo.origShwVtxY = oldshw->vtx.y;
    fMuonRemovalInfo.origShwVtxZ = oldshw->vtx.z;
    }
    else {
    fMuonRemovalInfo.shwe = 0;
    fMuonRemovalInfo.origShwPlanes = 0;
    fMuonRemovalInfo.origShwBegPlane = 0;
    fMuonRemovalInfo.origShwEndPlane = 0;
    fMuonRemovalInfo.origShwStrips = 0;
    fMuonRemovalInfo.origShwVtxPlane = 0;
    fMuonRemovalInfo.origShwVtxX = 0;
    fMuonRemovalInfo.origShwVtxY = 0;
    fMuonRemovalInfo.origShwVtxZ = 0;
    }
  }
  else {
    fMuonRemovalInfo.shwe = 0;
    fMuonRemovalInfo.origShwPlanes = 0;
    fMuonRemovalInfo.origShwBegPlane = 0;
    fMuonRemovalInfo.origShwEndPlane = 0;
    fMuonRemovalInfo.origShwStrips = 0;
    fMuonRemovalInfo.origShwVtxPlane = 0;
    fMuonRemovalInfo.origShwVtxX = 0;
    fMuonRemovalInfo.origShwVtxY = 0;
    fMuonRemovalInfo.origShwVtxZ = 0;

  }
  fMuonRemovalInfo.nrmstp = oldevt->nstrip;
  for(int i=0;i<ev->nstrip;i++){
    for(int j=0;j<oldevt->nstrip;j++){
      if(ev->stp[i]==oldevt->stp[j]) {
	fMuonRemovalInfo.nrmstp -= 1;
	break;
      }
    }
  }

  if(oldstobj->GetHeader().GetVldContext().GetSimFlag()==4){
    NtpTHEvent *oldthev = 
      dynamic_cast<NtpTHEvent *>((*oldstobj->thevt)[ev->orig_event]);
    if(oldthev) {
      fMuonRemovalInfo.orig_evt_purity   = oldthev->purity;
      fMuonRemovalInfo.orig_evt_complete = oldthev->completeall;
    }
  }

}


float MuonRemovalInfoAna::LoadROPID(int evtn){

//   vector<string> list =  Anp::Factory<StorekNNData>::Instance().List();
//   cout<<list.size()<<endl;
 
   //Loading up Rustems variables
   Anp::Handle<StorekNNData> data = Anp::Factory<StorekNNData>::Instance().Get("kNNData");
   if(!data.valid())
   {
      MAXMSG("NueModule", Msg::kError, 1)
          << "NueModule::Reco - Handle<StorekNNData> is invalid, assuming no Rustem variable to run"
          << endl;
      return ANtpDefVal::kFloat;
   }

   float knn_pid;

   data -> SetPrefix("OldSNTP");
   data -> Get(evtn, "knn_pid", knn_pid);
   return knn_pid;
}                                                                              

