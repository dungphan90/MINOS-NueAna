#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TClonesArray.h"
#include "TMinuit.h"
#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "CandNtupleSR/NtpSRSlice.h"
#include "MessageService/MsgService.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/TimingVarsAna.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "Calibrator/CalMIPCalibration.h"
#include "Plex/PlexStripEndId.h"
#include "RecoBase/PropagationVelocity.h"
#include "UgliGeometry/UgliStripHandle.h"
#include "UgliGeometry/UgliGeomHandle.h"
//#include "Midad/Base/Mint.h"

CVSID("$Id: TimingVarsAna.cxx,v 1.9 2009/07/03 14:45:34 vahle Exp $");

#include "DatabaseInterface/DbiResultPtr.tpl"

TimingVarsAna::TimingVarsAna(TimingVars &tv):
   fTimingVars(tv)
{

//    TimeHst = new THStack("TimeHist","Digit Times(ns) Red - Track, Blue - Shower.");
    TimeHstTrk = new TH1F();
    TimeHstTrk->SetName("TimeHstTrk");
    TimeHstShw = new TH1F();
    TimeHstShw->SetName("TimeHstShw");

    TotalHst = new TH1F();
    TotalHst->SetName("TotalHst");

//    plot = new TCanvas("plot","plot",200,10,900,600);

}


TimingVarsAna::~TimingVarsAna()
{
  if(TimeHstTrk) delete TimeHstTrk;
  if(TimeHstShw) delete TimeHstShw;
  if(TotalHst) delete TotalHst;
}

void TimingVarsAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj){
  if(srobj==0){
    return;
  }
  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }

  //Reset(srobj->GetHeader().GetSnarl(),evtn);
  NtpStRecord *st = dynamic_cast<NtpStRecord *>(srobj);
  NtpSREvent *event = SntpHelpers::GetEvent(evtn,srobj);

 
  if(!event){
      MSG("ShwfitAna",Msg::kError)<<"Couldn't get event "<<evtn
				   <<" from Snarl "<<srobj->GetHeader().GetSnarl()<<endl;
      return;
   }

   VldContext vc=st->GetHeader().GetVldContext();

   //VldTimeStamp vldts = vc.GetTimeStamp();
   UgliGeomHandle ugh(vc);
   //UgliGeomHandle ugh = gMint->GetUgliGeomHandle();
   if (! ugh.IsValid()) {
     MSG("NueDisplayModule",Msg::kWarning) << "Got invalid Ugli\n";
     //return;
   }

  fDetectorType = vc.GetDetector();
  ReleaseType::Release_t fRel;
  string relName = st->GetTitle();
  string reco = relName.substr(0,relName.find_first_of("("));
  if(reco == "CEDAR") fRel = ReleaseType::kCedar;
  if(reco == "DOGWOOD") fRel = ReleaseType::kDogwood;
  else fRel = ReleaseType::kBirch;    
  int ntrks = event->ntrack;
      //timing histogram
    //find range
    double tmin=0, tmax=0;
    bool first = true;
    NtpSRSlice *slice = SntpHelpers::GetSlice(event->slc,st);

    //shamelessly stolen from Mad
    if (slice){//slice
      float highest_plane = 0;
      float lowest_plane = 500;
      float highest_strip0 = 0;
      float lowest_strip0 = 192;
      float highest_strip1 = 0;
      float lowest_strip1 = 192;

      float highest_z = 0.;
      float lowest_z = 30.;
      float highest_t0 = -4.0;
      float lowest_t0 = 4.0;
      float highest_t1 = -4.0;
      float lowest_t1 = 4.0;

      for (int i = 0; i<slice->nstrip; i++){//loop over strips
	NtpSRStrip *strip = SntpHelpers::GetStrip(slice->stp[i],srobj);
        if(strip == 0) continue;
	double t = strip->time0;
	if (t>-100){
	  if (first){
	    tmin = tmax = t;
	    first = false;
	  }
	  if (t < tmin) tmin = t;
	  if (t > tmax) tmax = t;
	}
	t = strip->time1;
	if (t>-100){
	  if (first){
	    tmin = tmax = t;
	    first = false;
	  }
	  if (t < tmin) tmin = t;
	  if (t > tmax) tmax = t;
	}
	int tempo_pln = strip->plane;
	int tempo_stp = strip->strip;
	float tempo_tpos = strip->tpos;
	if(tempo_pln<lowest_plane) {
	  lowest_plane=tempo_pln;
	  lowest_z=strip->z;
	}
	if(tempo_pln>highest_plane) {
	  highest_plane=tempo_pln;
	  highest_z=strip->z;
	}

	if (strip->planeview==PlaneView::kU){
	 // fSlcUZ->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + strip->ph1.sigcor)/SIGMAPMEU);
	  if(tempo_tpos<lowest_t0) {
	    lowest_strip0=tempo_stp;
	    lowest_t0=tempo_tpos;
	  }
	  if(tempo_tpos>highest_t0) {
	    highest_strip0=tempo_stp;
	    highest_t0=tempo_tpos;
	  }
	}
	else if (strip->planeview==PlaneView::kV){
	  //fSlcVZ->Fill(strip->z,strip->tpos,(strip->ph0.sigcor + strip->ph1.sigcor)/SIGMAPMEU);
	  if(tempo_tpos<lowest_t1) {
	    lowest_strip1=tempo_stp;
	    lowest_t1=tempo_tpos;
	  }
	  if(tempo_tpos>highest_t1) {
	    highest_strip1=tempo_stp;
	    highest_t1=tempo_tpos;
	  }
	}
      }
      if(lowest_plane-10>=0) {
	lowest_plane-=10;
	lowest_z-=10.*0.06;
      }
      else {
	lowest_plane=0;
	lowest_z=0.;
      }

      if(lowest_strip0-5>=0) {
	lowest_strip0-=5;
	lowest_t0-=5.*0.041;
      }
      else {
	lowest_strip0=0;
	lowest_t0=-4.0;
      }

      if(lowest_strip1-5>=0) {
	lowest_strip1-=5;
	lowest_t1-=5.*0.041;
      }
      else {
	lowest_strip1=0;
	lowest_t1=-4.0;
      }

      if(highest_plane+10<=485) {
	highest_plane+=10;
	highest_z+=10.*0.06;
      }
      else {
	highest_plane=485;
	highest_z=30.;
      }

      if(highest_strip0+5<=191) {
	highest_strip0+=5;
	highest_t0+=5.*0.041;
      }
      else {
	highest_strip0=191;
	highest_t0=4.0;
      }

      if(highest_strip1+5<=191) {
	highest_strip1+=5;
	highest_t1+=5.*0.041;
      }
      else {
	highest_strip1=191;
	highest_t1=4.0;
      }
    }
    // give some buffer at either end...
    tmin -= 50e-9;
    tmax += 20e-9;

    double eps = 1.0e-8;
    if (tmin == tmax) { tmin -= eps; tmax += eps; }

    // for far detector, suppress display of pre-trigger time interval
    if (fDetectorType == Detector::kFar){
      if ((tmax - tmin)*1e9>1500) tmin = tmax - 1500e-9;
    }

    TimeHstTrk->SetBins(50,0,100);
    TimeHstShw->SetBins(50,0,100);
    TotalHst->SetBins(50,0,100);
    //record track hits
    //int ntrks = event->ntrack;
    if (ntrks){//if(ntrks)
      int trkidx = -1;
      int trkplanes = -1;
      for (int i = 0; i<ntrks; i++){
	int index = SntpHelpers::GetTrackIndex(i,event);
	NtpSRTrack *track = SntpHelpers::GetTrack(index,st);
        if(track == 0) continue;
	if (track->plane.n>trkplanes){
	  trkplanes = track->plane.n;
	  trkidx = index;
	}
	for (int istp = 0; istp<track->nstrip; istp++){
	  //ftrkshw->Fill(track->stpx[istp],track->stpy[istp],1);
	  NtpSRStrip *strip = SntpHelpers::GetStrip(track->stp[istp],st);
          if(strip == 0) continue;
	  TimeHstTrk->Fill((track->stpt0[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph0.pe);
	  TimeHstTrk->Fill((track->stpt1[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph1.pe);
	  TotalHst->Fill((track->stpt0[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph0.pe);
	  TotalHst->Fill((track->stpt1[istp]-tmin-(track->ds-track->stpds[istp])/3e8)/1e-9,strip->ph1.pe);
	  //record track strip information
	  if(strip->planeview==PlaneView::kU){
	    //fTrkUZ->Fill(strip->z,strip->tpos,100000);//paranoia
	  }
	  else if(strip->planeview==PlaneView::kV){
	    //fTrkVZ->Fill(strip->z,strip->tpos,100000);
	  }
	}
      }
      vector<double> spathLength;
      vector<double> st0;
      if (trkidx != -1 && fDetectorType == Detector::kNear){
	NtpSRTrack *track = SntpHelpers::GetTrack(trkidx,st);
	for (int istp = 0; istp<track->nstrip; istp++){
	  NtpSRStrip *strip = SntpHelpers::GetStrip(track->stp[istp],st);
          if(strip == 0) continue;
	  spathLength.push_back(track->ds-track->stpds[istp]);
	  PlexStripEndId seid(Detector::kNear,strip->plane,strip->strip,StripEnd::kWest);
	  UgliStripHandle stripHandle = ugh.GetStripHandle(seid);
	  float halfLength = stripHandle.GetHalfLength();
	  const TVector3 ghitxyz(track->stpx[istp],track->stpy[istp],track->stpz[istp]);
	  TVector3 lhitxyz = stripHandle.GlobalToLocal(ghitxyz);
	  float fiberDist = (halfLength - lhitxyz.x() + stripHandle.ClearFiber(StripEnd::kWest) + stripHandle.WlsPigtail(StripEnd::kWest));
	  //using strip time, I don't understand track time TJ
	  //st0.push_back(track->stpt1[istp]-fiberDist/PropagationVelocity::Velocity());
	  st0.push_back(strip->time1-fiberDist/PropagationVelocity::Velocity());
	  //cout<<strip->plane<<" "<<strip->strip<<" "<<track->ds-track->stpds[istp]<<Form(" %.9f %.9f %.9f",track->stpt1[istp],strip->time1,strip->time1-fiberDist/PropagationVelocity::Velocity())<<endl;
	}
	float trms1 = 0;
	float trms2 = 0;

	trms1/=st0.size();
	trms2/=st0.size();
	trms1=sqrt(trms1)*1e9;
	trms2=sqrt(trms2)*1e9;
      }
    }

    //record shower hits
    int nshws = event->nshower;
    if (nshws){
      for (int i = 0; i<nshws; i++){
	int index = SntpHelpers::GetShowerIndex(i,event);
	NtpSRShower *shower = SntpHelpers::GetShower(index,st);
        if(shower == 0) continue;
	for (int istp = 0; istp<shower->nstrip; istp++){
	  NtpSRStrip *strip = SntpHelpers::GetStrip(shower->stp[istp],st);
          if(strip == 0) continue;
	  if(ReleaseType::IsCedar(fRel)||ReleaseType::IsDogwood(fRel)){
	    TimeHstShw->Fill((shower->stpt0[istp]-tmin)/1e-9,strip->ph0.pe);
	    TimeHstShw->Fill((shower->stpt1[istp]-tmin)/1e-9,strip->ph1.pe);
            TotalHst->Fill((shower->stpt0[istp]-tmin)/1e-9,strip->ph0.pe);
            TotalHst->Fill((shower->stpt1[istp]-tmin)/1e-9,strip->ph1.pe);
            //cout<<"Bin: "<<(shower->stpt1[istp]-tmin)/1e-9<<" Weight: "<<strip->ph0.pe+strip->ph1.pe<<endl;
          }else {
	    TimeHstShw->Fill((strip->time0-tmin)/1e-9,strip->ph0.pe);
	    TimeHstShw->Fill((strip->time1-tmin)/1e-9,strip->ph1.pe);
	    TotalHst->Fill((strip->time0-tmin)/1e-9,strip->ph0.pe);
	    TotalHst->Fill((strip->time1-tmin)/1e-9,strip->ph1.pe);
	  }
	  //record shower strip information
	  if(strip->planeview==PlaneView::kU){
	  //  fShwUZ->Fill(strip->z,strip->tpos,100000);//paranoia
	  }
	  else if(strip->planeview==PlaneView::kV){
	   // fShwVZ->Fill(strip->z,strip->tpos,100000);
	  }
	}
      }
    }
//Max bins and contents
    Int_t ShwMaxBin = TimeHstShw->GetMaximumBin();
    Float_t ShwMaxBinCont = TimeHstShw->GetBinContent(ShwMaxBin);

    Int_t TrkMaxBin = TimeHstTrk->GetMaximumBin();
    Float_t TrkMaxBinCont = TimeHstTrk->GetBinContent(TrkMaxBin);

    Int_t TotalMaxBin = TotalHst->GetMaximumBin();
    Float_t TotalMaxBinCont = TotalHst->GetBinContent(TotalMaxBin);

//Total Pulse height
    Float_t ShwEventPulseHeight = TimeHstShw->Integral();
    Float_t TrkEventPulseHeight = TimeHstTrk->Integral();
    Float_t TotalEventPulseHeight = TotalHst->Integral();

//Second max bin contents
    Float_t ShwSecondMaxBin;
    if(TimeHstShw->GetBinContent(ShwMaxBin+1)>TimeHstShw->GetBinContent(ShwMaxBin-1)){
      ShwSecondMaxBin = TimeHstShw->GetBinContent(ShwMaxBin+1);
    }else{
      ShwSecondMaxBin = TimeHstShw->GetBinContent(ShwMaxBin-1);
    }
    Float_t TrkSecondMaxBin;
    if(TimeHstTrk->GetBinContent(TrkMaxBin+1)>TimeHstTrk->GetBinContent(TrkMaxBin-1)){
      TrkSecondMaxBin = TimeHstTrk->GetBinContent(TrkMaxBin+1);
    }else{
      TrkSecondMaxBin = TimeHstTrk->GetBinContent(TrkMaxBin-1);
    }
    Float_t TotalSecondMaxBin;
    if(TotalHst->GetBinContent(TotalMaxBin+1)>TotalHst->GetBinContent(TotalMaxBin-1)){
      TotalSecondMaxBin = TotalHst->GetBinContent(TotalMaxBin+1);
    }else{
      TotalSecondMaxBin = TotalHst->GetBinContent(TotalMaxBin-1);
    }

//Percent of pulse height in Max 2 bins
    Float_t TrkPercentofTotal = -9999.99;
    Float_t ShwPercentofTotal = -9999.99;
    Float_t TotalPercentofTotal = -9999.99;

    if (ShwEventPulseHeight>0){
      ShwPercentofTotal = (ShwSecondMaxBin + ShwMaxBinCont)/(ShwEventPulseHeight);
    }
    if (TrkEventPulseHeight>0){
      TrkPercentofTotal = (TrkSecondMaxBin + TrkMaxBinCont)/(TrkEventPulseHeight);
    }
    if (TotalEventPulseHeight>0){
      TotalPercentofTotal = (TotalSecondMaxBin + TotalMaxBinCont)/(TotalEventPulseHeight);
    }

//how many of the two bins prior to the max bin are at least 20% of the max bin ph
    Int_t ShwbinsPrior = 0;
    if(TimeHstShw->GetBinContent(ShwMaxBin-1)>(0.2*ShwMaxBinCont)){ShwbinsPrior++;}
    if(TimeHstShw->GetBinContent(ShwMaxBin-2)>(0.2*ShwMaxBinCont)){ShwbinsPrior++;}
    Int_t TrkbinsPrior = 0;
    if(TimeHstTrk->GetBinContent(TrkMaxBin-1)>(0.2*TrkMaxBinCont)){TrkbinsPrior++;}
    if(TimeHstTrk->GetBinContent(TrkMaxBin-2)>(0.2*TrkMaxBinCont)){TrkbinsPrior++;}
    Int_t TotalbinsPrior = 0;
    if(TotalHst->GetBinContent(TotalMaxBin-1)>(0.2*TotalMaxBinCont)){TotalbinsPrior++;}
    if(TotalHst->GetBinContent(TotalMaxBin-2)>(0.2*TotalMaxBinCont)){TotalbinsPrior++;}
    
///how many bins are over 40 ph
    Int_t ShwBinsOver40 = 0;
    Int_t TrkBinsOver40 = 0;
    Int_t TotalBinsOver40 = 0;
    for (int k=0;k<50;k++){
      if(TimeHstShw->GetBinContent(k)>20)ShwBinsOver40++;
      if(TimeHstTrk->GetBinContent(k)>20)TrkBinsOver40++;
      if(TotalHst->GetBinContent(k)>20)TotalBinsOver40++;
    }

    fTimingVars.ShwMaxBin = ShwMaxBin;
    fTimingVars.TrkMaxBin = TrkMaxBin;
    fTimingVars.TotalMaxBin = TotalMaxBin;

    fTimingVars.ShwMaxBinCont = ShwMaxBinCont;
    fTimingVars.TrkMaxBinCont = TrkMaxBinCont;
    fTimingVars.TotalMaxBinCont = TotalMaxBinCont;

    fTimingVars.ShwSecondMaxBin = ShwSecondMaxBin;
    fTimingVars.TrkSecondMaxBin = TrkSecondMaxBin;
    fTimingVars.TotalSecondMaxBin = TotalSecondMaxBin;

    fTimingVars.ShwEventPulseHeight = ShwEventPulseHeight;
    fTimingVars.TrkEventPulseHeight = TrkEventPulseHeight;
    fTimingVars.TotalEventPulseHeight = TotalEventPulseHeight;

    fTimingVars.ShwPercentofTotal = ShwPercentofTotal;
    fTimingVars.TrkPercentofTotal = TrkPercentofTotal;
    fTimingVars.TotalPercentofTotal = TotalPercentofTotal;

    fTimingVars.ShwbinsPrior = ShwbinsPrior;
    fTimingVars.TrkbinsPrior = TrkbinsPrior;
    fTimingVars.TotalbinsPrior = TotalbinsPrior;

    fTimingVars.ShwBinsOver40 = ShwBinsOver40;
    fTimingVars.TrkBinsOver40 = TrkBinsOver40;
    fTimingVars.TotalBinsOver40 = TotalBinsOver40;
    /*
    char name[100];
    sprintf(name,"/home/danche/ViewScan5/plot/ShwHst_snarl_%d.png",st->GetHeader().GetSnarl());

    TimeHstShw->Draw();
    plot->Print(name);
    */
}

