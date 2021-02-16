#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"
#include "FracVarAna.h"
#include "NueAna/mlpANN.h"
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

CVSID("$Id: FracVarAna.cxx,v 1.38 2009/08/05 01:50:59 jjling Exp $");

TMultiLayerPerceptron* FracVarAna::fneuralNet = 0;

static Int_t WeightedFit
(const Int_t n, const Double_t *x, const Double_t *y, const Double_t *w,
 Double_t *parm)
{
  Double_t sumx=0.;
  Double_t sumx2=0.;
  Double_t sumy=0.;
  Double_t sumy2=0.;
  Double_t sumxy=0.;
  Double_t sumw=0.;
  Double_t eparm[2];

  parm[0]  = 0.;
  parm[1]  = 0.;
  eparm[0] = 0.;
  eparm[1] = 0.;

  for (Int_t i=0; i<n; i++) {
    sumx += x[i]*w[i];
    sumx2 += x[i]*x[i]*w[i];
    sumy += y[i]*w[i]; 
    sumy2 += y[i]*y[i]*w[i];
    sumxy += x[i]*y[i]*w[i];
    sumw += w[i];
  }
  
  if (sumx2*sumw-sumx*sumx==0.) return 1;
  if (sumx2-sumx*sumx/sumw==0.) return 1;
  
  parm[0] = (sumy*sumx2-sumx*sumxy)/(sumx2*sumw-sumx*sumx);
  parm[1] = (sumxy-sumx*sumy/sumw)/(sumx2-sumx*sumx/sumw);
  
  eparm[0] = sumx2*(sumx2*sumw-sumx*sumx);
  eparm[1] = (sumx2-sumx*sumx/sumw);
  
  if (eparm[0]<0. || eparm[1]<0.) return 1;
  
  eparm[0] = sqrt(eparm[0])/(sumx2*sumw-sumx*sumx);
  eparm[1] = sqrt(eparm[1])/(sumx2-sumx*sumx/sumw);
  
  return 0;

}

const Float_t STP_WIDTH = 0.0412; //should be the same for the near and far detectors


FracVarAna::FracVarAna(FracVar &fv):
  fFracVar(fv),
  fDisplay(0)
{
  display = new TNtuple("display","display","tpos:z:ph:uv:core");
  static bool first = true;
                                                                                
  if(first){
    char *srt_dir = getenv("SRT_PRIVATE_CONTEXT");
    char annfile[10000];
    sprintf(annfile,"%s/NueAna/data/fvann032007.root",srt_dir);
    ifstream Test(annfile);
    if (!Test){
      srt_dir = getenv("SRT_PUBLIC_CONTEXT");
      sprintf(annfile,"%s/NueAna/data/fvann032007.root",srt_dir);
      ifstream Test_again(annfile);
      if (!Test_again){
        cout<<"Couldn't find ANN object for FracVarAna, blame Tingjun"<<endl;
        exit(0);
      }
    }
                                                                                
    static TFile *f = TFile::Open(annfile);
    fneuralNet = (TMultiLayerPerceptron*) f->Get("mlp");
    cout<<"Reading FracVar ann from : "<<annfile<<endl;
    first = false;
  }
                                                                                
  if (!fneuralNet) {
    cout<<"Couldn't find ANN object for FracVarAna, blame Tingjun"<<endl;
    exit(0);
  }
  
}

FracVarAna::~FracVarAna()
{
  if(display){
    delete display;
    display=0;
  }

}

void FracVarAna::Analyze(int evtn,  RecRecordImp<RecCandHeader> *srobj)
{
  if(srobj==0){
    return;
  }

  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }

  fFracVar.Reset();
  if (fDisplay) display->Reset();

  fFracVar.passcuts = 1;
  NtpSREvent *event = SntpHelpers::GetEvent(evtn,srobj);
  
  if ((event->ph.sigcor+event->ph.sigcor)<2e4)  fFracVar.passcuts = 0;
  if ((event->ph.sigcor+event->ph.sigcor)>1e5)  fFracVar.passcuts = 0;

  Int_t ntrks = event->ntrack;
  for (int itrk = 0; itrk<ntrks; itrk++){
    NtpSRTrack *track = SntpHelpers::GetTrack(itrk,srobj);
    if (track->plane.n>13) fFracVar.passcuts = 0;
  }

  //Int_t nshws = event->nshower;
  if (true){ //nshws){//we should have at lease one shower in the event for the rest of the analysis

    if (fDisplay){
      for (int istp = 0; istp<event->nstrip; istp++){
        Int_t index = event->stp[istp];
        NtpSRStrip*strip = SntpHelpers::GetStrip(index,srobj);
        if(!strip) continue;
//        if(!evtstp0mip){
//          MSG("FracVarAna",Msg::kError)<<"No mip strip information"<<endl;
//          continue;
//        }

        float x1 = strip->tpos;
        float x2 = strip->z;
        float x3 = (strip->ph0.sigcor+strip->ph1.sigcor)/sigcormeu;
                   //evtstp0mip[index] + evtstp1mip[index];;
        int x4 = int(strip->planeview);
        display->Fill(x1,x2,x3,x4,0);
      }
    }

//    //find the primary shower, now only require the maximum ph
//    Int_t pshw = 0; //primary shower index in the shower array
//    Float_t maxshwph = 0;
//    for (int ishw = 0; ishw<nshws; ishw++){
//      Int_t shwindex = event->shw[ishw];
//      NtpSRShower *shower = SntpHelpers::GetShower(shwindex,srobj);
//      if (shower->ph.sigcor>maxshwph){
//	pshw = shwindex;
//	maxshwph = shower->ph.sigcor;
//      }
//    }

    //Get the shower close to the track vertex or the biggest shower if
    //there is no track. Pull out from Mad.
/*
    int track_index   = -1;
    int shower_index  = -1;
    if(ReleaseType::IsBirch(release)){
      if (LoadLargestTrackFromEvent(event,srobj,track_index)){
	if(!LoadShowerAtTrackVertex(event,srobj,track_index,shower_index)){
	  LoadLargestShowerFromEvent(event,srobj,shower_index);
	}
      }
      else LoadLargestShowerFromEvent(event,srobj,shower_index);
    }
    else {//cedar ...
      if (ntrks) track_index = event->trk[0];
      //not using event->primshw, the concern is the 2GeV cut off
      if (nshws) shower_index = event->shw[0];
    }
    if (shower_index == -1) return; //no decent shower was found

    NtpSRShower *shower = SntpHelpers::GetShower(shower_index,srobj);
    if (!shower) return;     //no decent shower was found
*/

    if (event->plane.n<5) fFracVar.passcuts = 0;
    if (event->plane.n>15) fFracVar.passcuts = 0;
    //get rid of cosmics
    if (event->plane.beg>event->plane.end) return;
    //get rid of events with shw.ph>evt.ph
//    if (event->ph.mip>event->ph.mip){
//      MAXMSG("FracVarAna",Msg::kWarning,5)
//	<<"event ph("<<event->ph.mip<<"mip)>event ph("<<event->ph.mip<<"mip)"<<endl;
//      return;
//    }
      
    float vtx_ph = 0;
    //Int_t shwnstp = event->nstrip;
    Float_t plane[14];
    for (int i = 0; i<14; i++){
      plane[i] = 0;
    }

    vector<Float_t> stpph;
    vector<Float_t> stpph2; //for sorting purpose
    vector<Float_t>::iterator p;
    vector<Int_t> stppl;
    vector<Int_t> stpstp;
    vector<Float_t> stpt;
    vector<Float_t> stpz;
    vector<Int_t> planev;

    Float_t total_ph = 0;
    Int_t vtxpl = event->vtx.plane;
//    if(ReleaseType::IsCedar(release)){
//      NtpVtxFinder vtxf;
//      NtpStRecord *st = dynamic_cast<NtpStRecord *>(srobj);
//      vtxf.SetTargetEvent(evtn, st);
//      if(vtxf.FindVertex() > 0){
//	vtxpl = vtxf.VtxPlane();
//      }
//    }

    for (Int_t istp = 0; istp<event->nstrip; istp++){
      Int_t index = event->stp[istp];
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
      if(!strip) continue;
      //float charge = strip->ph0.sigcor+strip->ph1.sigcor;
 
      float charge = 0;
      //if (event->stpph1mip[istp]>0) charge += event->stpph1mip[istp];
      //if (event->stpph0mip[istp]>0) charge += event->stpph0mip[istp];

      if (evtstp0mip[index] > 0) charge += evtstp0mip[index];
      if (evtstp1mip[index] > 0) charge += evtstp1mip[index];

      if (abs(strip->plane-vtxpl)<3){
	vtx_ph+=charge;
      }
      stpph.push_back(charge);
      stpph2.push_back(charge);
      stppl.push_back(strip->plane);
      stpstp.push_back(strip->strip);
      stpt.push_back(strip->tpos);
      stpz.push_back(strip->z);
      planev.push_back(int(strip->planeview));
      total_ph += charge;
    }

    Float_t maxph = 0;

    Float_t maxph1 = -1;
    Float_t maxph2 = -1;
    int maxpl1 = -1;
    int maxpl2 = -1;

    for(unsigned int istp=0; istp<stpph.size(); istp++){//Loop over all strips
      if (stppl[istp]>=vtxpl-2&&stppl[istp]<vtxpl+12)
	plane[stppl[istp]-vtxpl+2] += stpph[istp];//fill plane energy information
      if (stpph[istp]>maxph){
	maxph = stpph[istp];
	fFracVar.shw_max = stppl[istp] - vtxpl;
      }
      if (stpph[istp]>maxph1){
	maxph1 = stpph[istp];
	maxpl1 = stppl[istp];
      }
      else if (stpph[istp]>maxph2){
	maxph2 = stpph[istp];
	maxpl2 = stppl[istp];
      }
    }
    if (maxpl1!=-1&&maxpl2!=-1) fFracVar.dis2stp = abs(maxpl1-maxpl2);
    if (fFracVar.shw_max>100) fFracVar.shw_max = -2;
    //I don't understand why this is happening. This should be a temporary solution.
    if (fFracVar.shw_max<-100) fFracVar.shw_max = -2;

    if (fFracVar.shw_max<0){
      MAXMSG("FracVarAna",Msg::kWarning,5)
	<<"fFracVar.shw_max<0"<<endl;
      fFracVar.shw_max = 0;
    }

    sort(stpph2.begin(), stpph2.end());
    
    Float_t sum_1_plane = 0,sum_2_planes = 0,sum_3_planes = 0;
    Float_t sum_4_planes = 0,sum_5_planes = 0,sum_6_planes = 0;
    Float_t sum_2_counts = 0, sum_4_counts = 0, sum_6_counts = 0;
    Float_t sum_8_counts = 0, sum_10_counts = 0, sum_12_counts = 0;
    Float_t sum_14_counts = 0, sum_16_counts = 0, sum_18_counts = 0;
    Float_t sum_20_counts = 0;

    //calculate fFracVar.fract_1_plane ...
    
    if (total_ph>0.){
      for (Int_t i = 0; i<14; i++){
	if (i<14){
	  sum_1_plane = plane[i];
	  if (sum_1_plane/total_ph > fFracVar.fract_1_plane) {
	    fFracVar.fract_1_plane = sum_1_plane/total_ph;
	  }
	}
	if (i<13){
	  sum_2_planes = plane[i]+plane[i+1];
	  if (sum_2_planes/total_ph > fFracVar.fract_2_planes) 
	    fFracVar.fract_2_planes = sum_2_planes/total_ph;
	}
	if (i<12){
	  sum_3_planes = plane[i]+plane[i+1]+plane[i+2];
	  if (sum_3_planes/total_ph > fFracVar.fract_3_planes) 
	    fFracVar.fract_3_planes = sum_3_planes/total_ph;
	}
	if (i<11){
	  sum_4_planes = plane[i]+plane[i+1]+plane[i+2]+plane[i+3];
	  if (sum_4_planes/total_ph > fFracVar.fract_4_planes) 
	    fFracVar.fract_4_planes = sum_4_planes/total_ph;
	}
	if (i<10) {
	  sum_5_planes = plane[i]+plane[i+1]+plane[i+2]+plane[i+3]+plane[i+4];
	  if (sum_5_planes/total_ph > fFracVar.fract_5_planes) 
	    fFracVar.fract_5_planes = sum_5_planes/total_ph;
	}
	if (i<9) {
	  sum_6_planes = plane[i]+plane[i+1]+plane[i+2]+plane[i+3]+plane[i+4]+plane[i+5];
	  if (sum_6_planes/total_ph > fFracVar.fract_6_planes) 
	    fFracVar.fract_6_planes = sum_6_planes/total_ph;
	}
      }
    }//if(total_ph>0)
    
    p = stpph2.end();
    while (p!=stpph2.begin()&&p>stpph2.end()-20){
      p--;
      if (p>=stpph2.end()-2) sum_2_counts += *p;
      if (p>=stpph2.end()-4) sum_4_counts += *p;
      if (p>=stpph2.end()-6) sum_6_counts += *p;
      if (p>=stpph2.end()-8) sum_8_counts += *p;
      if (p>=stpph2.end()-10) sum_10_counts += *p;
      if (p>=stpph2.end()-12) sum_12_counts += *p;
      if (p>=stpph2.end()-14) sum_14_counts += *p;
      if (p>=stpph2.end()-16) sum_16_counts += *p;
      if (p>=stpph2.end()-18) sum_18_counts += *p;
      if (p>=stpph2.end()-20) sum_20_counts += *p;
      
      //cout<<*p<<endl;
    }
    
    if (total_ph>0) {
      fFracVar.fract_2_counters = sum_2_counts/total_ph;
      fFracVar.fract_4_counters = sum_4_counts/total_ph;
      fFracVar.fract_6_counters = sum_6_counts/total_ph;
      fFracVar.fract_8_counters = sum_8_counts/total_ph;
      fFracVar.fract_10_counters = sum_10_counts/total_ph;
      fFracVar.fract_12_counters = sum_12_counts/total_ph;
      fFracVar.fract_14_counters = sum_14_counts/total_ph;
      fFracVar.fract_16_counters = sum_16_counts/total_ph;
      fFracVar.fract_18_counters = sum_18_counts/total_ph;
      fFracVar.fract_20_counters = sum_20_counts/total_ph;
    }

    vector<Double_t> upos;
    vector<Double_t> zupos;
    vector<Double_t> vpos;
    vector<Double_t> zvpos;
    vector<Double_t> ucharge;
    vector<Double_t> vcharge;
    vector<Int_t> ustrip;
    vector<Int_t> uplane;
    vector<Int_t> vstrip;
    vector<Int_t> vplane;
    
    vector<Double_t> ufit;
    vector<Double_t> zufit;
    vector<Double_t> vfit;
    vector<Double_t> zvfit;
    vector<Double_t> ufitcharge;
    vector<Double_t> vfitcharge;

    vector<Int_t> ifitu;
    vector<Int_t> ifitv;
    
    vector<Double_t> ushwstrip;
    vector<Double_t> ushwplane;
    vector<Double_t> vshwstrip;
    vector<Double_t> vshwplane;
    vector<Int_t> shwplane;
    
    Double_t total_charge = 0.;
    
    Float_t toler1 = 3*0.0412;
    Float_t toler2 = 1.5*0.0412;
    
    Double_t thr1 = 0.05; //threshold to determine shw_beg and shw_end
    
    Int_t begplane = event->plane.beg;
    Int_t endplane = event->plane.end;

    if (begplane>endplane) {
      MAXMSG("FracVarAna",Msg::kWarning,5)
	<<"event->plane.beg= "<<begplane
	<<"event->plane.end= "<<endplane<<endl;
      return;
    }
    vector<Double_t> planeph(endplane-begplane+1);
    for (int i = 0; i<endplane-begplane+1; i++){
      planeph[i] = 0;
    }
    
    for (unsigned int i = 0; i<stpph.size(); i++){
      if (planev[i] == 2){//u view
	upos.push_back(stpt[i]);
	zupos.push_back(stpz[i]);
	ustrip.push_back(stpstp[i]);
	uplane.push_back(stppl[i]);
	ucharge.push_back(stpph[i]);
	ifitu.push_back(0);
      }
      
      if (planev[i] == 3){//v view
	vpos.push_back(stpt[i]);
	zvpos.push_back(stpz[i]);
	vstrip.push_back(stpstp[i]);
	vplane.push_back(stppl[i]);
	vcharge.push_back(stpph[i]);
	ifitv.push_back(0);
      }
      
      if (stppl[i]>=begplane&&stppl[i]<=endplane){
	planeph[stppl[i]-begplane] += stpph[i];
      }
      
      total_charge += stpph[i];
      
    }
   
    // calculate the u and v charge, then asymmetry
    Float_t ucharge_tot = 0, vcharge_tot = 0;
  
    for (unsigned int i = 0; i<ucharge.size(); i++) {
      ucharge_tot += ucharge[i];
    }   

    for (unsigned int i = 0; i<vcharge.size(); i++) {
      vcharge_tot += vcharge[i];
    }

    if (ucharge_tot + vcharge_tot > 0) fFracVar.fract_asym = TMath::Abs(ucharge_tot-vcharge_tot)/(ucharge_tot+vcharge_tot); 
 
    //find event range
    Int_t shw_beg = 485;
    Int_t shw_end = 1;
    
    for (Int_t ipl = 0; ipl<endplane-begplane+1; ipl++){
      if (planeph[ipl]>total_charge*thr1/2 && ipl+begplane<shw_beg
	  && ipl<endplane-begplane && planeph[ipl+1]>total_charge*thr1)
	shw_beg = ipl + begplane;
      if (planeph[ipl]>total_charge*thr1/2&&ipl+begplane>shw_end)
	shw_end = ipl + begplane;
    }
    
    if (shw_beg == 485) shw_beg = begplane;
    if (shw_end == 1) shw_end = endplane;
    //cout<<begplane<<" "<<endplane<<endl;
    //cout<<shw_beg<<" "<<shw_end<<endl;
    
    for (unsigned int i = 0; i<ifitu.size(); i++){
      if (uplane[i]>=shw_beg&&uplane[i]<=shw_end){
	ufit.push_back(upos[i]);
	zufit.push_back(zupos[i]);
	ufitcharge.push_back(ucharge[i]);
      }
    }
    
    for (unsigned int i = 0; i<ifitv.size(); i++){
      if (vplane[i]>=shw_beg&&vplane[i]<=shw_end){
	vfit.push_back(vpos[i]);
	zvfit.push_back(zvpos[i]);
	vfitcharge.push_back(vcharge[i]);
      }
    }
    
    //perform weighted fit of u, v vs z positions
    Double_t uparm[2];
    Double_t vparm[2];
    Int_t fituok = 0;
    Int_t fitvok = 0;
    
    if (ufit.size()) {
      fituok = WeightedFit(ufit.size(),&zufit[0],&ufit[0],&ufitcharge[0],&uparm[0]);
      //cout<<"u "<<fituok<<" "<<uparm[0]<<" "<<uparm[1]<<endl;
    }
    if (vfit.size()) {
      fitvok = WeightedFit(vfit.size(),&zvfit[0],&vfit[0],&vfitcharge[0],&vparm[0]);
      //cout<<"v "<<fitvok<<" "<<vparm[0]<<" "<<vparm[1]<<endl;
    }
    
    ufit.clear();
    zufit.clear();
    ufitcharge.clear();
    vfit.clear();
    zvfit.clear();
    vfitcharge.clear();
    
    for (unsigned int i = 0; i<ifitu.size(); i++){
        if (!fituok&&TMath::Abs((upos[i]-(uparm[0]+zupos[i]*uparm[1]))*cos(atan(uparm[1])))<toler1&&uplane[i]>=shw_beg&&uplane[i]<=shw_end){
	ifitu[i] = 1;
	ushwstrip.push_back(ustrip[i]+0.5);
	ushwplane.push_back(uplane[i]+0.5);
	ufit.push_back(upos[i]);
	zufit.push_back(zupos[i]);
	ufitcharge.push_back(ucharge[i]);
      }
    }
    
    for (unsigned int i = 0; i<ifitv.size(); i++){
        if (!fitvok&&TMath::Abs((vpos[i]-(vparm[0]+zvpos[i]*vparm[1]))*cos(atan(vparm[1])))<toler1&&vplane[i]>=shw_beg&&vplane[i]<=shw_end){
	ifitv[i] = 1;
	vshwstrip.push_back(vstrip[i]+0.5);
	vshwplane.push_back(vplane[i]+0.5);
	vfit.push_back(vpos[i]);
	zvfit.push_back(zvpos[i]);
	vfitcharge.push_back(vcharge[i]);
      }
    }
    
    for(int itr = 0; itr<2; itr++){
      
      if (ufit.size()) {
	fituok = WeightedFit(ufit.size(),&zufit[0],&ufit[0],&ufitcharge[0],&uparm[0]);
	//cout<<"u "<<fituok<<" "<<uparm[0]<<" "<<uparm[1]<<endl;
      }
      if (vfit.size()) {
	fitvok = WeightedFit(vfit.size(),&zvfit[0],&vfit[0],&vfitcharge[0],&vparm[0]);
	//cout<<"v "<<fitvok<<" "<<vparm[0]<<" "<<vparm[1]<<endl;
      }
      
      ushwstrip.clear();
      ushwplane.clear();
      vshwstrip.clear();
      vshwplane.clear();
      ufit.clear();
      zufit.clear();
      ufitcharge.clear();
      vfit.clear();
      zvfit.clear();
      vfitcharge.clear();
      
      
      for (unsigned int i = 0; i<ifitu.size(); i++){
	ifitu[i] = 0;
          if (!fituok&&TMath::Abs((upos[i]-(uparm[0]+zupos[i]*uparm[1]))*cos(atan(uparm[1])))<toler2&&uplane[i]>=shw_beg&&uplane[i]<=shw_end){
	  ifitu[i] = 1;
	  ushwstrip.push_back(ustrip[i]+0.5);
	  ushwplane.push_back(uplane[i]+0.5);
	  ufit.push_back(upos[i]);
	  zufit.push_back(zupos[i]);
	  ufitcharge.push_back(ucharge[i]);
	  if (itr==1) {
	    shwplane.push_back(uplane[i]);
	    display->Fill(upos[i],zupos[i],100000,2,1);
	  }
	}
      }
      
      for (unsigned int i = 0; i<ifitv.size(); i++){
	ifitv[i] = 0;
          if (!fitvok&&TMath::Abs((vpos[i]-(vparm[0]+zvpos[i]*vparm[1]))*cos(atan(vparm[1])))<toler2&&vplane[i]>=shw_beg&&vplane[i]<=shw_end){
	  ifitv[i] = 1;
	  vshwstrip.push_back(vstrip[i]+0.5);
	  vshwplane.push_back(vplane[i]+0.5);
	  vfit.push_back(vpos[i]);
	  zvfit.push_back(zvpos[i]);
	  vfitcharge.push_back(vcharge[i]);
	  if (itr==1) {
	    shwplane.push_back(vplane[i]);
	    display->Fill(vpos[i],zvpos[i],100000,3,1);
	  }
	}
      }
    }
    
    Double_t shwcoreph = 0;
    fFracVar.shw_nstp = 0;
    for (unsigned int i = 0; i<ufitcharge.size(); i++){
      shwcoreph += ufitcharge[i];
      fFracVar.shw_nstp ++;
    }
    for (unsigned int i = 0; i<vfitcharge.size(); i++){
      shwcoreph += vfitcharge[i];
      fFracVar.shw_nstp ++;
    }
    if (total_ph) {
      fFracVar.fract_road = shwcoreph/total_ph;
      fFracVar.vtxph = vtx_ph/total_ph;
    }
    Int_t shw_begpl = 1000;
    Int_t shw_endpl = 0;
    for (unsigned int i = 0; i<shwplane.size(); i++){
      if (shw_begpl>shwplane[i]) shw_begpl = shwplane[i];
      if (shw_endpl<shwplane[i]) shw_endpl = shwplane[i];
   }
    if (shw_endpl>=shw_begpl) fFracVar.shw_npl = shw_endpl - shw_begpl + 1;

    double slope = -1;
    if (!fituok&&!fitvok&&TMath::Abs(uparm[1])<100&&TMath::Abs(vparm[1])<100){
      slope = sqrt(pow(uparm[1],2)+pow(vparm[1],2));
    }
    fFracVar.shw_slp = slope;
    //calculate ANN pid
    mlpANN ann;
    Double_t params[14];
    params[0] = fFracVar.fract_1_plane;
    params[1] = fFracVar.fract_2_planes;
    params[2] = fFracVar.fract_3_planes;
    params[3] = fFracVar.fract_4_planes;
    params[4] = fFracVar.fract_5_planes;
    params[5] = fFracVar.fract_6_planes;
    params[6] = fFracVar.fract_2_counters;
    params[7] = fFracVar.fract_4_counters;
    params[8] = fFracVar.fract_6_counters;
    params[9] = fFracVar.fract_8_counters;
    params[10] = fFracVar.fract_10_counters;
    params[11] = fFracVar.fract_12_counters;
    params[12] = fFracVar.fract_road;
    params[13] = fFracVar.shw_max;
    
    fFracVar.pid = ann.value(0,params[0],params[1],params[2],params[3],params[4],params[5],params[6],params[7],params[8],params[9],params[10],params[11],params[12],params[13]);
    fFracVar.pid1 = fneuralNet->Evaluate(0,params);
    
  }//if (nshws)
  else {
    fFracVar.passcuts = 0;
  }
}


//......................................................................

bool FracVarAna::LoadLargestTrackFromEvent(NtpSREvent* event,
				      RecRecordImp<RecCandHeader>* record,
				      int& trkidx)
{
  bool status = false;
  float longest_trk = 0;
  for (int i=0; i<event->ntrack; i++){
    NtpSRTrack *track = SntpHelpers::GetTrack(event->trk[i],record);
    if (!track) continue;
    int trklen = abs(track->plane.end-track->plane.beg)+1;
    if (trklen>longest_trk){
      trkidx = event->trk[i];
      longest_trk = trklen;
    }
  }
  if (longest_trk>0){
    status = true;
  }
  return status;
}

//.....................................................................

bool FracVarAna::LoadShowerAtTrackVertex(NtpSREvent* event,
				    RecRecordImp<RecCandHeader>* record,
				    int  trkidx,
				    int& shwidx)
{
  bool status = false;
  NtpSRTrack *track = SntpHelpers::GetTrack(trkidx,record);
  if (!track) return false;
  float closest_ph = 0.;
  float closest_vtx = 1000.;
  const float close_enough = 7*0.06; //about 1 interaction length
  for (int i=0; i<event->nshower; i++){
    NtpSRShower *shower = SntpHelpers::GetShower(event->shw[i],record);
    if (!shower) continue;
    float vtx_sep = TMath::Abs(shower->vtx.z-track->vtx.z);

    //if this is the closest shower
    if (vtx_sep<closest_vtx){
      shwidx = event->shw[i];
      closest_vtx = vtx_sep;
      if (shower->shwph.linCCgev>0) closest_ph = shower->shwph.linCCgev;
      else closest_ph = shower->ph.gev;
    }
    // if this isn't the closest shower but it's within close_enough m of the
    // track vertex and has a larger pulseheight
    else if(vtx_sep<close_enough && shower->ph.gev>closest_ph){
      shwidx = event->shw[i];
      closest_vtx = vtx_sep;

      if (shower->shwph.linCCgev>0) closest_ph = shower->shwph.linCCgev;
      else closest_ph = shower->ph.gev;
    }
  }
  if (closest_vtx<close_enough){
    status = true;
  }
  else shwidx = -1;
  return status;
}

//.........................................................................

bool FracVarAna::LoadLargestShowerFromEvent(NtpSREvent* event,
				       RecRecordImp<RecCandHeader>* record,
				       int& shwidx)
{
  bool status = false;
  float largest_ph = 0;
  for (int i=0; i<event->nshower; i++){
    NtpSRShower *shower = SntpHelpers::GetShower(event->shw[i],record);
    if (!shower) continue;
    if (shower->ph.gev>largest_ph){
      shwidx = event->shw[i];
      largest_ph = shower->ph.gev;
    }
  }
  if (largest_ph>0){
    status = true;
  }
  return status;
}

//........................................................................

void FracVarAna::Draw(TPad *pad){
  pad->cd(1);
  //TLatex *t1 = new TLatex(0.5,0.5,Form("%.3f",fFracVar.fract_road));
  //t1->Draw();
  display->Draw("tpos:z","ph*(uv==2&&!core)","colz");
  display->Draw("tpos:z","ph*(uv==2&&core)","box same");
  pad->cd(2);
  display->Draw("tpos:z","ph*(uv==3&&!core)","colz");
  display->Draw("tpos:z","ph*(uv==3&&core)","box same");
  pad->Modified();
}

void FracVarAna::Print(TPad */*pad*/){
}
			  
