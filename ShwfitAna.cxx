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
#include "MessageService/MsgService.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/ShwfitAna.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"

CVSID("$Id: ShwfitAna.cxx,v 1.39 2011/10/07 20:04:33 wingmc Exp $");

static Double_t shwfunc(Double_t *x, Double_t *par)
{

  // function to be fitted for em showers par[0]=a par[1]=E0
  Float_t R = 1.46676;
  Float_t xx=R*x[0];
//  cout<<"In shwfunc "<<xx<<" par[0] "<<par[0]<<" par[1] "<<par[1]<<" par[2] "<<par[2]<<endl;

  Double_t lnf = TMath::Log(R*par[2]*par[1])+(par[0]-1)*TMath::Log(xx*par[1])-
     par[1]*xx-TMath::LnGamma(par[0]);

  Double_t f = exp(lnf);

//  Double_t f = 1.46676*par[2]*TMath::Power(par[1],par[0])*
//               TMath::Power(xx,par[0]-1)*
//               TMath::Exp(-par[1]*(xx))/TMath::Gamma(par[0]);
  return f;
}

static Double_t hadfunc(Double_t *x, Double_t *par)
{
  // function to be fitted for em showers par[0]=a par[1]=b par[2]=E0
  // added hadronic part, new params par[3]=c par[4]=d par[5]=fzero
  Float_t xx=x[0];
  Double_t f1 = par[2]*TMath::Power(par[1],par[0])*
                TMath::Power(xx*1.46676,par[0]-1)*
                TMath::Exp(-par[1]*(xx*1.46676))/TMath::Gamma(par[0]);
  Double_t f2 = par[2]*TMath::Power(par[4],par[3])*
                TMath::Power(xx*0.1694,par[3]-1)*
                TMath::Exp(-par[4]*(xx*0.1694))/TMath::Gamma(par[3]);
  Double_t f=(par[5]*f1+(1-par[5])*0.1154*f2)*1.46676;
  return f;
}

//----------------------------------------------------------------------------------
//DDC Mawell fit Code
static Double_t maxwell(Double_t *x, Double_t *par)
{
  //Fit to Maxwell function N(E*b) x^2 exp{-b*x^2}
  //par[0] = b, par[1] = E = Energy
  Float_t xx = x[0];
                                                                                                   
  Double_t f = par[1]*(4/TMath::Sqrt(TMath::Pi()))*TMath::Power(par[0],(3/2))*TMath::Power(xx,2)*TMath::Exp(-par[0]*TMath::Power(xx,2));
                                                                                                   
  return f;
}
static Double_t maxwell3(Double_t *x, Double_t *par)
{
  //Fit to Maxwell-like function N(E*b) x^3 exp{-b*x^2}
  //par[0] = b, par[1] = E = Energy
  Float_t xx = x[0];
                                                                                                   
  Double_t f = par[1]*2*TMath::Power(par[0],2)*TMath::Power(xx,3)*TMath::Exp(-par[0]*TMath::Power(xx,2));
                                                                                                   
  return f;
}
//--------------------------------------------------------------------------------------

ShwfitAna::ShwfitAna(Shwfit &sf):
   fShwfit(sf)
{}

ShwfitAna::~ShwfitAna()
{}

//void ShwfitAna::Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord */*mc*/, NtpTHRecord */*th*/)

void ShwfitAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
  if(srobj==0){
    return;
  }
  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
     ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
    return;
  }

  MSG("ShwfitAna",Msg::kDebug)<<"In ShwfitAna::Analyze"<<endl;
  MSG("ShwfitAna",Msg::kDebug)<<"On Snarl "<<srobj->GetHeader().GetSnarl()
			      <<" event "<<evtn<<endl;

  Reset(srobj->GetHeader().GetSnarl(),evtn);

  NtpSREvent *event = SntpHelpers::GetEvent(evtn,srobj);

  if(!event){
      MSG("ShwfitAna",Msg::kError)<<"Couldn't get event "<<evtn
				   <<" from Snarl "<<srobj->GetHeader().GetSnarl()<<endl;
      return;
   }

  Int_t vtxPlane = event->vtx.plane;
  Float_t vtxU = event->vtx.u;
  Float_t vtxV = event->vtx.v;

  if(ReleaseType::IsCedar(release)){
    NtpStRecord* st = dynamic_cast<NtpStRecord *>(srobj);
    NtpVtxFinder vtxf(evtn, st);
    if(vtxf.FindVertex() > 0){
       vtxPlane = vtxf.VtxPlane();
       vtxU = vtxf.VtxU();
       vtxV = vtxf.VtxV();
    }
  }

   //loop over strips to fill histograms
  fShwfit.hiPhStripCountM4=0;
  fShwfit.hiPhPlaneCountM4=0;
  fShwfit.hiPhStripCountM2=0;
  fShwfit.hiPhPlaneCountM2=0;
  fShwfit.hiPhStripCount=0;
  fShwfit.hiPhPlaneCount=0;
  fShwfit.hiPhStripCountP2=0;
  fShwfit.hiPhPlaneCountP2=0;
  fShwfit.hiPhStripCountP4=0;
  fShwfit.hiPhPlaneCountP4=0;

  // new contPlaneCount - Minerba
  fShwfit.contPlaneCount=0;

  float VertexEnergy = 0.0;
  Float_t vtxEnergy0 = 0.0; 
  Float_t vtxEnergy1 = 0.0; 
  Float_t vtxEnergy2 = 0.0;


   doSlopes(event,srobj);


  for(int i=0;i<event->nstrip;i++){
      Int_t index = SntpHelpers::GetStripIndex(i,event);
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
      if(!strip){
	 continue;
      }
      if(!evtstp0mip){
        MSG("ShwfitAna",Msg::kError)<<"No mip strip information"<<endl;
        continue;
      }

      Float_t deltaplanes = strip->plane-event->vtx.plane;
      Float_t stripPh = evtstp0mip[index] + evtstp1mip[index];
      double strippe = strip->ph0.pe+strip->ph1.pe;

      if(deltaplanes == 0 || deltaplanes == 1) VertexEnergy += stripPh;
      if(deltaplanes == 0) vtxEnergy0 += stripPh;
      if(deltaplanes == 1) vtxEnergy1 += stripPh;
      if(deltaplanes == 2) vtxEnergy2 += stripPh;

      const Float_t STRIPWIDTH=0.041;//Munits::meters
      MSG("ShwfitAna",Msg::kDebug)<< "plane " << deltaplanes << " stripPh " << stripPh <<endl;
      //fill longitudnal energy deposition histogram
      if(deltaplanes>=0.){
	 fShwfit.lenepl->Fill(deltaplanes-0.5,stripPh);
         fShwfit.ph_hist->Fill(stripPh);        
 
	 if(sfDPlaneCut>0&&deltaplanes<sfDPlaneCut-4&&stripPh>sfPhStripCut) fShwfit.hiPhStripCountM4++;
	 if(sfDPlaneCut>0&&deltaplanes<sfDPlaneCut-2&&stripPh>sfPhStripCut) fShwfit.hiPhStripCountM2++;
	 if(sfDPlaneCut>0&&deltaplanes<sfDPlaneCut&&stripPh>sfPhStripCut) fShwfit.hiPhStripCount++;
	 if(sfDPlaneCut>0&&deltaplanes<sfDPlaneCut+2&&stripPh>sfPhStripCut) fShwfit.hiPhStripCountP2++;
	 if(sfDPlaneCut>0&&deltaplanes<sfDPlaneCut+4&&stripPh>sfPhStripCut) fShwfit.hiPhStripCountP4++;

      }

      //fill transverse energy deposition histograms
      //u view
      if(strip->planeview==PlaneView::kU){
         double dist = (strip->tpos-event->vtx.u)/STRIPWIDTH;
  
	 fShwfit.tenestu->Fill(dist,stripPh);
         if(dist < 9.0 && strippe > 2.0){
               fShwfit.tenestu_9s_2pe->Fill(dist,stripPh);
               fShwfit.tenestu_9s_2pe_dw->Fill(dist,stripPh*stripPh);
         }
      }

      //v view
      else if(strip->planeview==PlaneView::kV){
         double dist = (strip->tpos-event->vtx.v)/STRIPWIDTH;
                                                                                                             
         fShwfit.tenestv->Fill(dist,stripPh);
         if(dist < 9.0 && strippe > 2.0){
               fShwfit.tenestv_9s_2pe->Fill(dist,stripPh);
               fShwfit.tenestv_9s_2pe_dw->Fill(dist,stripPh*stripPh);
         }
      }
      //unknown view
      else{
	 MSG("ShwfitAna",Msg::kError)<<"Don't know what to do with a PlaneView "
				      <<strip->planeview<<" skipping"<<endl;
	 continue;
      }     
   }


  fShwfit.vtxEnergy = VertexEnergy;
  fShwfit.energyPlane0 =  vtxEnergy0;
  fShwfit.energyPlane1 =  vtxEnergy1;
  fShwfit.energyPlane2 =  vtxEnergy2;


  int planebins=fShwfit.lenepl->GetNbinsX();

  if(sfDPlaneCut<planebins) planebins=sfDPlaneCut;
  for(int i=1;i<=planebins+4;i++){
    if(fShwfit.lenepl->GetBinContent(i)>sfPhPlaneCut&&i<=planebins-4) fShwfit.hiPhPlaneCountM4++;
    if(fShwfit.lenepl->GetBinContent(i)>sfPhPlaneCut&&i<=planebins-2) fShwfit.hiPhPlaneCountM2++;
    if(fShwfit.lenepl->GetBinContent(i)>sfPhPlaneCut&&i<=planebins) fShwfit.hiPhPlaneCount++;
    if(fShwfit.lenepl->GetBinContent(i)>sfPhPlaneCut&&i<=planebins+2) fShwfit.hiPhPlaneCountP2++;
    if(fShwfit.lenepl->GetBinContent(i)>sfPhPlaneCut&&i<=planebins+4) fShwfit.hiPhPlaneCountP4++;
  }

  //  contPlaneCount 

  if(sfContPhPlaneCut>0){
    bool PlaneCut=true;
    fShwfit.contPlaneCount=0;
    if(fShwfit.lenepl->GetBinContent(1)>sfContPhPlaneCut){
      fShwfit.contPlaneCount++;
    }
    for(int i=2;i<=30; i++){
      if(fShwfit.lenepl->GetBinContent(i)>sfContPhPlaneCut&&PlaneCut){
	fShwfit.contPlaneCount++;
      }else {
	PlaneCut=false;
      }   
    }
  }


  bool PlaneCut015=true; // ~0.75PE
  bool PlaneCut030=true; // ~1.50PE
  bool PlaneCut050=true; // ~2.00PE
  bool PlaneCut075=true; // ~3.00PE
  bool PlaneCut100=true; // ~4.00PE
  bool PlaneCut200=true; // ~8.00PE

  fShwfit.contPlaneCount015=0;
  fShwfit.contPlaneCount030=0;
  fShwfit.contPlaneCount050=0;
  fShwfit.contPlaneCount075=0;
  fShwfit.contPlaneCount100=0;
  fShwfit.contPlaneCount200=0;

  if(fShwfit.lenepl->GetBinContent(1)>0.15) fShwfit.contPlaneCount015++;
  if(fShwfit.lenepl->GetBinContent(1)>0.30) fShwfit.contPlaneCount030++;
  if(fShwfit.lenepl->GetBinContent(1)>0.50) fShwfit.contPlaneCount050++;
  if(fShwfit.lenepl->GetBinContent(1)>0.75) fShwfit.contPlaneCount075++;
  if(fShwfit.lenepl->GetBinContent(1)>1.00) fShwfit.contPlaneCount100++;
  if(fShwfit.lenepl->GetBinContent(1)>2.00) fShwfit.contPlaneCount200++;
  
  for(int i=2;i<=30; i++){
    if(fShwfit.lenepl->GetBinContent(i)>0.15&&PlaneCut015){
      fShwfit.contPlaneCount015++;
    }else {
      PlaneCut015=false;
    }   
    if(fShwfit.lenepl->GetBinContent(i)>0.30&&PlaneCut030){
	fShwfit.contPlaneCount030++;
    }else {
	PlaneCut030=false;
    }   
    if(fShwfit.lenepl->GetBinContent(i)>0.50&&PlaneCut050){
      fShwfit.contPlaneCount050++;
    }else {
      PlaneCut050=false;
    }   
    if(fShwfit.lenepl->GetBinContent(i)>0.75&&PlaneCut075){
      fShwfit.contPlaneCount075++;
    }else {
      PlaneCut075=false;
    }   
    if(fShwfit.lenepl->GetBinContent(i)>1.00&&PlaneCut100){
      fShwfit.contPlaneCount100++;
    }else {
      PlaneCut100=false;
    }   
    if(fShwfit.lenepl->GetBinContent(i)>2.00&&PlaneCut200){
      fShwfit.contPlaneCount200++;
    }else {
      PlaneCut200=false;
    }   
  }
  

  MSG("ShwfitAna",Msg::kDebug)<<"StripCount "<< fShwfit.hiPhStripCount 
			      << " PlaneCount " << fShwfit.hiPhPlaneCount 
                              << " contPlaneCount "<<fShwfit.contPlaneCount<< endl;
  //////////////////// 


 MSG("ShwfitAna",Msg::kDebug)<<"In ShwfitAna, trying to fit"<<endl;
   //fit EM shower
   Int_t trk_plane_num=0;
   Int_t ntrack=event->ntrack;
   if(ntrack>0){  
      Int_t trkNum=SntpHelpers::GetTrackIndex(0,event);
      NtpSRTrack *track = SntpHelpers::GetTrack(trkNum,srobj);
      trk_plane_num=track->plane.n;
   }

   if(event->plane.n > 4){
      FitLShower(event->ph.mip);
      FitLShower_Dan(event->ph.mip);
   }
   MSG("ShwfitAna",Msg::kDebug)<<"Fit shower "<<fShwfit.par_a<<" "
			       <<fShwfit.par_b<<" "<<fShwfit.par_e0<<endl;
   
   //fill transverse variables (the following variables filled inside
   //   function and then stored
   FitTShower(event->ph.mip);

   TransVar(fShwfit.tenestu, PlaneView::kU);
   fShwfit.u_asym_peak=asym_peak;
   fShwfit.u_asym_vert=asym_vert;
   fShwfit.u_molrad_peak=molrad_peak;
   fShwfit.u_molrad_vert=molrad_vert;
   fShwfit.u_mean=mean;
   fShwfit.u_rms=rms;
   fShwfit.u_skew=skew;
   fShwfit.u_kurt=kurt;
  
   TransVar(fShwfit.tenestv, PlaneView::kV);
   fShwfit.v_asym_peak=asym_peak;
   fShwfit.v_asym_vert=asym_vert;
   fShwfit.v_molrad_peak=molrad_peak;
   fShwfit.v_molrad_vert=molrad_vert;
   fShwfit.v_mean=mean;
   fShwfit.v_rms=rms;
   fShwfit.v_skew=skew;
   fShwfit.v_kurt=kurt;

   //fill uv variables

   fShwfit.uv_asym_peak = BuildUVVar(fShwfit.u_asym_peak, fShwfit.v_asym_peak);
   fShwfit.uv_asym_vert = BuildUVVar(fShwfit.u_asym_vert, fShwfit.v_asym_vert);
   fShwfit.uv_molrad_peak = BuildUVVar(fShwfit.u_molrad_peak, fShwfit.v_molrad_peak);
   fShwfit.uv_molrad_vert = BuildUVVar(fShwfit.u_molrad_vert, fShwfit.v_molrad_vert);
   fShwfit.uv_mean = BuildUVVar(fShwfit.u_mean, fShwfit.v_mean);
   fShwfit.uv_rms = BuildUVVar(fShwfit.u_rms, fShwfit.v_rms);
   fShwfit.uv_skew = BuildUVVar(fShwfit.u_skew, fShwfit.v_skew);
   fShwfit.uv_kurt = BuildUVVar(fShwfit.u_kurt, fShwfit.v_kurt);

   if(fShwfit.tenestu->Integral()>0&&fShwfit.tenestv->Integral()>0){
      fShwfit.uv_ratio=fShwfit.tenestu->Integral()/fShwfit.tenestv->Integral();
      if(fShwfit.uv_ratio>100.){
	 fShwfit.uv_ratio= ANtpDefVal::kFloat;
      }
   }

   TransVar(fShwfit.tenestu_9s_2pe, PlaneView::kU);
   fShwfit.u_asym_peak_9s_2pe=asym_peak;
   fShwfit.u_asym_vert_9s_2pe=asym_vert;
   fShwfit.u_molrad_peak_9s_2pe=molrad_peak;
   fShwfit.u_molrad_vert_9s_2pe=molrad_vert;
   fShwfit.u_mean_9s_2pe=mean;
   fShwfit.u_rms_9s_2pe=rms;
   fShwfit.u_skew_9s_2pe=skew;
   fShwfit.u_kurt_9s_2pe=kurt;

   TransVar(fShwfit.tenestv_9s_2pe, PlaneView::kV);
   fShwfit.v_asym_peak_9s_2pe=asym_peak;
   fShwfit.v_asym_vert_9s_2pe=asym_vert;
   fShwfit.v_molrad_peak_9s_2pe=molrad_peak;
   fShwfit.v_molrad_vert_9s_2pe=molrad_vert;
   fShwfit.v_mean_9s_2pe=mean;
   fShwfit.v_rms_9s_2pe=rms;
   fShwfit.v_skew_9s_2pe=skew;
   fShwfit.v_kurt_9s_2pe=kurt;
     
   fShwfit.uv_asym_peak_9s_2pe = BuildUVVar(fShwfit.u_asym_peak_9s_2pe, fShwfit.v_asym_peak_9s_2pe);
   fShwfit.uv_asym_vert_9s_2pe = BuildUVVar(fShwfit.u_asym_vert_9s_2pe, fShwfit.v_asym_vert_9s_2pe);
   fShwfit.uv_molrad_peak_9s_2pe = BuildUVVar(fShwfit.u_molrad_peak_9s_2pe, fShwfit.v_molrad_peak_9s_2pe);
   fShwfit.uv_molrad_vert_9s_2pe = BuildUVVar(fShwfit.u_molrad_vert_9s_2pe, fShwfit.v_molrad_vert_9s_2pe);
   fShwfit.uv_mean_9s_2pe = BuildUVVar(fShwfit.u_mean_9s_2pe, fShwfit.v_mean_9s_2pe);
   fShwfit.uv_rms_9s_2pe = BuildUVVar(fShwfit.u_rms_9s_2pe, fShwfit.v_rms_9s_2pe);
   fShwfit.uv_skew_9s_2pe = BuildUVVar(fShwfit.u_skew_9s_2pe, fShwfit.v_skew_9s_2pe);
   fShwfit.uv_kurt_9s_2pe = BuildUVVar(fShwfit.u_kurt_9s_2pe, fShwfit.v_kurt_9s_2pe);
   if(fShwfit.tenestu_9s_2pe->Integral()>0&&fShwfit.tenestv_9s_2pe->Integral()>0){
     fShwfit.uv_ratio_9s_2pe=fShwfit.tenestu_9s_2pe->Integral()/fShwfit.tenestv_9s_2pe->Integral();
     if(fShwfit.uv_ratio_9s_2pe>100.){
         fShwfit.uv_ratio_9s_2pe= ANtpDefVal::kFloat;
      }
   }

  TransVar(fShwfit.tenestu_9s_2pe_dw, PlaneView::kU);
   fShwfit.u_asym_peak_9s_2pe_dw=asym_peak;
   fShwfit.u_asym_vert_9s_2pe_dw=asym_vert;
   fShwfit.u_molrad_peak_9s_2pe_dw=molrad_peak;
   fShwfit.u_molrad_vert_9s_2pe_dw=molrad_vert;
   fShwfit.u_mean_9s_2pe_dw=mean;
   fShwfit.u_rms_9s_2pe_dw=rms;
   fShwfit.u_skew_9s_2pe_dw=skew;
   fShwfit.u_kurt_9s_2pe_dw=kurt;
                                                                                
                                                                                
   TransVar(fShwfit.tenestv_9s_2pe_dw, PlaneView::kV);
   fShwfit.v_asym_peak_9s_2pe_dw=asym_peak;
   fShwfit.v_asym_vert_9s_2pe_dw=asym_vert;
   fShwfit.v_molrad_peak_9s_2pe_dw=molrad_peak;
   fShwfit.v_molrad_vert_9s_2pe_dw=molrad_vert;
   fShwfit.v_mean_9s_2pe_dw=mean;
   fShwfit.v_rms_9s_2pe_dw=rms;
   fShwfit.v_skew_9s_2pe_dw=skew;
   fShwfit.v_kurt_9s_2pe_dw=kurt;
                                                                                
                                                                                
   fShwfit.uv_asym_peak_9s_2pe_dw = BuildUVVar(fShwfit.u_asym_peak_9s_2pe_dw, fShwfit.v_asym_peak_9s_2pe_dw);
   fShwfit.uv_asym_vert_9s_2pe_dw = BuildUVVar(fShwfit.u_asym_vert_9s_2pe_dw, fShwfit.v_asym_vert_9s_2pe_dw);
   fShwfit.uv_molrad_peak_9s_2pe_dw = BuildUVVar(fShwfit.u_molrad_peak_9s_2pe_dw, fShwfit.v_molrad_peak_9s_2pe_dw);
   fShwfit.uv_molrad_vert_9s_2pe_dw = BuildUVVar(fShwfit.u_molrad_vert_9s_2pe_dw, fShwfit.v_molrad_vert_9s_2pe_dw);
   fShwfit.uv_mean_9s_2pe_dw = BuildUVVar(fShwfit.u_mean_9s_2pe_dw, fShwfit.v_mean_9s_2pe_dw);
   fShwfit.uv_rms_9s_2pe_dw = BuildUVVar(fShwfit.u_rms_9s_2pe_dw, fShwfit.v_rms_9s_2pe_dw);
   fShwfit.uv_skew_9s_2pe_dw = BuildUVVar(fShwfit.u_skew_9s_2pe_dw, fShwfit.v_skew_9s_2pe_dw);
   fShwfit.uv_kurt_9s_2pe_dw = BuildUVVar(fShwfit.u_kurt_9s_2pe_dw, fShwfit.v_kurt_9s_2pe_dw);
                                                                                
                                                                                
  if(fShwfit.tenestu_9s_2pe_dw->Integral()>0&&fShwfit.tenestv_9s_2pe_dw->Integral()>0){
     fShwfit.uv_ratio_9s_2pe_dw=fShwfit.tenestu_9s_2pe_dw->Integral()/fShwfit.tenestv_9s_2pe_dw->Integral();
     if(fShwfit.uv_ratio_9s_2pe_dw>100.){
         fShwfit.uv_ratio_9s_2pe_dw= ANtpDefVal::kFloat;
     }
  }
}


void ShwfitAna::FitLShower(Float_t pulseheight)
{
   fShwfit.efit->SetParameters(3.,0.5,pulseheight);
   fShwfit.lenepl->Fit(fShwfit.efit,"RLLQ0+");
   //   fShwfit.lenepl->Print("all");
   //   fShwfit.efit->Print("all");

   MSG("ShwfitAna",Msg::kDebug)<<" STATUS: "<<gMinuit->fCstatu
				<<" "<<fShwfit.efit->GetParameter(0)
				<<" "<<fShwfit.efit->GetNDF()<<" "<<pulseheight<<endl;

   string fitstatus = (string)(gMinuit->fCstatu);
   MSG("ShwfitAna",Msg::kDebug)<<" STATUS: "<<fitstatus
				<<", "<<fShwfit.efit->GetParameter(0)
				<<", "<<fShwfit.efit->GetNDF()<<" "<<pulseheight<<endl;

   if(fitstatus=="CONVERGED "&&fShwfit.efit->GetParameter(0)<29.9&&
      fShwfit.efit->GetNDF()>0&&pulseheight>0){
      MSG("ShwfitAna",Msg::kDebug)<<" filling vars"<<endl;

      fShwfit.par_a=fShwfit.efit->GetParameter(0);
      fShwfit.par_b=fShwfit.efit->GetParameter(1);
      fShwfit.par_e0=fShwfit.efit->GetParameter(2);
      fShwfit.chisq=fShwfit.efit->GetChisquare();
      fShwfit.shwmax=1.*GetMaximumX(fShwfit.efit);
      //      fShwfit.shwmax=1.*fShwfit.lenepl->GetFunction(fShwfit.efit->GetName())->GetMaximumX();


      fShwfit.shwmaxplane=(int)(fShwfit.lenepl->
				GetBinContent(fShwfit.lenepl->FindBin(fShwfit.shwmax)));
      fShwfit.conv=1;
      if(fShwfit.efit->GetNDF()>0){
	 fShwfit.chisq_ndf=fShwfit.chisq/fShwfit.efit->GetNDF();
      }
      if(pulseheight!=0){
	 fShwfit.e0_pe_ratio=fShwfit.par_e0/pulseheight;
      }
      fShwfit.caldet_comp=ANtpDefVal::kFloat;
      fShwfit.max_pe_plane=fShwfit.lenepl->GetMaximumBin();
      fShwfit.shwmaxplane_diff=fShwfit.shwmaxplane-fShwfit.max_pe_plane;
   }

   return;
}

void ShwfitAna::FitLShower_Dan(Float_t pulseheight)
{
   //DDC Mawell fit Code
   fShwfit.efit_maxwell->SetParameters(1.0,pulseheight);
   fShwfit.lenepl->Fit(fShwfit.efit_maxwell,"RLQ0+");
   fShwfit.Beta_Maxwell = fShwfit.efit_maxwell->GetParameter(0);
   int np = 1000;
   double *x = new double[np];
   double *w = new double[np];
   fShwfit.efit_maxwell->CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
   fShwfit.Energy_Maxwell = fShwfit.efit_maxwell->IntegralFast(np,x,w,0,1000);
   fShwfit.chisq_Maxwell=fShwfit.efit_maxwell->GetChisquare();
   fShwfit.ndf_Maxwell=fShwfit.efit_maxwell->GetNDF();
   
   fShwfit.efit_maxwell3->SetParameters(1.0,pulseheight);
   fShwfit.lenepl->Fit(fShwfit.efit_maxwell3,"RLQ0+");
   fShwfit.Beta_Maxwell3 = fShwfit.efit_maxwell3->GetParameter(0);
   int np3 = 1000;
   double *x3 = new double[np3];
   double *w3 = new double[np3];
   fShwfit.efit_maxwell3->CalcGaussLegendreSamplingPoints(np3,x3,w3,1e-15);
   fShwfit.Energy_Maxwell3 = fShwfit.efit_maxwell3->IntegralFast(np3,x3,w3,0,1000);
   fShwfit.chisq_Maxwell3=fShwfit.efit_maxwell3->GetChisquare();
   fShwfit.ndf_Maxwell3=fShwfit.efit_maxwell3->GetNDF();
 
                                                                                                  
   float E_half_num = 0;
   float E_half_den = 0;
   float n_half_num = 0;
   float n_half_den = 0;
   float E_2_num = 0;
   float E_2_den = 0;
   float n_2_num = 0;
   float n_2_den = 0;
   float E_split_val = 0;
   float E_split_point = 0;
   float n_split_val = 0;
   float n_split_point = 0;
   float n_integral = 0;
   float E_integral = 0;
   float ph_hist_length = 0;
   float ph_hist_counter = 0;
   for (int q = 1; q<=fShwfit.ph_hist->GetNbinsX(); q++){
     n_integral += fShwfit.ph_hist->GetBinContent(q);
     E_integral += fShwfit.ph_hist->GetBinContent(q)*(fShwfit.ph_hist->GetBinLowEdge(q)+(1.0/2.0)*fShwfit.ph_hist->GetBinWidth(q));
     ph_hist_counter++;
     if(fShwfit.ph_hist->GetBinContent(q)>0){
       ph_hist_length += ph_hist_counter;
       ph_hist_counter = 0;
     }
   }
   for (int q = 1; q<=fShwfit.ph_hist->GetNbinsX(); q++){
     if (q<=2){
       n_2_num += fShwfit.ph_hist->GetBinContent(q);
       E_2_num += fShwfit.ph_hist->GetBinContent(q)*(fShwfit.ph_hist->GetBinLowEdge(q)+(1.0/2.0)*fShwfit.ph_hist->GetBinWidth(q));
     }else{
       n_2_den += fShwfit.ph_hist->GetBinContent(q);
       E_2_den += fShwfit.ph_hist->GetBinContent(q)*(fShwfit.ph_hist->GetBinLowEdge(q)+(1.0/2.0)*fShwfit.ph_hist->GetBinWidth(q));
     }
     if (q<=ph_hist_length/2){
       n_half_num += fShwfit.ph_hist->GetBinContent(q);
       E_half_num += fShwfit.ph_hist->GetBinContent(q)*(fShwfit.ph_hist->GetBinLowEdge(q)+(1.0/2.0)*fShwfit.ph_hist->GetBinWidth(q));
     }else{
       n_half_den += fShwfit.ph_hist->GetBinContent(q);
       E_half_den += fShwfit.ph_hist->GetBinContent(q)*(fShwfit.ph_hist->GetBinLowEdge(q)+(1.0/2.0)*fShwfit.ph_hist->GetBinWidth(q));
     }
     if(n_split_val<(1.0/2.0)*n_integral) n_split_val += fShwfit.ph_hist->GetBinContent(q);
     else if (n_split_point==0) n_split_point = (fShwfit.ph_hist->GetBinLowEdge(q)+fShwfit.ph_hist->GetBinWidth(q));
     if(E_split_val<(1.0/2.0)*E_integral) E_split_val += fShwfit.ph_hist->GetBinContent(q)*(fShwfit.ph_hist->GetBinLowEdge(q)+(1.0/2.0)*fShwfit.ph_hist->GetBinWidth(q));
     else if (E_split_point==0) E_split_point = (fShwfit.ph_hist->GetBinLowEdge(q)+fShwfit.ph_hist->GetBinWidth(q));
     //cout<<"n_split_point = "<<n_split_point<<" n_split_val = "<<n_split_val<<" 1/2 Integral = "<<n_integral<<endl;
     //cout<<"E_split_point = "<<E_split_point<<" E_split_val = "<<E_split_val<<" 1/2 Integral = "<<E_integral<<endl;
   }
   if(E_half_den!=0) fShwfit.E_ratio_half = E_half_num/E_half_den;
   else fShwfit.E_ratio_half = 0;
   if(n_half_den!=0) fShwfit.n_ratio_half = n_half_num/n_half_den;
   else fShwfit.n_ratio_half = 0;
   if(E_2_den!=0) fShwfit.E_ratio_2 = E_2_num/E_2_den;
   else fShwfit.E_ratio_2 = 0;
   if(n_2_den!=0) fShwfit.n_ratio_2 = n_2_num/n_2_den;
   else fShwfit.n_ratio_2 = 0;
   fShwfit.pos_E_split = E_split_point;
   fShwfit.pos_n_split = n_split_point;
   

   delete [] x;
   delete [] x3;
   delete [] w;
   delete [] w3;
 
   return;
}


float ShwfitAna::BuildUVVar(float u, float v)
{
   float uv = ANtpDefVal::kFloat;
   if(!ANtpDefVal::IsDefault(u)&& !ANtpDefVal::IsDefault(v)&&
          u<1.e15 && v<1.e15){
      uv=sqrt(u*u+v*v)/sqrt(2.0);
   }
                                                                                
   return uv;
}                                                         

//--------------------------------------------
//DDC Transverse fit Code
void ShwfitAna::FitTShower(Float_t /*pulseheight*/)
{
  fShwfit.tenestu->Fit(fShwfit.ufit,"RLQ0+");
  fShwfit.trans_u_mean = fShwfit.ufit->GetParameter(1);
  fShwfit.trans_u_sigma = fShwfit.ufit->GetParameter(2);
  fShwfit.trans_u_chisq = fShwfit.ufit->GetChisquare();
  fShwfit.trans_u_ndf = fShwfit.ufit->GetNDF();
                                                                                                   
  fShwfit.tenestv->Fit(fShwfit.vfit,"RLQ0+");
  fShwfit.trans_v_mean = fShwfit.vfit->GetParameter(1);
  fShwfit.trans_v_sigma = fShwfit.vfit->GetParameter(2);
  fShwfit.trans_v_chisq = fShwfit.vfit->GetChisquare();
  fShwfit.trans_v_ndf = fShwfit.vfit->GetNDF();
                                                                                                   
}
//----------------------------------------------------
                                                                                               
void ShwfitAna::TransVar(TH1F *histo, PlaneView::EPlaneView pv)
{
   MSG("ShwfitAna",Msg::kDebug)<<"In ShwfitAna::TransVar"<<endl;
   MSG("ShwfitAna",Msg::kDebug)<<"plane view is "<<pv<<endl;
   MSG("ShwfitAna",Msg::kDebug)<<"Hist name is "<<histo->GetName()<<endl;
   MSG("ShwfitAna",Msg::kDebug)<<"Hist has "<<histo->GetEntries()<<" entries"<<endl;
   
   Int_t THISTBINS = (Int_t)(histo->GetNbinsX());

  Int_t binmax=histo->GetMaximumBin();
  Int_t binvert=Int_t(((THISTBINS-1)/2)+1);
  Float_t peak_bin=histo->Integral(binmax,binmax);
  Float_t vert_bin=histo->Integral(binvert,binvert);
  Float_t tot=histo->Integral(1,THISTBINS);

  molrad_peak = molrad_vert = 0;
  mean =  rms =  skew = kurt = 0;

  asym_peak= ANtpDefVal::kFloat;
  asym_vert= ANtpDefVal::kFloat;

  if(tot-peak_bin){
    asym_peak=(TMath::Abs(histo->Integral(binmax+1,THISTBINS)-histo->Integral(1,binmax-1)))
    /(tot-peak_bin);
  }

  if(tot-vert_bin){
    asym_vert=(TMath::Abs(histo->Integral(binvert+1,THISTBINS)-histo->Integral(1,binvert-1)))
    /(tot-vert_bin);
  }

  Float_t ratio;

  if(tot){
    for(Int_t i=0; i<=binvert-1;i++){
      ratio=histo->Integral(binmax-i>0?binmax-i:1,
                            binmax+i<THISTBINS?binmax+i:THISTBINS)/tot;
      if(ratio>0.90){molrad_peak=i+1; break;}
    }
    for(Int_t i=0; i<=binvert-1;i++){
      ratio=histo->Integral(binvert-i,binvert+i)/tot;
      if(ratio>0.90){molrad_vert=i+1; break;}
    }

  }

  mean=histo->GetMean();
  rms=histo->GetRMS();
  Int_t n_count=0;

  for(Int_t i=1;i<=THISTBINS;i++){

    if(histo->GetBinContent(i)){
      skew=skew+(TMath::Power((histo->GetBinCenter(i)-mean),3)
		     *histo->GetBinContent(i));
      kurt=kurt+(TMath::Power((histo->GetBinCenter(i)-mean),4)
		     *histo->GetBinContent(i));
      n_count++;
    }

  }
  if(rms>0 && n_count>1){
    skew=skew/((Float_t)(n_count-1)*TMath::Power((rms),3));
    kurt=(kurt/((Float_t)(n_count-1)*TMath::Power((rms),4)))-3;
  }
  else{skew= ANtpDefVal::kFloat; kurt= ANtpDefVal::kFloat;}

}

void ShwfitAna::Reset(int snarl, int event)
{

//putting histogram creators here so we can name them according to
//event and snarl number, and we won't get errors about replacing 
//existing histograms when we have multiple events in a snarl
   const Int_t LHISTBINS=30;
   const Int_t THISTBINS=41;
   fShwfit.Reset();

   char ln[100];
   char tun[100];
   char tvn[100];
   sprintf(ln,"lenepl_%d_%d",snarl,event);
   sprintf(tun,"tenestu_%d_%d",snarl,event);
   sprintf(tvn,"tenestv_%d_%d",snarl,event);
   fShwfit.lenepl = new TH1F(ln,"longitudinal energy by plane",
			     LHISTBINS,0.0,LHISTBINS);
   
   fShwfit.tenestu = new TH1F(tun,"trasverse energy by strip (U view)",
			      THISTBINS,-THISTBINS/2.,THISTBINS/2.);
   fShwfit.tenestv = new TH1F(tvn,"trasverse energy by strip (V view)",
			      THISTBINS,-THISTBINS/2.,THISTBINS/2.);

   char tun_9s_2pe[100];
   char tvn_9s_2pe[100];
   char tun_9s_2pe_dw[100];
   char tvn_9s_2pe_dw[100];

   sprintf(tun_9s_2pe,"tenestu_9s_2pe_%d_%d",snarl,event);
   sprintf(tvn_9s_2pe,"tenestv_9s_2pe_%d_%d",snarl,event);
   sprintf(tun_9s_2pe_dw,"tenestu_9s_2pe_dw_%d_%d",snarl,event);
   sprintf(tvn_9s_2pe_dw,"tenestv_9s_2pe_dw_%d_%d",snarl,event);

   fShwfit.tenestu_9s_2pe = new TH1F(tun_9s_2pe,"trasverse energy by strip (U view)",
                                  THISTBINS,-THISTBINS/2.,THISTBINS/2.);
   fShwfit.tenestv_9s_2pe = new TH1F(tvn_9s_2pe,"trasverse energy by strip (V view)",
                                  THISTBINS,-THISTBINS/2.,THISTBINS/2.);
   fShwfit.tenestu_9s_2pe_dw = new TH1F(tun_9s_2pe_dw,"trasverse energy by strip (U view)",
                                  THISTBINS,-THISTBINS/2.,THISTBINS/2.);
   fShwfit.tenestv_9s_2pe_dw = new TH1F(tvn_9s_2pe_dw,"trasverse energy by strip (V view)",
                                  THISTBINS,-THISTBINS/2.,THISTBINS/2.);


    
   int hmin=0;
   int hmax=30;
   int npare=3;
   int nparh=6;
   char efn[100];
   char hfn[100];
   sprintf(efn,"efit_%d_%d",snarl,event);
   sprintf(hfn,"hfit_%d_%d",snarl,event);
   fShwfit.efit = new TF1(efn,shwfunc,hmin+0.001,hmax,npare);
   fShwfit.hfit = new TF1(hfn,hadfunc,hmin+0.001,hmax,nparh);
   fShwfit.efit->SetParNames("a","b","e0");
   fShwfit.efit->SetParLimits(0,hmin+0.001,hmax);
   fShwfit.efit->SetParLimits(1,0.001,20000);
   fShwfit.efit->SetParLimits(2,0+0.001,1000000);
   
   fShwfit.hfit->SetParNames("a","b","e0");
   fShwfit.hfit->SetParLimits(0,hmin+0.001,hmax);
   fShwfit.hfit->SetParLimits(1,0.001,20000);
   fShwfit.hfit->SetParLimits(2,0.001,1000000);

   //--------------------------------------------
   //DDC Mawell fit Code
   char phh[100];

   sprintf(phh,"ph_hist_%d_%d",snarl,event);
   fShwfit.ph_hist = new TH1F(phh,"energy per srip", 100,0.5,10.5);

   int npar_maxwell=2;
   char mfn[100];
   sprintf(mfn,"efit_maxwell_%d_%d",snarl,event);
   fShwfit.efit_maxwell = new TF1(mfn,maxwell,hmin+0.001,hmax,npar_maxwell);
   fShwfit.efit_maxwell->SetParNames("Beta","Energy");
   fShwfit.efit_maxwell->SetParLimits(0,hmin+0.001,hmax);
                                                                                                   
   int npar_maxwell3=2;
   char mfn3[100];
   sprintf(mfn3,"efit_maxwell3_%d_%d",snarl,event);
   fShwfit.efit_maxwell3 = new TF1(mfn3,maxwell3,hmin+0.001,hmax,npar_maxwell3);
   fShwfit.efit_maxwell3->SetParNames("Beta","Energy");
   fShwfit.efit_maxwell3->SetParLimits(0,hmin+0.001,hmax);
                                                                                                   
   char ufn[100];
   sprintf(ufn,"ufit_%d_%d",snarl,event);
   fShwfit.ufit = new TF1(ufn,"gaus",-20,20);
                                                                                                   
   char vfn[100];
   sprintf(vfn,"vfit_%d_%d",snarl,event);
   fShwfit.vfit = new TF1(vfn,"gaus",-20,20);
   
}
Bool_t ShwfitAna::PassCuts(int PhNStrips, int PhNPlanes )
{
  bool passstrips=false;
  bool passplanes=false;
  if(PhNStrips>0&&fShwfit.hiPhStripCount>PhNStrips||PhNStrips<=0) passstrips=true;
  if(PhNPlanes>0&&fShwfit.hiPhPlaneCount>PhNPlanes||PhNPlanes<=0) passplanes=true;
  cout << "StripCount " << fShwfit.hiPhStripCount << " PlaneCount " << fShwfit.hiPhPlaneCount << endl; 
  if(passstrips&&passplanes) return true;
  else return false;
}

Double_t ShwfitAna::GetMaximumX(TF1* efit, Double_t xmin, Double_t xmax)
{                           
   Double_t fXmin = 0.001;
   Double_t fXmax = 30;
   Double_t fNpx = 30;                                                    

   Double_t x,y;
   if (xmin >= xmax) {xmin = fXmin; xmax = fXmax;}
   Double_t dx = (xmax-xmin)/fNpx;
   Double_t xxmax = xmin;
   Double_t yymax = efit->Eval(xmin+dx);
   for (Int_t i=0;i<fNpx;i++) {
      x = xmin + (i+0.5)*dx;
      y = efit->Eval(x);
      if (y > yymax) {xxmax = x; yymax = y;}
   }
   if (dx < 1.e-9*(fXmax-fXmin)) return TMath::Min(xmax,xxmax);
   else return GetMaximumX(efit, TMath::Max(xmin,xxmax-dx), TMath::Min(xmax,xxmax+dx));
}


void ShwfitAna::doSlopes( NtpSREvent *event , RecRecordImp<RecCandHeader> *srobj)
 {
   
   

   Float_t Uslope=0;
   Float_t Uuncer=0;
   Float_t UOffset=0;
   Float_t Ugoodfit=0;

   Float_t Vslope=0;
   Float_t Vuncer=0;
   Float_t VOffset=0;
   Float_t Vgoodfit=0;

   Float_t USlopeMom=0;
   Float_t VSlopeMom=0;

   Float_t UBeamLike=0;
   Float_t VBeamLike=0;

   Float_t UBeamLikeOffset=0;
   Float_t VBeamLikeOffset=0;

   Float_t VWd=0;

   Float_t UWd=0;


   Float_t num[1000];
   Float_t den[1000];
   Float_t zpos[1000];
   Float_t wcenter[1000];
   
   Int_t stripcount[1000];
   Float_t stripe[1000];

   Int_t umin=10000;
   Int_t umax=0;
   Int_t vmin=10000;
   Int_t vmax=0;
   Int_t curplane;
   Float_t curstrip;

   Float_t ULongE=0;
   Float_t VLongE=0;

   Float_t tote=0.0;

   Float_t weight[1000];
   for(int i=0;i<1000;i++){
     num[i]=0;
     den[i]=0;
     wcenter[i]=0;
     weight[i]=0;
     zpos[i]=0.0;
     stripcount[i]=0;
     stripe[i]=0;
   }
   



   //the slope of the beam 3deg inc in Y for far and 3deg dec in Y for near
   Float_t beamslope=0;
   if(srobj->GetHeader().GetVldContext().GetDetector()== Detector::kFar)
     beamslope=0.037058;
   else if(srobj->GetHeader().GetVldContext().GetDetector()== Detector::kNear)
     beamslope=-0.037058;


     for(int i=0;i<event->nstrip;i++){
      Int_t index = SntpHelpers::GetStripIndex(i,event);
      NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
      if(!strip) continue;
      if(!evtstp0mip){
        MSG("ShwfitAna",Msg::kError)<<"No mip strip information"<<endl;
        continue;
      }

        Float_t stripPh = evtstp0mip[index] + evtstp1mip[index];
     
        curplane=strip->plane;

	stripcount[curplane]++;
	stripe[curplane]+=stripPh;

         curstrip=strip->tpos;
	 zpos[curplane]=strip->z;

        num[curplane]+=curstrip*stripPh;
        den[curplane]+=stripPh;

	tote+=stripPh;
     
        if(strip->planeview==PlaneView::kU)
       {
	 if(umin>curplane)umin=curplane;
	 if(umax<curplane)umax=curplane;
	 USlopeMom+=stripPh;
	 if((strip->z-event->vtx.z)>0.0)
	 ULongE+=stripPh*TMath::Cos(TMath::ATan((strip->tpos-event->vtx.u)/(strip->z-event->vtx.z)));
	 // else
	 //  ULongE+=stripPh;
       }
     else if(strip->planeview==PlaneView::kV)
       {
	 if(vmin>curplane)vmin=curplane;
	 if(vmax<curplane)vmax=curplane;
	 VSlopeMom+=stripPh;
	 if((strip->z-event->vtx.z)>0.0)//only want hits in front of vertex
	 VLongE+=stripPh*TMath::Cos(TMath::ATan((strip->tpos-event->vtx.v)/(strip->z-event->vtx.z)));
	 // else
	 //  VLongE+=stripPh;
       }
     
     }



   for(int i=umin;i<=umax;i+=2)
     {
       if(den[i]!=0){
	 wcenter[i]=num[i]/den[i];
	 UBeamLikeOffset+=den[i]*(wcenter[i]-beamslope*zpos[i]);  //den is the energy in each plane
       }else{
	 num[i]=0;
	 wcenter[i]=0;
       }

     }
  
   if(USlopeMom>0){
   UBeamLikeOffset=UBeamLikeOffset/USlopeMom; //USlopeMom currently contains total energy in UPlane
   }

   for(int i=vmin;i<=vmax;i+=2)
     {
       if(den[i]!=0){
	 wcenter[i]=num[i]/den[i];
	 VBeamLikeOffset+=den[i]*(wcenter[i]-beamslope*zpos[i]);  //den is the energy in each plane
       }else{
	 num[i]=0;
	 wcenter[i]=0;
       }
     }
   
   if(VSlopeMom>0){
   VBeamLikeOffset=VBeamLikeOffset/VSlopeMom; //VSlopeMom currently contains total energy in UPlane
   }



Int_t stripmin=umin; if(vmin<umin)stripmin=vmin;
Int_t stripmax=umax; if(vmax>umax)stripmax=vmax;

Float_t complex=0;
Float_t wcomplex=0;
 Float_t ncomplex=0;

for (int i=stripmin;i<=stripmax;i++)
{
        complex+=stripcount[i]*stripcount[i+1];
	wcomplex+=stripcount[i]*stripcount[i+1]*stripe[i]*stripe[i+1];
	ncomplex+=stripe[i]*stripe[i+1];
}   

 if(ncomplex>0) wcomplex=wcomplex/ncomplex;
 else wcomplex=0;



   Float_t VSw=0;
   Float_t VSwxy=0;
   Float_t VSwx=0;
   Float_t VSwy=0;
   Float_t VSwxx=0;
   
   Float_t USw=0;
   Float_t USwxy=0;
   Float_t USwx=0;
   Float_t USwy=0;
   Float_t USwxx=0;




   Float_t pos;




   for(int i=umin;i<=umax;i+=2){

     if (i>999)continue;

     UBeamLike+=den[i]*den[i]*(wcenter[i]-UBeamLikeOffset-beamslope*zpos[i])*(wcenter[i]-UBeamLikeOffset-beamslope*zpos[i]);

             
     if(den[i]>0){
       
       pos=num[i]/(den[i]);
       weight[i]=den[i]/tote;
       
       USw=USw+weight[i];
       USwxy=USwxy+weight[i]*zpos[i]*pos;
       USwx=USwx+weight[i]*zpos[i];
       USwy=USwy+weight[i]*pos;
       USwxx=USwxx+weight[i]*zpos[i]*zpos[i];
       

     }

	 
   }

   if(USlopeMom>0){
     UBeamLike=TMath::Sqrt(UBeamLike)/USlopeMom;    //USlopeMom still contains total energy in U plane
   } 

 
   for(int i=vmin;i<=vmax;i+=2){
     
     
     if (i>999)continue;
     
      VBeamLike+=den[i]*den[i]*(wcenter[i]-VBeamLikeOffset-beamslope*zpos[i])*(wcenter[i]-VBeamLikeOffset-beamslope*zpos[i]);

     
     if(den[i]>0){
       
       pos=num[i]/(den[i]);
       weight[i]=den[i]/tote;
       
       VSw=VSw+weight[i];
       VSwxy=VSwxy+weight[i]*zpos[i]*pos;
       VSwx=VSwx+weight[i]*zpos[i];
       VSwy=VSwy+weight[i]*pos;
       VSwxx=VSwxx+weight[i]*zpos[i]*zpos[i];
       

       
     }

       
       
   }

   if(VSlopeMom>0){
    VBeamLike=TMath::Sqrt(VBeamLike)/VSlopeMom;    //VSlopeMom still contains total energy in V plane
   }

   Float_t Udelta=(USw*USwxx-USwx*USwx);

   if (Udelta>0){
     Uslope=(USw*USwxy-USwx*USwy)/Udelta;
     Uuncer=TMath::Sqrt(USw/Udelta);
     
     
     
     UOffset=(USwxx*USwy-USwx*USwxy)/Udelta;
    
     
     Float_t Uvtx=event->vtx.u;

    
     
      Ugoodfit=0;
     Float_t SumWeight=0;
     
     for(int i=umin;i<=umax;i+=2){
       if (i>999)continue;
       if(den[i]>0){
	 Ugoodfit=Ugoodfit+(UOffset+Uslope*zpos[i])*(UOffset+Uslope*zpos[i])/(weight[i]*weight[i]);
	 SumWeight=SumWeight+1/weight[i];

	  UWd=UWd+(UOffset+Uslope*zpos[i]-Uvtx)*(UOffset+Uslope*zpos[i]-Uvtx)/(weight[i]*weight[i]);
	 
       }
     }
     
     Ugoodfit=TMath::Sqrt(Ugoodfit)/SumWeight;
     
     UWd=TMath::Sqrt(UWd)/SumWeight;
     
     
     

   }
   
   
   
   
   Float_t Vdelta=(VSw*VSwxx-VSwx*VSwx);

   if (Vdelta>0){
     Vslope=(VSw*VSwxy-VSwx*VSwy)/Vdelta;
     Vuncer=TMath::Sqrt(VSw/Vdelta);
     
     
     
     VOffset=(VSwxx*VSwy-VSwx*VSwxy)/Vdelta;

     
     Float_t Vvtx=event->vtx.v;
     
     Vgoodfit=0;
     Float_t SumWeight=0;
     for(int i=vmin;i<=vmax;i+=2){
       if (i>999)continue;
       if(den[i]>0){
	 Vgoodfit=Vgoodfit+(VOffset+Vslope*zpos[i])*(VOffset+Vslope*zpos[i])/(weight[i]*weight[i]);
	 SumWeight=SumWeight+1/weight[i];

	 VWd=VWd+(VOffset+Vslope*zpos[i]-Vvtx)*(VOffset+Vslope*zpos[i]-Vvtx)/(weight[i]*weight[i]);

       }
     }
     
     Vgoodfit=TMath::Sqrt(Vgoodfit)/SumWeight;
     
          VWd=TMath::Sqrt(VWd)/SumWeight;
   }
   

   if(USlopeMom>0){
     USlopeMom=1/USlopeMom;
     USlopeMom=USlopeMom*Uslope;
   }

   if(VSlopeMom>0){
     VSlopeMom=1/VSlopeMom;
     VSlopeMom=VSlopeMom*Vslope;
   }
   
   
   fShwfit.UBeamLike=UBeamLike;
   fShwfit.VBeamLike=VBeamLike;

   fShwfit.slopefix=beamslope;
   fShwfit.USlope=Uslope;
   fShwfit.VSlope=Vslope;
   fShwfit.UOffset=UOffset;
   fShwfit.VOffset=VOffset;
   fShwfit.UFitquality=Ugoodfit;
   fShwfit.VFitquality=Vgoodfit;
   fShwfit.UWDiff=UWd;
   fShwfit.VWDiff=VWd;
   fShwfit.USlopeMom=USlopeMom;
   fShwfit.VSlopeMom=VSlopeMom;
   fShwfit.UVSlope=Uslope*Uslope+Vslope*Vslope;

   fShwfit.ULongE=ULongE;
   fShwfit.VLongE=VLongE;
   fShwfit.LongE = ULongE + VLongE;
   

   fShwfit.complexity=complex;

   fShwfit.wcomplexity=wcomplex;

   return;
 }




