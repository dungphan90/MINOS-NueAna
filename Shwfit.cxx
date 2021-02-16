/**
 *
 * $Id: Shwfit.cxx,v 1.21 2008/07/28 17:56:42 boehm Exp $
 *
 * \class Shwfit
 *
 * \package NueAna
 *
 * \brief Hold variables related to the ShwfitAna package
 **/
#include <iostream>
#include "TH1F.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "NueAna/Shwfit.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "MessageService/MsgService.h"

CVSID("$Id: Shwfit.cxx,v 1.21 2008/07/28 17:56:42 boehm Exp $");

ClassImp(Shwfit)

Shwfit::Shwfit():
   par_a(ANtpDefaultValue::kFloat),
   par_b(ANtpDefaultValue::kFloat),
   par_e0(ANtpDefaultValue::kFloat),
   chisq(ANtpDefaultValue::kFloat),
   shwmax(ANtpDefaultValue::kFloat),
   shwmaxplane(ANtpDefaultValue::kInt),
   conv(ANtpDefaultValue::kInt),
   chisq_ndf(ANtpDefaultValue::kFloat),
   e0_pe_ratio(ANtpDefaultValue::kFloat),
   caldet_comp(ANtpDefaultValue::kFloat),
   max_pe_plane(ANtpDefaultValue::kFloat),
   shwmaxplane_diff(ANtpDefaultValue::kFloat),
   hiPhStripCountM4(ANtpDefaultValue::kInt),
   hiPhPlaneCountM4(ANtpDefaultValue::kInt),
   hiPhStripCountM2(ANtpDefaultValue::kInt),
   hiPhPlaneCountM2(ANtpDefaultValue::kInt),
   hiPhStripCount(ANtpDefaultValue::kInt),
   hiPhPlaneCount(ANtpDefaultValue::kInt),
   hiPhStripCountP2(ANtpDefaultValue::kInt),
   hiPhPlaneCountP2(ANtpDefaultValue::kInt),
   hiPhStripCountP4(ANtpDefaultValue::kInt),
   hiPhPlaneCountP4(ANtpDefaultValue::kInt),
   // Dan's code
   Beta_Maxwell(ANtpDefaultValue::kFloat),
   Energy_Maxwell(ANtpDefaultValue::kFloat),
   chisq_Maxwell(ANtpDefaultValue::kFloat),
   ndf_Maxwell(ANtpDefaultValue::kFloat),
   Beta_Maxwell3(ANtpDefaultValue::kFloat),
   Energy_Maxwell3(ANtpDefaultValue::kFloat),
   chisq_Maxwell3(ANtpDefaultValue::kFloat),
   ndf_Maxwell3(ANtpDefaultValue::kFloat),
   trans_u_mean(ANtpDefaultValue::kFloat),
   trans_u_sigma(ANtpDefaultValue::kFloat),
   trans_u_chisq(ANtpDefaultValue::kFloat),
   trans_u_ndf(ANtpDefaultValue::kFloat),
   trans_v_mean(ANtpDefaultValue::kFloat),
   trans_v_sigma(ANtpDefaultValue::kFloat),
   trans_v_chisq(ANtpDefaultValue::kFloat),
   trans_v_ndf(ANtpDefaultValue::kFloat),
   E_ratio_half(ANtpDefaultValue::kFloat),
   n_ratio_half(ANtpDefaultValue::kFloat),
   E_ratio_2(ANtpDefaultValue::kFloat),
   n_ratio_2(ANtpDefaultValue::kFloat),
   pos_E_split(ANtpDefaultValue::kFloat),
   pos_n_split(ANtpDefaultValue::kFloat),
   
  // new contPlaneCount - Minerba
   contPlaneCount(ANtpDefaultValue::kInt),
   contPlaneCount015(ANtpDefaultValue::kInt),
   contPlaneCount030(ANtpDefaultValue::kInt),
   contPlaneCount050(ANtpDefaultValue::kInt),
   contPlaneCount075(ANtpDefaultValue::kInt),
   contPlaneCount100(ANtpDefaultValue::kInt),
   contPlaneCount200(ANtpDefaultValue::kInt),

   u_asym_peak(ANtpDefaultValue::kFloat),
   u_asym_vert(ANtpDefaultValue::kFloat),
   u_molrad_peak(ANtpDefaultValue::kFloat),
   u_molrad_vert(ANtpDefaultValue::kFloat),
   u_mean(ANtpDefaultValue::kFloat),
   u_rms(ANtpDefaultValue::kFloat),
   u_skew(ANtpDefaultValue::kFloat),
   u_kurt(ANtpDefaultValue::kFloat),
   v_asym_peak(ANtpDefaultValue::kFloat),
   v_asym_vert(ANtpDefaultValue::kFloat),
   v_molrad_peak(ANtpDefaultValue::kFloat),
   v_molrad_vert(ANtpDefaultValue::kFloat),
   v_mean(ANtpDefaultValue::kFloat),
   v_rms(ANtpDefaultValue::kFloat),
   v_skew(ANtpDefaultValue::kFloat),
   v_kurt(ANtpDefaultValue::kFloat),
   uv_asym_peak(ANtpDefaultValue::kFloat),
   uv_asym_vert(ANtpDefaultValue::kFloat),
   uv_molrad_peak(ANtpDefaultValue::kFloat),
   uv_molrad_vert(ANtpDefaultValue::kFloat),
   uv_mean(ANtpDefaultValue::kFloat),
   uv_rms(ANtpDefaultValue::kFloat),
   uv_skew(ANtpDefaultValue::kFloat),
   uv_kurt(ANtpDefaultValue::kFloat),
   uv_ratio(ANtpDefaultValue::kFloat),
   vtxEnergy(ANtpDefaultValue::kFloat),
   energyPlane0(ANtpDefaultValue::kFloat),
   energyPlane1(ANtpDefaultValue::kFloat),
   energyPlane2(ANtpDefaultValue::kFloat),
USlope(ANtpDefaultValue::kFloat),
VSlope(ANtpDefaultValue::kFloat),
UOffset(ANtpDefaultValue::kFloat),
VOffset(ANtpDefaultValue::kFloat),
UFitquality(ANtpDefaultValue::kFloat),
VFitquality(ANtpDefaultValue::kFloat),
UWDiff(ANtpDefaultValue::kFloat),
VWDiff(ANtpDefaultValue::kFloat),
USlopeMom(ANtpDefaultValue::kFloat),
VSlopeMom(ANtpDefaultValue::kFloat),
slopefix(ANtpDefaultValue::kFloat),
UBeamLike(ANtpDefaultValue::kFloat),
VBeamLike(ANtpDefaultValue::kFloat),
UVSlope(ANtpDefaultValue::kFloat),
ULongE(ANtpDefaultValue::kFloat),
VLongE(ANtpDefaultValue::kFloat),
   LongE(ANtpDefaultValue::kFloat),
complexity(ANtpDefaultValue::kFloat),
   wcomplexity(ANtpDefaultValue::kFloat),
  u_asym_peak_9s_2pe(ANtpDefaultValue::kFloat),
   u_asym_vert_9s_2pe(ANtpDefaultValue::kFloat),
   u_molrad_peak_9s_2pe(ANtpDefaultValue::kFloat),
   u_molrad_vert_9s_2pe(ANtpDefaultValue::kFloat),
   u_mean_9s_2pe(ANtpDefaultValue::kFloat),
   u_rms_9s_2pe(ANtpDefaultValue::kFloat),
   u_skew_9s_2pe(ANtpDefaultValue::kFloat),
   u_kurt_9s_2pe(ANtpDefaultValue::kFloat),
   v_asym_peak_9s_2pe(ANtpDefaultValue::kFloat),
   v_asym_vert_9s_2pe(ANtpDefaultValue::kFloat),
   v_molrad_peak_9s_2pe(ANtpDefaultValue::kFloat),
   v_molrad_vert_9s_2pe(ANtpDefaultValue::kFloat),
   v_mean_9s_2pe(ANtpDefaultValue::kFloat),
   v_rms_9s_2pe(ANtpDefaultValue::kFloat),
   v_skew_9s_2pe(ANtpDefaultValue::kFloat),
   v_kurt_9s_2pe(ANtpDefaultValue::kFloat),
   uv_asym_peak_9s_2pe(ANtpDefaultValue::kFloat),
   uv_asym_vert_9s_2pe(ANtpDefaultValue::kFloat),
   uv_molrad_peak_9s_2pe(ANtpDefaultValue::kFloat),
   uv_molrad_vert_9s_2pe(ANtpDefaultValue::kFloat),
   uv_mean_9s_2pe(ANtpDefaultValue::kFloat),
   uv_rms_9s_2pe(ANtpDefaultValue::kFloat),
   uv_skew_9s_2pe(ANtpDefaultValue::kFloat),
   uv_kurt_9s_2pe(ANtpDefaultValue::kFloat),
   uv_ratio_9s_2pe(ANtpDefaultValue::kFloat),
                                                                                
   u_asym_peak_9s_2pe_dw(ANtpDefaultValue::kFloat),
   u_asym_vert_9s_2pe_dw(ANtpDefaultValue::kFloat),
   u_molrad_peak_9s_2pe_dw(ANtpDefaultValue::kFloat),
   u_molrad_vert_9s_2pe_dw(ANtpDefaultValue::kFloat),
   u_mean_9s_2pe_dw(ANtpDefaultValue::kFloat),
   u_rms_9s_2pe_dw(ANtpDefaultValue::kFloat),
   u_skew_9s_2pe_dw(ANtpDefaultValue::kFloat),
   u_kurt_9s_2pe_dw(ANtpDefaultValue::kFloat),
   v_asym_peak_9s_2pe_dw(ANtpDefaultValue::kFloat),
   v_asym_vert_9s_2pe_dw(ANtpDefaultValue::kFloat),
   v_molrad_peak_9s_2pe_dw(ANtpDefaultValue::kFloat),
   v_molrad_vert_9s_2pe_dw(ANtpDefaultValue::kFloat),
   v_mean_9s_2pe_dw(ANtpDefaultValue::kFloat),
   v_rms_9s_2pe_dw(ANtpDefaultValue::kFloat),
   v_skew_9s_2pe_dw(ANtpDefaultValue::kFloat),
   v_kurt_9s_2pe_dw(ANtpDefaultValue::kFloat),
   uv_asym_peak_9s_2pe_dw(ANtpDefaultValue::kFloat),
   uv_asym_vert_9s_2pe_dw(ANtpDefaultValue::kFloat),
   uv_molrad_peak_9s_2pe_dw(ANtpDefaultValue::kFloat),
   uv_molrad_vert_9s_2pe_dw(ANtpDefaultValue::kFloat),
   uv_mean_9s_2pe_dw(ANtpDefaultValue::kFloat),
   uv_rms_9s_2pe_dw(ANtpDefaultValue::kFloat),
   uv_skew_9s_2pe_dw(ANtpDefaultValue::kFloat),
   uv_kurt_9s_2pe_dw(ANtpDefaultValue::kFloat),
   uv_ratio_9s_2pe_dw(ANtpDefaultValue::kFloat),

   lenepl(0),
   ph_hist(0),
   tenestu(0),
   tenestv(0),
   tenestu_9s_2pe(0),
   tenestv_9s_2pe(0),
   tenestu_9s_2pe_dw(0),
   tenestv_9s_2pe_dw(0),
   efit(0),
   hfit(0),
   efit_maxwell(0),
   efit_maxwell3(0),
   ufit(0),
   vfit(0),
   info1(0)
{}

Shwfit::Shwfit(const Shwfit &s):
  TObject(),
  par_a(s.par_a),
  par_b(s.par_b),
  par_e0(s.par_e0),
  chisq(s.chisq),
  shwmax(s.shwmax),
  shwmaxplane(s.shwmaxplane),
  conv(s.conv),
  chisq_ndf(s.chisq_ndf),
  e0_pe_ratio(s.e0_pe_ratio),
  caldet_comp(s.caldet_comp),
  max_pe_plane(s.max_pe_plane),
  shwmaxplane_diff(s.shwmaxplane_diff),
  hiPhStripCountM4(s.hiPhStripCountM4),
  hiPhPlaneCountM4(s.hiPhPlaneCountM4),
  hiPhStripCountM2(s.hiPhStripCountM2),
  hiPhPlaneCountM2(s.hiPhPlaneCountM2),
  hiPhStripCount(s.hiPhStripCount),
  hiPhPlaneCount(s.hiPhPlaneCount),
  hiPhStripCountP2(s.hiPhStripCountP2),
  hiPhPlaneCountP2(s.hiPhPlaneCountP2),
  hiPhStripCountP4(s.hiPhStripCountP4),
  hiPhPlaneCountP4(s.hiPhPlaneCountP4),
  // dan var
  Beta_Maxwell(s.Beta_Maxwell),
  Energy_Maxwell(s.Energy_Maxwell),
  chisq_Maxwell(s.chisq_Maxwell),
  ndf_Maxwell(s.ndf_Maxwell),
  Beta_Maxwell3(s.Beta_Maxwell3),
  Energy_Maxwell3(s.Energy_Maxwell3),
  chisq_Maxwell3(s.chisq_Maxwell3),
  ndf_Maxwell3(s.ndf_Maxwell3),
  trans_u_mean(s.trans_u_mean),
  trans_u_sigma(s.trans_u_sigma),
  trans_u_chisq(s.trans_u_chisq),
  trans_u_ndf(s.trans_u_ndf),
  trans_v_mean(s.trans_v_mean),
  trans_v_sigma(s.trans_v_sigma),
  trans_v_chisq(s.trans_v_chisq),
  trans_v_ndf(s.trans_v_ndf),
  E_ratio_half(s.E_ratio_half),
  n_ratio_half(s.n_ratio_half),
  E_ratio_2(s.E_ratio_2),
  n_ratio_2(s.n_ratio_2),
  pos_E_split(s.pos_E_split),
  pos_n_split(s.pos_n_split),

  // contPlaneCount 
  contPlaneCount(s.contPlaneCount),
  contPlaneCount015(s.contPlaneCount015),
  contPlaneCount030(s.contPlaneCount030),
  contPlaneCount050(s.contPlaneCount050),
  contPlaneCount075(s.contPlaneCount075),
  contPlaneCount100(s.contPlaneCount100),
  contPlaneCount200(s.contPlaneCount200),

  u_asym_peak(s.u_asym_peak),
  u_asym_vert(s.u_asym_vert),
  u_molrad_peak(s.u_molrad_peak),
  u_molrad_vert(s.u_molrad_vert),
  u_mean(s.u_mean),
  u_rms(s.u_rms),
  u_skew(s.u_skew),
  u_kurt(s.u_kurt),
  v_asym_peak(s.v_asym_peak),
  v_asym_vert(s.v_asym_vert),
  v_molrad_peak(s.v_molrad_peak),
  v_molrad_vert(s.v_molrad_vert),
  v_mean(s.v_mean),
  v_rms(s.v_rms),
  v_skew(s.v_skew),
  v_kurt(s.v_kurt),
  uv_asym_peak(s.uv_asym_peak),
  uv_asym_vert(s.uv_asym_vert),
  uv_molrad_peak(s.uv_molrad_peak),
  uv_molrad_vert(s.uv_molrad_vert),
  uv_mean(s.uv_mean),
  uv_rms(s.uv_rms),
  uv_skew(s.uv_skew),
  uv_kurt(s.uv_kurt),
  uv_ratio(s.uv_ratio),
  vtxEnergy(s.vtxEnergy),
  energyPlane0(s.energyPlane0),
  energyPlane1(s.energyPlane1),
  energyPlane2(s.energyPlane2),
  USlope(s.USlope),
  VSlope(s.VSlope),
  UOffset(s.UOffset),
  VOffset(s.VOffset),
  UFitquality(s.UFitquality),
  VFitquality(s.VFitquality),
  UWDiff(s.UWDiff),
  VWDiff(s.VWDiff),
  USlopeMom(s.USlopeMom),
  VSlopeMom(s.VSlopeMom),
  slopefix(s.slopefix),
  UBeamLike(s.UBeamLike),
  VBeamLike(s.VBeamLike),
  UVSlope(s.UVSlope),
  ULongE(s.ULongE),
  VLongE(s.VLongE),
  LongE(s.LongE),
  complexity(s.complexity),
  wcomplexity(s.wcomplexity),
   u_asym_peak_9s_2pe(s.u_asym_peak_9s_2pe),
   u_asym_vert_9s_2pe(s.u_asym_vert_9s_2pe),
   u_molrad_peak_9s_2pe(s.u_molrad_peak_9s_2pe),
   u_molrad_vert_9s_2pe(s.u_molrad_vert_9s_2pe),
   u_mean_9s_2pe(s.u_mean_9s_2pe),
   u_rms_9s_2pe(s.u_rms_9s_2pe),
   u_skew_9s_2pe(s.u_skew_9s_2pe),
   u_kurt_9s_2pe(s.u_kurt_9s_2pe),
   v_asym_peak_9s_2pe(s.v_asym_peak_9s_2pe),
   v_asym_vert_9s_2pe(s.v_asym_vert_9s_2pe),
   v_molrad_peak_9s_2pe(s.v_molrad_peak_9s_2pe),
   v_molrad_vert_9s_2pe(s.v_molrad_vert_9s_2pe),
   v_mean_9s_2pe(s.v_mean_9s_2pe),
   v_rms_9s_2pe(s.v_rms_9s_2pe),
   v_skew_9s_2pe(s.v_skew_9s_2pe),
   v_kurt_9s_2pe(s.v_kurt_9s_2pe),
   uv_asym_peak_9s_2pe(s.uv_asym_peak_9s_2pe),
   uv_asym_vert_9s_2pe(s.uv_asym_vert_9s_2pe),
   uv_molrad_peak_9s_2pe(s.uv_molrad_peak_9s_2pe),
   uv_molrad_vert_9s_2pe(s.uv_molrad_vert_9s_2pe),
   uv_mean_9s_2pe(s.uv_mean_9s_2pe),
   uv_rms_9s_2pe(s.uv_rms_9s_2pe),
   uv_skew_9s_2pe(s.uv_skew_9s_2pe),
   uv_kurt_9s_2pe(s.uv_kurt_9s_2pe),
   uv_ratio_9s_2pe(s.uv_ratio_9s_2pe),
                                                                                
   u_asym_peak_9s_2pe_dw(s.u_asym_peak_9s_2pe_dw),
   u_asym_vert_9s_2pe_dw(s.u_asym_vert_9s_2pe_dw),
   u_molrad_peak_9s_2pe_dw(s.u_molrad_peak_9s_2pe_dw),
   u_molrad_vert_9s_2pe_dw(s.u_molrad_vert_9s_2pe_dw),
   u_mean_9s_2pe_dw(s.u_mean_9s_2pe_dw),
   u_rms_9s_2pe_dw(s.u_rms_9s_2pe_dw),
   u_skew_9s_2pe_dw(s.u_skew_9s_2pe_dw),
   u_kurt_9s_2pe_dw(s.u_kurt_9s_2pe_dw),
   v_asym_peak_9s_2pe_dw(s.v_asym_peak_9s_2pe_dw),
   v_asym_vert_9s_2pe_dw(s.v_asym_vert_9s_2pe_dw),
   v_molrad_peak_9s_2pe_dw(s.v_molrad_peak_9s_2pe_dw),
   v_molrad_vert_9s_2pe_dw(s.v_molrad_vert_9s_2pe_dw),
   v_mean_9s_2pe_dw(s.v_mean_9s_2pe_dw),
   v_rms_9s_2pe_dw(s.v_rms_9s_2pe_dw),
   v_skew_9s_2pe_dw(s.v_skew_9s_2pe_dw),
   v_kurt_9s_2pe_dw(s.v_kurt_9s_2pe_dw),
   uv_asym_peak_9s_2pe_dw(s.uv_asym_peak_9s_2pe_dw),
   uv_asym_vert_9s_2pe_dw(s.uv_asym_vert_9s_2pe_dw),
   uv_molrad_peak_9s_2pe_dw(s.uv_molrad_peak_9s_2pe_dw),
   uv_molrad_vert_9s_2pe_dw(s.uv_molrad_vert_9s_2pe_dw),
   uv_mean_9s_2pe_dw(s.uv_mean_9s_2pe_dw),
   uv_rms_9s_2pe_dw(s.uv_rms_9s_2pe_dw),
   uv_skew_9s_2pe_dw(s.uv_skew_9s_2pe_dw),
   uv_kurt_9s_2pe_dw(s.uv_kurt_9s_2pe_dw),
   uv_ratio_9s_2pe_dw(s.uv_ratio_9s_2pe_dw)

{
  if(s.lenepl!=0){
    lenepl=new TH1F(*(s.lenepl));
  }
  else{
    lenepl=0;
  }
  if(s.tenestu!=0){
    tenestu=new TH1F(*(s.tenestu));
  }
  else{
    tenestu=0;
  }
  if(s.tenestv!=0){
    tenestv=new TH1F(*(s.tenestv));
  }
  else{
    tenestv=0;
  }
  if(s.efit!=0){
    efit=new TF1(*(s.efit));
  }
  else{
    efit=0;
  }
  if(s.hfit!=0){
    hfit=new TF1(*(s.hfit));
  }
  else{
    hfit=0;
  }
  if(s.info1!=0){
    info1=new TPaveText(*(s.info1));
  }
  else{
    info1=0;
  }

  if(s.efit_maxwell!=0) { efit_maxwell=new TF1(*(s.efit_maxwell));  }
  else{efit_maxwell=0;  }
  if(s.efit_maxwell3!=0){ efit_maxwell3=new TF1(*(s.efit_maxwell3));}
  else{   efit_maxwell3=0;  }
  if(s.ufit!=0){  ufit=new TF1(*(s.ufit)); }
  else{           ufit=0;  }
  if(s.vfit!=0){  vfit=new TF1(*(s.vfit)); }
  else{           vfit=0;  }

  if(s.tenestv_9s_2pe!=0){    tenestv_9s_2pe=new TH1F(*(s.tenestv_9s_2pe));}
  else{ tenestv_9s_2pe=0;  }
  if(s.tenestv_9s_2pe_dw!=0){
       tenestv_9s_2pe_dw=new TH1F(*(s.tenestv_9s_2pe_dw));}
  else{ tenestv_9s_2pe_dw=0;  }
   
  if(s.tenestu_9s_2pe!=0){    tenestu_9s_2pe=new TH1F(*(s.tenestu_9s_2pe));}
  else{ tenestu_9s_2pe=0;  }
  if(s.tenestu_9s_2pe_dw!=0){ 
       tenestu_9s_2pe_dw=new TH1F(*(s.tenestu_9s_2pe_dw));}
  else{ tenestu_9s_2pe_dw=0;  }

  if(s.ph_hist!=0){
    ph_hist = new TH1F(*(s.ph_hist));
  }
  else{    ph_hist=0;  }
}


Shwfit::~Shwfit()
{
   if(lenepl!=0){
      delete lenepl;
      lenepl=0;
   }
   if(tenestu!=0){
      delete tenestu;
      tenestu=0;
   }
   if(tenestv!=0){
      delete tenestv;
      tenestv=0;
   }

   if(efit!=0){
      delete efit;
      efit=0;
   }
   if(hfit!=0){
      delete hfit;
      hfit=0;
   }
   if(info1!=0){
      delete info1;
      info1=0;
   }

   if(efit_maxwell!=0){      delete efit_maxwell;      efit_maxwell=0;   }
   if(efit_maxwell3!=0){     delete efit_maxwell3;     efit_maxwell3=0;  }
   if(ufit!=0){		     delete ufit;	       ufit=0;		 }
   if(vfit!=0){		     delete vfit;	       vfit=0;		 }

   if(tenestv_9s_2pe!=0)   {   delete tenestv_9s_2pe;  tenestv_9s_2pe=0;}
   if(tenestv_9s_2pe_dw!=0){   delete tenestv_9s_2pe_dw;  tenestv_9s_2pe_dw=0;}
   if(tenestu_9s_2pe!=0)   {   delete tenestu_9s_2pe;  tenestu_9s_2pe=0;}
   if(tenestu_9s_2pe_dw!=0){   delete tenestu_9s_2pe_dw;  tenestu_9s_2pe_dw=0;}

   if(ph_hist !=0) { delete ph_hist; ph_hist = 0; }

}


void Shwfit::Draw(TPad *pad)
{

   pad->cd(1);
   Print("test");

//   TPaveText *info1;
   if(info1) info1->Delete();
   info1= new TPaveText(.05,.05,0.95,0.95);

   info1->Clear();
   info1->SetFillColor(kWhite);
   info1->SetBorderSize(1);
   info1->SetTextAlign(13);
   info1->SetTextFont(42);
   char temp[100];

   sprintf(temp," ");
   info1->AddText(temp);   
   
   sprintf(temp," Long Fit:");
   info1->AddText(temp);   
   sprintf(temp,"          a =  %2.2f, b = %2.2f, E_0 = %2.2f ",par_a,par_b,par_e0); 
   info1->AddText(temp);
   sprintf(temp,"          chisq/ndf =  %2.2f, E0/MEU ratio = %2.2f, Conv = %d, ContPlane = %d ",chisq_ndf,e0_pe_ratio,conv,contPlaneCount); 
   info1->AddText(temp);

   sprintf(temp,"U  Trans: U RMS = %2.2f, U kurt = %2.2f, U MolRad = %2.2f", u_rms, u_kurt, u_molrad_vert); 
   info1->AddText(temp);

   sprintf(temp,"V  Trans: V RMS = %2.2f, V kurt = %2.2f, V MolRad = %2.2f", v_rms, v_kurt, v_molrad_vert); 
   info1->AddText(temp);

   sprintf(temp,"UV Trans: UV RMS = %2.2f, UV kurt = %2.2f, UV MolRad = %2.2f", uv_rms, uv_kurt, uv_molrad_vert); 
   info1->AddText(temp);

   info1->Draw();

   pad->cd(3);
   lenepl->Draw();
   if(conv==1) efit->Draw("sames");
   pad->cd(2);

   tenestu->Draw();
   pad->cd(4);
   tenestv->Draw();
   pad->cd();
   pad->Modified();


}
void Shwfit::Draw(Option_t */*option*/)
{
/// useless function right now... 
}


void Shwfit::Print(Option_t */*option*/) const
{
    MSG("Shwfit::Print",Msg::kDebug) << "par_a = " << par_a   
             << " par_b = " << par_b
             << " par_e0 = " << par_e0
             << " chisq_ndf = " << chisq_ndf
             << " e0_pe_ratio = " << e0_pe_ratio
             << " conv = " << conv
             << endl;
}

void Shwfit::Clear(Option_t */*option*/)
{
  Reset();
  if(lenepl!=0) { delete lenepl;  lenepl=0;  }
  if(tenestu!=0){ delete tenestu; tenestu=0; }
  if(tenestv!=0){ delete tenestv; tenestv=0; }
  if(efit!=0){    delete efit;    efit=0;    }
  if(hfit!=0){    delete hfit;    hfit=0;    }
  if(info1!=0){   delete info1;   info1=0;   }
  if(tenestv_9s_2pe!=0)   {   delete tenestv_9s_2pe;  tenestv_9s_2pe=0;}
  if(tenestv_9s_2pe_dw!=0){   delete tenestv_9s_2pe_dw;  tenestv_9s_2pe_dw=0;}
  if(tenestu_9s_2pe!=0)   {   delete tenestu_9s_2pe;  tenestu_9s_2pe=0;}
  if(tenestu_9s_2pe_dw!=0){   delete tenestu_9s_2pe_dw;  tenestu_9s_2pe_dw=0;}

  if(efit_maxwell!=0){    delete efit_maxwell;    efit_maxwell=0;    }
  if(efit_maxwell3!=0){    delete efit_maxwell3;    efit_maxwell3=0;    }
  if(ufit!=0){    delete ufit;    ufit=0;    }
  if(vfit!=0){    delete vfit;    vfit=0;    }
   if(ph_hist !=0) { delete ph_hist; ph_hist = 0; }
  
}

void Shwfit::Reset()
{
   par_a= ANtpDefaultValue::kFloat;
   par_b= ANtpDefaultValue::kFloat;
   par_e0=ANtpDefaultValue::kFloat;
   chisq= ANtpDefaultValue::kFloat;
   shwmax=ANtpDefaultValue::kFloat;
   shwmaxplane= ANtpDefaultValue::kInt;
   conv= ANtpDefaultValue::kInt;
   chisq_ndf=ANtpDefaultValue::kFloat;
   e0_pe_ratio=ANtpDefaultValue::kFloat;
   caldet_comp=ANtpDefaultValue::kFloat;
   max_pe_plane=ANtpDefaultValue::kFloat;
   shwmaxplane_diff=ANtpDefaultValue::kFloat;

   Beta_Maxwell= ANtpDefaultValue::kFloat;
   Energy_Maxwell= ANtpDefaultValue::kFloat;
   chisq_Maxwell= ANtpDefaultValue::kFloat;
   ndf_Maxwell= ANtpDefaultValue::kFloat;
   Beta_Maxwell3= ANtpDefaultValue::kFloat;
   Energy_Maxwell3= ANtpDefaultValue::kFloat;
   chisq_Maxwell3= ANtpDefaultValue::kFloat;
   ndf_Maxwell3= ANtpDefaultValue::kFloat;
   trans_u_mean= ANtpDefaultValue::kFloat;
   trans_u_sigma= ANtpDefaultValue::kFloat;
   trans_u_chisq= ANtpDefaultValue::kFloat;
   trans_u_ndf= ANtpDefaultValue::kFloat;
   trans_v_mean= ANtpDefaultValue::kFloat;
   trans_v_sigma= ANtpDefaultValue::kFloat;
   trans_v_chisq= ANtpDefaultValue::kFloat;
   trans_v_ndf= ANtpDefaultValue::kFloat;
   E_ratio_half= ANtpDefaultValue::kFloat;
   n_ratio_half= ANtpDefaultValue::kFloat;
   E_ratio_2= ANtpDefaultValue::kFloat;
   n_ratio_2= ANtpDefaultValue::kFloat;
   pos_E_split= ANtpDefaultValue::kFloat;
   pos_n_split= ANtpDefaultValue::kFloat;

   hiPhStripCountM4=ANtpDefaultValue::kInt;
   hiPhPlaneCountM4=ANtpDefaultValue::kInt;
   hiPhStripCountM2=ANtpDefaultValue::kInt;
   hiPhPlaneCountM2=ANtpDefaultValue::kInt;
   hiPhStripCount=ANtpDefaultValue::kInt;
   hiPhPlaneCount=ANtpDefaultValue::kInt;
   hiPhStripCountP2=ANtpDefaultValue::kInt;
   hiPhPlaneCountP2=ANtpDefaultValue::kInt;
   hiPhStripCountP4=ANtpDefaultValue::kInt;
   hiPhPlaneCountP4=ANtpDefaultValue::kInt;
   contPlaneCount=ANtpDefaultValue::kInt;
   contPlaneCount015=ANtpDefaultValue::kInt;
   contPlaneCount030=ANtpDefaultValue::kInt;
   contPlaneCount050=ANtpDefaultValue::kInt;
   contPlaneCount075=ANtpDefaultValue::kInt;
   contPlaneCount100=ANtpDefaultValue::kInt;
   contPlaneCount200=ANtpDefaultValue::kInt;
   u_asym_peak=ANtpDefaultValue::kFloat;
   u_asym_vert=ANtpDefaultValue::kFloat;
   u_molrad_peak=ANtpDefaultValue::kFloat;
   u_molrad_vert=ANtpDefaultValue::kFloat;
   u_mean=ANtpDefaultValue::kFloat;
   u_rms= ANtpDefaultValue::kFloat;
   u_skew=ANtpDefaultValue::kFloat;
   u_kurt=ANtpDefaultValue::kFloat;
   v_asym_peak=ANtpDefaultValue::kFloat;
   v_asym_vert=ANtpDefaultValue::kFloat;
   v_molrad_peak=ANtpDefaultValue::kFloat;
   v_molrad_vert=ANtpDefaultValue::kFloat;
   v_mean=ANtpDefaultValue::kFloat;
   v_rms= ANtpDefaultValue::kFloat;
   v_skew=ANtpDefaultValue::kFloat;
   v_kurt=ANtpDefaultValue::kFloat;
   uv_asym_peak=ANtpDefaultValue::kFloat;
   uv_asym_vert=ANtpDefaultValue::kFloat;
   uv_molrad_peak=ANtpDefaultValue::kFloat;
   uv_molrad_vert=ANtpDefaultValue::kFloat;
   uv_mean= ANtpDefaultValue::kFloat;
   uv_rms = ANtpDefaultValue::kFloat;
   uv_skew= ANtpDefaultValue::kFloat;
   uv_kurt= ANtpDefaultValue::kFloat;
   uv_ratio=ANtpDefaultValue::kFloat;
   vtxEnergy=ANtpDefaultValue::kFloat;
   energyPlane0 = ANtpDefaultValue::kFloat;
   energyPlane1 = ANtpDefaultValue::kFloat;
   energyPlane2 = ANtpDefaultValue::kFloat;
   USlope = ANtpDefaultValue::kFloat;
   VSlope = ANtpDefaultValue::kFloat;
   UOffset = ANtpDefaultValue::kFloat;
   VOffset = ANtpDefaultValue::kFloat;
   UFitquality = ANtpDefaultValue::kFloat;
   VFitquality = ANtpDefaultValue::kFloat;
   UWDiff = ANtpDefaultValue::kFloat;
   VWDiff = ANtpDefaultValue::kFloat;
   USlopeMom = ANtpDefaultValue::kFloat;
   VSlopeMom = ANtpDefaultValue::kFloat;
   slopefix = ANtpDefaultValue::kFloat;
   UBeamLike = ANtpDefaultValue::kFloat;
   VBeamLike = ANtpDefaultValue::kFloat;
UVSlope = ANtpDefaultValue::kFloat;
ULongE = ANtpDefaultValue::kFloat;
VLongE = ANtpDefaultValue::kFloat;
LongE = ANtpDefaultValue::kFloat;
complexity =  ANtpDefaultValue::kFloat;
wcomplexity = ANtpDefaultValue::kFloat;

 u_asym_peak_9s_2pe = ANtpDefaultValue::kFloat;
   u_asym_vert_9s_2pe = ANtpDefaultValue::kFloat;
   u_molrad_peak_9s_2pe = ANtpDefaultValue::kFloat;
   u_molrad_vert_9s_2pe = ANtpDefaultValue::kFloat;
   u_mean_9s_2pe = ANtpDefaultValue::kFloat;
   u_rms_9s_2pe = ANtpDefaultValue::kFloat;
   u_skew_9s_2pe = ANtpDefaultValue::kFloat;
   u_kurt_9s_2pe = ANtpDefaultValue::kFloat;
   v_asym_peak_9s_2pe = ANtpDefaultValue::kFloat;
   v_asym_vert_9s_2pe = ANtpDefaultValue::kFloat;
   v_molrad_peak_9s_2pe = ANtpDefaultValue::kFloat;
   v_molrad_vert_9s_2pe = ANtpDefaultValue::kFloat;
   v_mean_9s_2pe = ANtpDefaultValue::kFloat;
   v_rms_9s_2pe = ANtpDefaultValue::kFloat;
   v_skew_9s_2pe = ANtpDefaultValue::kFloat;
   v_kurt_9s_2pe = ANtpDefaultValue::kFloat;
   uv_asym_peak_9s_2pe = ANtpDefaultValue::kFloat;
   uv_asym_vert_9s_2pe = ANtpDefaultValue::kFloat;
   uv_molrad_peak_9s_2pe = ANtpDefaultValue::kFloat;
   uv_molrad_vert_9s_2pe = ANtpDefaultValue::kFloat;
   uv_mean_9s_2pe = ANtpDefaultValue::kFloat;
   uv_rms_9s_2pe = ANtpDefaultValue::kFloat;
   uv_skew_9s_2pe = ANtpDefaultValue::kFloat;
   uv_kurt_9s_2pe = ANtpDefaultValue::kFloat;
   uv_ratio_9s_2pe = ANtpDefaultValue::kFloat;
                                                                                
   u_asym_peak_9s_2pe_dw = ANtpDefaultValue::kFloat;
   u_asym_vert_9s_2pe_dw = ANtpDefaultValue::kFloat;
   u_molrad_peak_9s_2pe_dw = ANtpDefaultValue::kFloat;
   u_molrad_vert_9s_2pe_dw = ANtpDefaultValue::kFloat;
   u_mean_9s_2pe_dw = ANtpDefaultValue::kFloat;
   u_rms_9s_2pe_dw = ANtpDefaultValue::kFloat;
   u_skew_9s_2pe_dw = ANtpDefaultValue::kFloat;
   u_kurt_9s_2pe_dw = ANtpDefaultValue::kFloat;
   v_asym_peak_9s_2pe_dw = ANtpDefaultValue::kFloat;
   v_asym_vert_9s_2pe_dw = ANtpDefaultValue::kFloat;
   v_molrad_peak_9s_2pe_dw = ANtpDefaultValue::kFloat;
   v_molrad_vert_9s_2pe_dw = ANtpDefaultValue::kFloat;
   v_mean_9s_2pe_dw = ANtpDefaultValue::kFloat;
   v_rms_9s_2pe_dw = ANtpDefaultValue::kFloat;
   v_skew_9s_2pe_dw = ANtpDefaultValue::kFloat;
   v_kurt_9s_2pe_dw = ANtpDefaultValue::kFloat;
   uv_asym_peak_9s_2pe_dw = ANtpDefaultValue::kFloat;
   uv_asym_vert_9s_2pe_dw = ANtpDefaultValue::kFloat;
   uv_molrad_peak_9s_2pe_dw = ANtpDefaultValue::kFloat;
   uv_molrad_vert_9s_2pe_dw = ANtpDefaultValue::kFloat;
   uv_mean_9s_2pe_dw = ANtpDefaultValue::kFloat;
   uv_rms_9s_2pe_dw = ANtpDefaultValue::kFloat;
   uv_skew_9s_2pe_dw = ANtpDefaultValue::kFloat;
   uv_kurt_9s_2pe_dw = ANtpDefaultValue::kFloat;
   uv_ratio_9s_2pe_dw = ANtpDefaultValue::kFloat;

/*
   if(lenepl==0){ lenepl=new TH1F;} 
   if(tenestu!=0){ tenestu=new TH1F(*(s.tenestu)); }
   if(tenestv!=0){ tenestv=new TH1F(*(s.tenestv)); }
   if(efit!=0){    efit=new TF1(*(s.efit)); }
   if(hfit!=0){    hfit=new TF1(*(s.hfit)); }
   if(info1!=0){ info1=new TPaveText(*(s.info1)); }
*/
}

