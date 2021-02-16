#ifndef SHWFITANA_H
#define SHWFITANA_H

#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/Shwfit.h"
#include "Conventions/PlaneView.h"

class NtpSRRecord;

class ShwfitAna : public NueAnaBase
{

public:
   ShwfitAna(Shwfit &sf);
   virtual ~ShwfitAna();

   //   void Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord *mc=0, NtpTHRecord *th=0);

   void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
   void doSlopes( NtpSREvent *event , RecRecordImp<RecCandHeader> *srobj);
   void FitLShower(Float_t pulseheight);
   void FitLShower_Dan(Float_t pulseheight);
   void TransVar(TH1F *h, PlaneView::EPlaneView pv);
   void Reset(int snarl, int event);
   void SetCutParams(int planes, float striphcut, float planephcut, float contplanephcut){sfDPlaneCut=planes; sfPhStripCut=striphcut; sfPhPlaneCut=planephcut; sfContPhPlaneCut=contplanephcut; }
   Bool_t PassCuts(int PhNStrips, int PhNPlanes );
   int sfDPlaneCut;

   void FitTShower(Float_t pulseheight);



   //  contPlaneCount
   float sfContPhPlaneCut;

   float sfPhStripCut;
   float sfPhPlaneCut;
private:
   Shwfit &fShwfit;
   Double_t GetMaximumX(TF1* efit, Double_t xmin=0, Double_t xmax=0);
   float BuildUVVar(float u, float v);

   float asym_peak;
   float asym_vert;
   float molrad_peak;
   float molrad_vert;
   float mean;
   float rms;
   float skew;
   float kurt;
};

#endif// SHWFITANA_H
