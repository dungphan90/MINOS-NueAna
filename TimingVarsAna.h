#ifndef TIMINGVARSANA_H
#define TIMINGVARSANA_H

#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/Shwfit.h"
#include "Conventions/PlaneView.h"
#include "TimingVars.h"

class NtpSRRecord;

class TimingVarsAna : public NueAnaBase
{

public:
   TimingVarsAna(TimingVars &tv);
   virtual ~TimingVarsAna();
   
   void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);

    TH1F *TimeHstShw;
    TH1F *TimeHstTrk;
    TH1F *TotalHst;
//    TCanvas *plot;

private:

   Detector::Detector_t fDetectorType;
   TimingVars &fTimingVars;
};

#endif// TIMINGVARSANA_H
