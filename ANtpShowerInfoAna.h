#ifndef ANTPSHOWERINFOANA_H
#define ANTPSHOWERINFOANA_H

#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/ANtpShowerInfoNue.h"
#include "AnalysisNtuples/Module/ANtpInfoObjectFiller.h"
#include "DataUtil/EnergyCorrections.h"

class NtpSRRecord;

class ANtpShowerInfoAna : public NueAnaBase
{

public:
   ANtpShowerInfoAna(ANtpShowerInfoNue &ansi);
   virtual ~ANtpShowerInfoAna();

   //   void Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord *mc=0, NtpTHRecord *th=0);
   void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
   void FillNueShowerInformation(NtpSRShower *ntpShower, NtpSREvent *ntpEvent, ANtpShowerInfoNue *showerInfoNue,RecRecordImp<RecCandHeader> *srobj);

    Float_t RecoShwEnergy(Float_t linearCCGeV, SimFlag::SimFlag_t s, const Detector::Detector_t& det, int mode);
    Float_t RecoShwEnergy(NtpSRShower * ntpShower, Int_t opt, Int_t det);
    Float_t RecoShwEnergyNew(NtpSRShower * ntpShower,
            Int_t opt, VldContext cx);

    Float_t GetShwEnergy(NtpSRShower* ntpShower, Int_t opt);    
    CandShowerHandle::EShowerType GetShwHandleType(Int_t opt);

    void FillGapInformation(NtpSREvent *evt, NtpSRShower *shw, NtpSRTrack* trk, NtpStRecord * st);

private:
    ANtpShowerInfoNue &fANtpShowerInfo;
    ANtpInfoObjectFiller *fInfoFiller;
};

#endif// ANTPSHOWERINFOANA_H
