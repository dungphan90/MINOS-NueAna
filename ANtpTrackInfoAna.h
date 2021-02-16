#ifndef ANTPTRACKINFOANA_H
#define ANTPTRACKINFOANA_H

#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/ANtpTrackInfoNue.h"
#include "AnalysisNtuples/Module/ANtpInfoObjectFiller.h"
#include "DataUtil/EnergyCorrections.h"

class NtpSRRecord;

class ANtpTrackInfoAna : public NueAnaBase
{

public:
   ANtpTrackInfoAna(ANtpTrackInfoNue &anti);
   virtual ~ANtpTrackInfoAna();

   //   void Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord *mc=0, NtpTHRecord *th=0);
   void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);

   void FillNueTrackInformation(NtpSRTrack *ntpTrack, NtpSREvent *ntpEvent, ANtpTrackInfoNue *trackInfoNue);

   void DetermineSigInOut(NtpSRTrack *ntpTrack, RecRecordImp<RecCandHeader> *srobj);
   Bool_t IsFidAll(NtpSRTrack *ntpTrack);
   Float_t RecoMuEnergy(SimFlag::SimFlag_t s, const Detector::Detector_t det);
   Float_t RecoMuEnergyNew(VldContext cx, EnergyCorrections::WhichCorrection_t corrver = EnergyCorrections::kDefault);


private:
     Int_t fDetectorType;
  
    ANtpTrackInfoNue &fANtpTrackInfo;
    ANtpInfoObjectFiller *fInfoFiller;
};

#endif// ANTPTRACKINFOANA_H
