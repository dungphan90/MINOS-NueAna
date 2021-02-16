#ifndef ANTPEVENTINFOANA_H
#define ANTPEVENTINFOANA_H

#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/ANtpEventInfoNue.h"
#include "AnalysisNtuples/Module/ANtpInfoObjectFiller.h"


class NtpSRRecord;

class ANtpEventInfoAna : public NueAnaBase
{

public:
    ANtpEventInfoAna(ANtpEventInfoNue &anei);
    virtual ~ANtpEventInfoAna();

    //    void Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord *mc=0, NtpTHRecord *th=0);
    void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
    void FillNueEventInformation(NtpSREvent *ntpEvent, ANtpEventInfoNue *eventInfoNue);
    Float_t GetNDCoilCurrent(const VldContext& vc);
    Int_t FDRCBoundary(NtpSREventSummary *eventSummary);

    void FillStripVariables(NtpSREvent *ntpEvent, RecRecordImp<RecCandHeader> *srobj);
 
 
private:
    ANtpEventInfoNue &fANtpEventInfo;
    ANtpInfoObjectFiller *fInfoFiller;

    Int_t fDetectorType;
};

#endif// ANTPEVENTINFOANA_H
