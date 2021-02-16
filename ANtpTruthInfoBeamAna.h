#ifndef ANTPTRUTHINFOBEAMANA_H
#define ANTPTRUTHINFOBEAMANA_H

#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/ANtpTruthInfoBeamNue.h"
#include "AnalysisNtuples/Module/ANtpInfoObjectFillerBeam.h"


class NtpSRRecord;
class NtpMCRecord;
class NtpTHRecord;
class NtpStRecord;

class ANtpTruthInfoBeamAna : public NueAnaBase
{

public:
    ANtpTruthInfoBeamAna(ANtpTruthInfoBeamNue &antib);
    virtual ~ANtpTruthInfoBeamAna();


    //need two analyze functions, one for NtpStRecord, and one for 
    // when the three individual objects are read.
    void Analyze(int evtn, NtpStRecord *srobj);
    void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
    void Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord *mcobj, NtpTHRecord *thobj);

    Int_t GetNueClass(Int_t inu, Int_t inunoosc, Int_t iaction);
    Float_t GetNueWeight(Int_t inu, Int_t inunoosc);
    Float_t GetOscProb();

    Float_t TrueLepDCosNeu(NtpMCTruth *ntpTruth);
    Float_t TrueLepDCosZ(NtpMCTruth *ntpTruth);
    Float_t Get3Momenta(NtpMCTruth *ntpTruth,
                           Float_t &p4_0, Float_t &p4_1, Float_t &p4_2);


private:
    ANtpTruthInfoBeamNue &fANtpTruthInfoBeam;
    ANtpInfoObjectFillerBeam *fInfoFiller;

};

#endif// ANTPTRUTHINFOBEAMANA_H
