/**
 *
 * $Id: ANtpAnalysisInfoAna.h,v 1.2 2008/11/19 18:22:51 rhatcher Exp $
 *
 * \class ANtpAnalysisInfoAna
 *
 * \package NueAna
 *
 * \brief 
 *
 * Contact: M. Sanchez
 *
 * Created on: Mon Apr 18 16:59:24 2005
 *
 */

#ifndef ANTPANALYSISINFOANA_H
#define ANTPANALYSISINFOANA_H

                                                                                
#include <deque>
#include <vector>
#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/ANtpAnalysisInfoNue.h"
#include "AnalysisNtuples/Module/ANtpInfoObjectFiller.h"
#include "AnalysisNtuples/Module/ANtpRecoNtpManipulator.h"
#include "Mad/MadNsID.h"
#include "Mad/MadDpID.h"

class NtpSRRecord;
//class NtpSRShower;
//class NtpSRTrack;
//class NtpSREvent;

class ANtpAnalysisInfoAna : public NueAnaBase
{

public:
    ANtpAnalysisInfoAna(ANtpAnalysisInfoNue &anai);
    ~ANtpAnalysisInfoAna();

    void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
	
//    void Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord *mc=0, NtpTHRecord *th=0);
    void FillNueAnalysisInformation(NtpSREvent *ntpEvent, NtpSRTrack *ntpTrack,                                   NtpSRShower *ntpShower, ANtpAnalysisInfoNue *analysisInfoNue, ANtpInfoObjectFiller * filla);
    Float_t RecoMuEnergy(NtpSRTrack *ntpTrack);
    Float_t RecoShwEnergy(NtpSRShower *ntpShower);

    typedef std::deque<Float_t> DeqFloat_t;
    typedef DeqFloat_t::iterator IterDeqFloat_t;   

    void Set3DHit(DeqFloat_t &x, DeqFloat_t &y, DeqFloat_t &z, DeqFloat_t &e);
    Int_t IsFidVtxEvt(NtpSREvent *ntpEvent, Int_t detType);

private:

    ANtpAnalysisInfoNue &fANtpAnalysisInfo;

    Float_t RecoQENuEnergy(NtpSRTrack *ntpTrack);
    Float_t RecoQEQ2(NtpSRTrack *ntpTrack);
    Float_t RecoMuDCosZ(NtpSRTrack *ntpTrack);
    Float_t RecoMuDCosNeu(NtpSRTrack *ntpTrack);
    Float_t RecoMuQP(NtpSRTrack *ntpTrack);
    Int_t  GetChargeSign(NtpSRTrack *ntpTrack);
    Bool_t IsFidVtx(NtpSRTrack *ntpTrack);
    Int_t IsFidAll(Float_t vtxX, Float_t vtxY, Float_t vtxZ, NtpSREvent *event = 0);
    Int_t IsFidVtxEvt(NtpSREvent *ntpEvent);
    Int_t InsideNearFiducial(Float_t x, Float_t y, Float_t z);
    Int_t InsideFarFiducial(Float_t x, Float_t y, Float_t z);
	    
    DeqFloat_t fX;
    DeqFloat_t fY;
    DeqFloat_t fZ;
    DeqFloat_t fE;

    Detector::Detector_t fDetectorType;

    static MadNsID nsid;
    static MadDpID dpid;	    
// Note that AnalysisInfo does not have a filler object in the AnalysisPackage at this moment 
//    ANtpInfoObjectFiller *fInfoFiller;   

};                              // end of class ANtpAnalysisInfoAna

#endif  // ANTPANALYSISINFOANA_H
