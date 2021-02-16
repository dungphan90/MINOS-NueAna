/**
 *
 * $Id: AnalysisInfoAna.h,v 1.9 2008/11/19 18:22:51 rhatcher Exp $
 *
 * \class AnalysisInfoAna
 *
 * \package NueAna
 *
 * \brief
 *
 * Contact: Josh Boehm
 *
 * Created on: Mon Apr 18 16:59:24 2005
 *
 */
                                                                                
#ifndef ANALYSISINFOANA_H
#define ANALYSISINFOANA_H
      
#include <deque>
#include <vector>
#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/AnalysisInfoNue.h"
#include "Mad/MadNsID.h"
#include "Mad/MadDpID.h"
#include "Mad/MadAbID.h"
#include "PhysicsNtuple/Store/Interface.h"                                                                         
class NtpSRRecord;
//class NtpSRShower;
//class NtpSRTrack;
//class NtpSREvent;
                                                                                
class AnalysisInfoAna : public NueAnaBase
{
                                                                                
public:
    AnalysisInfoAna(AnalysisInfoNue &anai);
    ~AnalysisInfoAna();

    void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj); 

    //    void FillNueAnalysisInformation(NtpSREvent *ntpEvent, NtpSRTrack *ntpTrack,
    //           NtpSRShower *ntpShower, AnalysisInfoNue *analysisInfoNue);
                                                
    typedef std::deque<Float_t> DeqFloat_t;
    typedef DeqFloat_t::iterator IterDeqFloat_t;
                                                                                
    void Set3DHit(DeqFloat_t &x, DeqFloat_t &y, DeqFloat_t &z, DeqFloat_t &e);
    //Int_t IsFidVtxEvt(NtpSREvent *ntpEvent, Int_t detType);
                                                                                
private:
                                                                                
    AnalysisInfoNue &fAnalysisInfo;
    Int_t IsFidAll(Float_t vtxX, Float_t vtxY, Float_t vtxZ, NtpSREvent *event = 0);

    string BuildABPIDFile();
    string BuildROPIDFile();
    DeqFloat_t fX;
    DeqFloat_t fY;
    DeqFloat_t fZ;
    DeqFloat_t fE;
                                                                                
    Detector::Detector_t fDetectorType;
//    BeamType::BeamType_t fBeam;
                                                                                
    static MadNsID nsid;
    static MadDpID dpid;
    static MadAbID abid;
// Note that AnalysisInfo does not have a filler object in the AnalysisPackage at this moment

    static bool readabidfile;
};       // end of class AnalysisInfoAna
                                                                                
#endif  // ANALYSISINFOANA_H

