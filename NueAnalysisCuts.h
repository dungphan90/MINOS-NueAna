#ifndef NUEANALYSISCUTS_H
#define NUEANALYSISCUTS_H
                                                                                           
#include "CandNtupleSR/NtpSREvent.h"
#include "StandardNtuple/NtpStRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include "MCNtuple/NtpMCStdHep.h"
#include "TruthHelperNtuple/NtpTHEvent.h"
#include "AnalysisNtuples/Module/ANtpRecoNtpManipulator.h"
#include "NueAna/NueRecord.h"
#include <vector>
#include "Registry/Registry.h"
#include "NueAna/NueAnaTools/NueConvention.h"

using namespace std;

class NtpStRecord;

namespace FiducialAlg
{
    static const int kUseEvtVtx = 1;
    static const int kUseTrkVtx = 2;                                                                                  
    static const int kDefaultVolume = 1;
    static const int kExtendedVolume = 2;
}

class NueAnalysisCuts : public TObject
{
   public:
     NueAnalysisCuts();
     virtual ~NueAnalysisCuts() {};

     void Reset();
     const Registry& DefaultConfig() const;
     void Config(const Registry &r);
//     void LoadDefaultRegistry(Registry &r) const;
//     void SetRegistry(const Registry &r);

     void SetInfoObject(NueRecord* nr);
     void SetInfoObject(int evtn, NtpStRecord* st);    
     void SetInfoObject(int evtn, NtpSRRecord* st, NtpMCRecord *mcobj=0, NtpTHRecord * thobj=0);
     void SetNtpInfoObjects(int evtn, int run, VldContext &vc);

     void Report();  //Sends a report of set cut values to MSGServernalysisCuts::SetNtpInfoObjects(ANtpRecoNtpManipulator &ntpManipulator, int evtn)
     
     // Conditional Checks
     bool PassesTrackPlaneCut();
     bool PassesTrackLikePlaneCut();
     bool PassesHighShowerEnergyCut();
     bool PassesLowShowerEnergyCut();
     bool PassesHighEnergyCut();
     bool PassesLowEnergyCut();
     bool PassesEventPlaneCut();
     bool PassesHotChannel();
     bool PassesFiducialVolume();
     bool PassesFullContainment();
     bool PassesPhProngCut();
     bool PassesCutOnClasses();
     bool IsGoodBeam();
     bool PassesBeamCut() {return IsGoodBeam();};
     bool PassesHornCurrent();   
     bool PassesDataQuality();
  
     bool PassesLowPhNPlaneCut();
     bool PassesLowPhNStripCut();
     bool PassesEMFraction();
     bool PassesResonanceCode();
     bool PassesFileCut();

     bool PassesAllCuts();

   private:
     bool IsValid();
     bool IsMC();

     void FillVertexPosition(float &x, float &y, float &z);
     bool IsInsideFarFiducial(float x, float y, float z);
     bool IsInsideNearFiducial(float x, float y, float z);
//     int IsInsideNearFiducial_Extended(float x, float y, float z);
//     int IsInsideFarFiducial_Extended(float x, float y, float z);

     bool PassesLowEMFraction();
     bool PassesHighEMFraction();

     bool IsSignal();
     bool IsBackground();

     void ReportOnRecord(NueRecord *nr, string ID);
     bool NeedsSpecialAttention(TString name);

     NueRecord LowCuts;
     NueRecord HighCuts;

     float kMEUPerGeV;
     float kSigCorrPerMEU;
     Detector::Detector_t fDetectorType;
     SimFlag::SimFlag_t fSimFlag;

     int fFiducialCut;
     int fFidVolumeChoice;
     int fFidVtxChoice;
     int fContainmentCut;

     int fBeamCut;
     int fCheckHotChannel;
     float fTargetCurrent;
     float fPhProngCut;

     int fFileTrim;
     std::string fFileCutName;

     std::string kOutputFile;

     int fCutClasses;
     vector<int> signalClasses;
     vector<int> bgClasses;

     int fCutResCode;
     vector<int> sigResonanceCodes;
     vector<int> bgResonanceCodes;  
     ANtpTruthInfoBeam HighSigCuts;
     ANtpTruthInfoBeam LowSigCuts;
     ANtpTruthInfoBeam HighBgCuts;
     ANtpTruthInfoBeam LowBgCuts;



    // Variables 
     NueRecord* nrInfo;  
     NtpStRecord* stInfo;
     NtpSRRecord* srInfo;

     ANtpRecoNtpManipulator ntpManipulator;
     NtpSREvent* eventInfo;
     NtpSRTrack* trackInfo;
     NtpSRShower* showerInfo;
     NtpMCTruth *mcInfo;
//     NtpMCStdHep *mcstdhep;
     NtpTHEvent *thEventInfo;

     int fEventNum;

     ClassDef(NueAnalysisCuts,4)
};

#endif
