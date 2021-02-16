#ifndef NUERECORDANA_H
#define NUERECORDANA_H
#include "NueAna/ShwfitAna.h"
#include "NueAna/HitCalcAna.h"
#include "NueAna/AngClusterAna.h"
#include "NueAna/AngClusterFitAna.h"
#include "NueAna/MdaDiscrimAna.h"
#include "NueAna/NueRecord.h"
#include "NueAna/AnalysisInfoAna.h"
#include "NueAna/ANtpEventInfoAna.h"
#include "NueAna/ANtpShowerInfoAna.h"
#include "NueAna/ANtpTrackInfoAna.h"
#include "NueAna/ANtpTruthInfoBeamAna.h"
#include "NueAna/MSTCalcAna.h"
#include "NueAna/FracVarAna.h"
#include "NueAna/SubShowerVarAna.h"
#include "NueAna/HighHitVarsAna.h"
#include "NueAna/ShieldRejVarAna.h"
#include "NueAna/AnnAna.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/NueRecord.h"
#include "NueAna/TreePIDAna.h"
#include "NueAna/MCFluxInfoAna.h"
#include "NueAna/NueFluxWeightsAna.h"
#include "NueAna/NueXsecWeightAna.h"
#include "NueAna/StdHepInfoAna.h"
#include "NueAna/MuonRemovalInfoAna.h"
#include "NueAna/CalDetInfoAna.h"
#include "CalDetDST/UberRecord.h"
#include "NueAna/TimingVarsAna.h"
#include "NueAna/BagTreeAna.h"
#include "NueAna/EventQualAna.h"

class NueRecordAna: public NueAnaBase
{
 public:
  NueRecordAna(NueRecord &nr);
  virtual ~NueRecordAna();

  typedef std::deque<Float_t> DeqFloat_t; 
  typedef std::deque< std::deque <Int_t> > DeqDeqInt_t; 
  
  //   void Analyze(int evtn, NtpSRRecord *srobj, 
  //		NtpMCRecord *mcobj=0, NtpTHRecord *thobj=0);
  void FillTrue(int evtn, NtpSRRecord *srobj, 
		NtpMCRecord *mcobj=0, NtpTHRecord *thobj=0);

  void FillTrue(int evtn, NtpStRecord *str);
  void FillReco(int evtn, RecRecordImp<RecCandHeader> *srobj);

  void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);

  //for filling muon removal info:
  void Analyze(int evtn, NtpMRRecord *mr,NtpStRecord *oldst); 
  void Analyze(UberRecord *ur); //for filling CalDet info

  void SetEventEnergyArray(float* array, float* arr2);
  void SetRelease(int rel);

  void SetBeamType(int bType);
  void SetBeamType(BeamType::BeamType_t bType);


   ShwfitAna sfa;
   HitCalcAna hca;
   AngClusterAna aca;
   AngClusterFitAna acf;
   MSTCalcAna msta;
   FracVarAna fva;
   SubShowerVarAna sva;
   HighHitVarsAna hhv;
   ShieldRejVarAna shrva;
   AnnAna anna;
   AnalysisInfoAna anaia;
   ANtpEventInfoAna aneia;
   ANtpShowerInfoAna ansia;
   ANtpTrackInfoAna antia;
   ANtpTruthInfoBeamAna antiba;
//   BeamMonAna bmona;
   MdaDiscrimAna mda;
   TreePIDAna tpa;
   MCFluxInfoAna fiana;
   NueFluxWeightsAna nuefwa;
   NueXsecWeightAna nuexsa;
   StdHepInfoAna shia;
   MuonRemovalInfoAna mria;
   CalDetInfoAna cdia;
   TimingVarsAna tva;
   BagTreeAna bta;
   EventQualAna equala;

  private:
   bool active3DHits;

};

#endif //NUERECORDANA_H
