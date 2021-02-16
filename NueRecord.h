#ifndef NUERECORD_H
#define NUERECORD_H

#include "Record/RecRecordImp.h"
#include "Record/RecHeader.h"
#include "NueAna/NueHeader.h"
#include "NueAna/Shwfit.h"
#include "NueAna/HitCalc.h"
#include "NueAna/AngCluster.h"
#include "NueAna/AngClusterFit.h"
#include "NueAna/MdaDiscrim.h"
#include "NueAna/MSTCalc.h"
#include "NueAna/FracVar.h"
#include "NueAna/SubShowerVar.h"
#include "NueAna/HighHitVars.h"
#include "NueAna/ShieldRejVar.h"
#include "NueAna/Ann.h"
#include "NueAna/AnalysisInfoNue.h"
#include "NueAna/ANtpEventInfoNue.h"
#include "NueAna/ANtpTrackInfoNue.h"
#include "NueAna/ANtpShowerInfoNue.h"
#include "NueAna/BeamMon.h"
#include "NueAna/ANtpTruthInfoBeamNue.h"
#include "NueAna/TreePID.h"
#include "NueAna/NueFluxWeights.h"
#include "NueAna/NueXsecWeight.h"
#include "MCNtuple/NtpMCFluxInfo.h"
#include "NueAna/StdHepInfo.h"
#include "NueAna/MuonRemovalInfo.h"
#include "NueAna/CalDetInfo.h"
#include "NueAna/TimingVars.h"
#include "NueAna/MCNNVars.h"
#include "NueAna/BagTree.h"
#include "NueAna/EventQual.h"
#include "NueAna/ParticlePID/ParticleAna/PRecord.h"


class NtpSRRecord;
class NtpMCRecord;
class NtpTHRecord;


class NueRecord : public RecRecordImp<NueHeader>
{
public:
   NueRecord();
   NueRecord(const NueHeader& head);
   NueRecord(const NueRecord& nr);
   virtual ~NueRecord();

   void Reset();
   void Clear(Option_t* /* option */);

   Shwfit shwfit;
   HitCalc hitcalc; 
   AngCluster angcluster;
   AngClusterFit angclusterfit;
//   RecoCalc reco;
//   VertexFinder vertfind;
//   VertexInfo vtx;
//   EmShw emvars;
   MSTCalc mstvars;
   FracVar fracvars;
   SubShowerVar subshowervars;
   HighHitVars highhitvars;
   ShieldRejVar shieldrejvars;
   Ann ann;
   AnalysisInfoNue anainfo;
   ANtpEventInfoNue srevent;
   ANtpShowerInfoNue srshower;
   ANtpTrackInfoNue srtrack;
   ANtpTruthInfoBeamNue mctrue;
   BeamMon bmon;
   MdaDiscrim mdadiscrim;
   TreePID treepid;
   NtpMCFluxInfo fluxinfo;
   NueFluxWeights fluxweights;
   NueXsecWeight xsecweights;
   StdHepInfo shi;
   MuonRemovalInfo mri;
   CalDetInfo cdi;
   TimingVars timingvars;
   MCNNVars mcnnv;
   BagTree dtree;
   EventQual eventq;

   PRecord precord;
   PRecord precordMRCC;
 
private:
   ClassDef(NueRecord,20)
};

#endif // NUERECORD_H
