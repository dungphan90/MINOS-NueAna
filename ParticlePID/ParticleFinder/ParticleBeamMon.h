#ifndef PARTICLEBEAMMON_H
#define PARTICLEBEAMMON_H

#include "TNamed.h"

#include "Conventions/BeamType.h"

//#ifndef RECRECORDIMP_H
#include "Record/RecRecordImp.h" // base class
//#endif
//#ifndef RECCANDHEADER_H
#include "Record/RecCandHeader.h" 
//#endif


#include "Validity/VldTimeStamp.h"

class ParticleBeamMon : public RecRecordImp<RecCandHeader>
{
 public:
    ParticleBeamMon();
    ParticleBeamMon(const RecCandHeader& header);
    virtual ~ParticleBeamMon();
  
    void ResetAll();
    void PrintInfo(Option_t *t) const;

    /// Returns the best value for pot
    double GetPot();
    

    /// Satisfies this spill the beam cuts?
    /// -1 means that it is not determined, 0 is bad, 1 is good
    int goodParticleBeamMon;   
    /// Filled with the best estimate of the pot, see ParticleBeamMon:GetPot()  
    double bI;
    /// Beam intensity from tortgt [1e12 ppp]
    double tortgt;
    /// Beam intensity from trtgtd [1e12 ppp]
    double trtgtd;
    /// Beam intensity from tor101 [1e12 ppp]
    double tor101;
    /// Beam intensity from tr101d [1e12 ppp]
    double tr101d;
    /// horiz beam pos from bpms
    double hpos2;
    /// vert beam pos from bpms
    double vpos2;
    /// beam position in X at target from BPM 121 and TGT
    /// for possibly 6 batches (mm)
    double batchposx[6];
    /// beam position in Y at target from BPM 121 and TGT
    /// for possibly 6 batches (mm)
    double batchposy[6]; 
    /// intensity in BPM TGT for possibly 6 batches
    double batchint[6];
    /// horiz beam pos from prof mons
    double hpos1;
    /// vert beam pos from prof mons
    double vpos1;	
    /// horiz beam width from prof mon
    double hbw;
    /// vert beam width from prof mon
    double vbw;		

    /// beam direction -- not filled
    double htan;
    /// beam direction -- not filled
    double vtan;	
    /// horn current
    double hornI;
    /// target pos -- not filled
    double nuTarZ;
    /// absolute beam monitor spill time -- not filled
    double time;
    /// absolute beam monitoring timestamp
    VldTimeStamp bmst_vts;
    /// absolute SpillTimeND time -- not filled
    double stnd_time;
    /// absolute SpillTimeND timestamp
    VldTimeStamp stnd_vts;
    /// detector time - BM Spill Time
    double dt_bmst;
    /// detector time - SpillTimeND 
    double dt_stnd;	

    /// Passes all of the NueStandard DataQuality Cuts
    int goodDataQual;  
    
    BeamType::BeamType_t beamtype;
    
    
private:
  ClassDef(ParticleBeamMon,1)
};
#endif //PARTICLEBEAMMON_H
