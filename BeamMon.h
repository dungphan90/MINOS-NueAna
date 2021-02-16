#ifndef BEAMMON_H
#define BEAMMON_H
#include "TObject.h"
#include "Validity/VldTimeStamp.h"

class BeamMon : public TObject
{
 public:
    BeamMon();
    virtual ~BeamMon();
  
    void Reset();
    void Print(Option_t *t) const;

    /// Returns the best value for pot
    double GetPot();
    

    /// Satisfies this spill the beam cuts?
    /// -1 means that it is not determined, 0 is bad, 1 is good
    int goodBeamMon;   
    /// Filled with the best estimate of the pot, see BeamMon:GetPot()  
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
    
    
private:
  ClassDef(BeamMon,4)
};
#endif //BEAMMON_H
