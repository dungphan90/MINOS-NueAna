#ifndef BEAMMONANA_H
#define BEAMMONANA_H

#include <vector>
#include <set>
#include "NueAna/NueAnaBase.h"
#include "NueAna/BeamMon.h"

class BeamSummary;

using std::vector;
using std::set;

class BeamMonAna: public NueAnaBase
{

public:
   BeamMonAna(BeamMon &bmon);
   virtual ~BeamMonAna();

   void SetBeamSummary(BeamSummary *b);
   //   void Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord *mc=0, NtpTHRecord *th=0);   
   void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);

private:
   BeamMon &fBmon;
   BeamSummary *bs;
};
#endif//BEAMMONANA_H
