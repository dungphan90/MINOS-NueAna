#ifndef HIGHHITVARSANA_H
#define HIGHHITVARSANA_H

#include "NueAna/NueAnaBase.h"
#include "NueAna/HighHitVars.h"

class TPad;
class TNtuple;

class NtpSRRecord;
class NtpSREvent;
class HighHitVarsAna : public NueAnaBase

{

 public:

  HighHitVarsAna(HighHitVars &hhv);
  virtual ~HighHitVarsAna();

  void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
     void Analyze(RecRecordImp<RecCandHeader> *srobj, NtpSRShower* ntpShower);

 private:
  HighHitVars &fHighHitVars;
 
};

#endif// HIGHHITVARSANA_H
