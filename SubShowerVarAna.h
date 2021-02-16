#ifndef SUBSHOWERVARANA_H
#define SUBSHOWERVARANA_H

#include "NueAna/NueAnaBase.h"
#include "CandNtupleSR/NtpSRCluster.h"
#include "NueAna/SubShowerVar.h"

class TPad;
class TNtuple;

class NtpSRRecord;
class NtpSRShower;

class SubShowerVarAna : public NueAnaBase
{

 public:

    SubShowerVarAna(SubShowerVar &sv);
    virtual ~SubShowerVarAna();

    void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
    void Analyze(RecRecordImp<RecCandHeader> *srobj, NtpSRShower* ntpShower);
    
 private:
  SubShowerVar &fSubShowerVar;
};

#endif// SUBSHOWERVARANA_H
