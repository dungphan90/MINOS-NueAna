#ifndef SHIELDREJVARANA_H
#define SHIELDREJVARANA_H

#include "NueAna/NueAnaBase.h"
#include "NueAna/ShieldRejVar.h"
#include "CandShield/ShieldGeom.h"

class TNtuple;

class NtpSRRecord;

class ShieldRejVarAna : public NueAnaBase
{

 public:

  ShieldRejVarAna(ShieldRejVar &sv);
  virtual ~ShieldRejVarAna();

  void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);

 private:
  ShieldRejVar &fShieldRejVar;


};

#endif// SHIELDREJVARANA_H
