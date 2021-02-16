#ifndef STDHEPINFOANA_H
#define STDHEPINFOANA_H

#include "NueAnaBase.h"
#include "StdHepInfo.h"

class StdHepInfoAna : public NueAnaBase
{

 public:

  StdHepInfoAna(StdHepInfo &sv);
  virtual ~StdHepInfoAna();

  void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
  void Analyze(int evtn, NtpStRecord *srobj);

 private:
  StdHepInfo &fStdHepInfo;
};

#endif// STDHEPINFOANA_H
