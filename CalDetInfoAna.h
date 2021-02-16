#ifndef CALDETINFOANA_H
#define CALDETINFOANA_H

#include "NueAna/NueAnaBase.h"
#include "NueAna/CalDetInfo.h"
#include "CalDetDST/UberRecHeader.h"

class CalDetInfoAna : public NueAnaBase
{

 public:
  
    CalDetInfoAna(CalDetInfo &sv);
    virtual ~CalDetInfoAna();
    void Analyze(int /*evtn*/, RecRecordImp<RecCandHeader> */*srobj*/) {};
    void Analyze(RecRecordImp<UberRecHeader> *uberrecord);
    
 private:
  CalDetInfo &fCalDetInfo;
};

#endif// CALDETINFOANA_H
