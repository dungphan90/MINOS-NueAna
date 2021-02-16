#ifndef EVENTQUALANA_H
#define EVENTQUALANA_H

#include "TObject.h"
#include "NueAna/NueAnaBase.h"
#include "NueAna/EventQual.h"

class NtpSRRecord;
class NueRecord;

class EventQualAna : public NueAnaBase
{

public:
    
    EventQualAna(NueRecord &nr);
    ~EventQualAna();

    void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
    Int_t FDRCBoundary(NtpSREventSummary *eventSummary);
 
private:
    NueRecord &nr;
    EventQualAna(const EventQualAna& rhs);

    ClassDef(EventQualAna, 1)
};

#endif  // EVENTQUALANA_H
