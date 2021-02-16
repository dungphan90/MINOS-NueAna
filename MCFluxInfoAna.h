#ifndef MCFLUXINFOANA_H
#define MCFLUXINFOANA_H
#include <string>
#include "NueAna/NueAnaBase.h"
#include "MCNtuple/NtpMCFluxInfo.h"

class NtpMCTruth;
class NtpStRecord;
class MuParentHelper;

class MCFluxInfoAna: public NueAnaBase
{

public:
   MCFluxInfoAna(NtpMCFluxInfo &fluxinfo);
   virtual ~MCFluxInfoAna();

   void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj);
   void Analyze(NtpMCTruth *mcrec);

   void ResetMuParentInfo();

   NtpMCFluxInfo* GetMCFluxInfo() {return &fMCFluxInfo;}

   void SetMuParentHelper(MuParentHelper *mp) {mupar=mp;}
   void FixMuParents(bool debugger) {fluxfileneedsdebuggering=debugger;}
private:
   NtpMCFluxInfo &fMCFluxInfo;
   bool fluxfileneedsdebuggering;
   MuParentHelper *mupar;

};
#endif//MCFLUXINFOANA_H
