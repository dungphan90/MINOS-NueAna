#ifndef NUEANABASE_H
#define NUEANABASE_H
#include "Record/RecCandHeader.h"
#include "Record/RecRecordImp.h"
#include "NueAnaTools/NueConvention.h"
#include "Conventions/ReleaseType.h"
#include "Conventions/BeamType.h"

class NtpSRRecord;
class NtpMCRecord;
class NtpTHRecord;

class NueAnaBase{

public:
   NueAnaBase(): sigcormeu(1.), MeuPerGeV(1.),
                 release(ReleaseType::kUnknown),
                 evtstp0mip(0), evtstp1mip(0) {};

   virtual ~NueAnaBase() {};
   // Sub class overrides to fill ntuple data
   //   virtual void Analyze(int evtn, NtpSRRecord *srobj=0,
   //		NtpMCRecord *mcobj=0,
   //		NtpTHRecord *thobj=0) = 0;
   virtual void Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)=0;
   void SetParams(float scm) {sigcormeu = scm; }
   void SetParams(float scm, float mpg) {sigcormeu = scm; MeuPerGeV = mpg;}
   void SetRelease(int rel) {release = rel;}
   void SetBeamType(int type) { beam = (BeamType::BeamType_t) type; };
   void SetBeamType(BeamType::BeamType_t type) { beam = type; };


   void SetEventEnergyArray(float* ph0, float* ph1) {evtstp0mip = ph0; evtstp1mip = ph1;}
   float sigcormeu;
   float MeuPerGeV;
   ReleaseType::Release_t release;
   BeamType::BeamType_t beam;

   float* evtstp0mip;
   float* evtstp1mip;
 
protected:
//   float sigcormeu;
      
};
#endif //NUEANABASE_H
