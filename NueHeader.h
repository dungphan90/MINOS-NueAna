#ifndef NUEHEADER_H
#define NUEHEADER_H

#include "Validity/VldContext.h"
#include "Record/RecHeader.h"
#include "Conventions/ReleaseType.h"
#include "Conventions/BeamType.h"

class NueHeader : public RecHeader
{
public:
   NueHeader();
   NueHeader(const VldContext &vld);
   virtual ~NueHeader();

   virtual std::ostream &Print(std::ostream &os) const;
   virtual void Print(Option_t *option="") const;

   int GetSnarl() const;
   int GetRun() const;
   int GetSubRun() const;
   int GetEventNo() const;
   int GetEvents() const;
   int GetTrackLength() const;
   ReleaseType::Release_t GetRelease() const;
   BeamType::BeamType_t GetBeamType() const;
   int GetNueRelease() const;
   bool FoundSR() const;
   bool FoundMC() const;
   bool FoundTH() const;
   bool FoundMR() const;

   void SetSnarl(int s);
   void SetRun(int r);
   void SetSubRun(int sr);
   void SetEventNo(int i);
   void SetEvents(int nevt);
   void SetTrackLength(int t);
   void SetFoundBits(bool SR, bool MC, bool TH, bool MR);
   void SetRelease(int r);
   void SetBeamType(int b);
   void SetNueRelease(int nr);

private:
   int fSnarl;
   int fRun;
   int fSubRun;
   int fEvtNo;
   int fEvents;
   int fTrackLength;
   ReleaseType::Release_t fRelease;
   BeamType::BeamType_t fBeamType;   
   int fNueRelease;

   bool foundSR;
   bool foundMC;
   bool foundTH;
   bool foundMR;

   ClassDef(NueHeader, 5)
};

inline int NueHeader::GetSnarl() const {return fSnarl;}
inline int NueHeader::GetRun() const {return fRun;}
inline int NueHeader::GetSubRun() const {return fSubRun;}
inline int NueHeader::GetEventNo() const {return fEvtNo;}
inline int NueHeader::GetEvents() const {return fEvents;}
inline bool NueHeader::FoundSR() const {return foundSR;}
inline bool NueHeader::FoundMC() const {return foundMC;}
inline bool NueHeader::FoundTH() const {return foundTH;}
inline bool NueHeader::FoundMR() const {return foundMR;}
inline int NueHeader::GetTrackLength() const {return fTrackLength;}
inline ReleaseType::Release_t NueHeader::GetRelease() const {return fRelease;}
inline BeamType::BeamType_t NueHeader::GetBeamType() const {return fBeamType;}
inline int NueHeader::GetNueRelease() const {return fNueRelease;}


inline void NueHeader::SetSnarl(int s)  {fSnarl = s;}
inline void NueHeader::SetRun(int r)  {fRun = r;}
inline void NueHeader::SetSubRun(int sr)  {fSubRun = sr;}
inline void NueHeader::SetEventNo(int i) {fEvtNo = i;}
inline void NueHeader::SetEvents(int n) {fEvents = n;}
inline void NueHeader::SetTrackLength(int t){fTrackLength=t;}
inline void NueHeader::SetFoundBits(bool SR, bool MC, bool TH, bool MR){foundSR=SR;foundMC=MC;foundTH=TH;foundMR=MR;}
inline void NueHeader::SetRelease(int r) {fRelease = r;}
inline void NueHeader::SetBeamType(int b) {fBeamType = (BeamType::BeamType_t) b;}
inline void NueHeader::SetNueRelease(int r) {fNueRelease = r;}

#endif  //NUEHEADER_H
