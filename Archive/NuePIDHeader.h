#ifndef NUEPIDHEADER_H
#define NUEPIDHEADER_H

#include "Validity/VldContext.h"
#include "Record/RecHeader.h"

class NuePIDHeader : public RecHeader
{
public:
   typedef enum EDecider {
      kUnknown = 0,
      kDT = 1, //mayly decision tree
      kMA = 2, //alex multivariate analysis
      kNN = 3, //tingjun's neural net
      kTH = 4, //truth (for testing)
      kCT = 5 //maylys cuts
   } Decider_t;

   static const char *AsString(Decider_t d){
      switch(d){
      case kUnknown: return "Unknown"; break;
      case kDT: return "Decision Tree"; break;
      case kMA: return "Multivariate Analysis"; break;
      case kNN: return "Neural Net"; break;
      case kTH: return "Truth"; break;
      case kCT: return "Cuts"; break;
      default: return "Unknown"; break;
      }
      return "Unknown";
   }

   NuePIDHeader();
   NuePIDHeader(const VldContext &vld);
   virtual ~NuePIDHeader();

   virtual std::ostream &Print(std::ostream &os) const;
   virtual void Print(Option_t *option="") const;

   const int GetSnarl() const;
   const int GetRun() const;
   const int GetSubRun() const;
   const int GetEventNo() const;
   const int GetEvents() const;
   const Decider_t GetDecider() const;

   void SetSnarl(int s);
   void SetRun(int r);
   void SetSubRun(int sr);
   void SetEventNo(int i);
   void SetEvents(int nevt);
   void SetDecider(Decider_t d);

private:
   int fSnarl;
   int fRun;
   int fSubRun;
   int fEvtNo;
   int fEvents;
   Decider_t fDecider;

   ClassDef(NuePIDHeader, 1)
};

inline const int NuePIDHeader::GetSnarl() const {return fSnarl;}
inline const int NuePIDHeader::GetRun() const {return fRun;}
inline const int NuePIDHeader::GetSubRun() const {return fSubRun;}
inline const int NuePIDHeader::GetEventNo() const {return fEvtNo;}
inline const int NuePIDHeader::GetEvents() const {return fEvents;}
inline const NuePIDHeader::Decider_t NuePIDHeader::GetDecider() const {return fDecider;}

inline void NuePIDHeader::SetSnarl(int s)  {fSnarl = s;}
inline void NuePIDHeader::SetRun(int r)  {fRun = r;}
inline void NuePIDHeader::SetSubRun(int sr)  {fSubRun = sr;}
inline void NuePIDHeader::SetEventNo(int i) {fEvtNo = i;}
inline void NuePIDHeader::SetEvents(int n) {fEvents = n;}
inline void NuePIDHeader::SetDecider(NuePIDHeader::Decider_t d) {fDecider = d;}


#endif  //NUEPIDHEADER_H
