#include <iostream>
#include "TObject.h"
#include "NueAna/NuePIDHeader.h"
#include "Record/RecRecordImp.h"
using namespace std;
#include "Record/RecRecordImp.cxx"

ClassImp(NuePIDHeader)
template class RecRecordImp<NuePIDHeader>;

NuePIDHeader::NuePIDHeader():
   RecHeader(),
   fSnarl(0),
   fRun(0),
   fSubRun(0),
   fEvtNo(0),
   fEvents(0),
   fDecider(NuePIDHeader::kUnknown)
{}

NuePIDHeader::NuePIDHeader(const VldContext &vld):
   RecHeader(vld),
   fSnarl(0),
   fRun(0),
   fSubRun(0),
   fEvtNo(0),
   fEvents(0),
   fDecider(NuePIDHeader::kUnknown)
{}

NuePIDHeader::~NuePIDHeader()
{}

void NuePIDHeader::Print(Option_t * /*option*/) const
{
   NuePIDHeader::Print(std::cout);
   return;
}

std::ostream &NuePIDHeader::Print(std::ostream &os) const
{
   os<<"Run: "<<fRun<<" SubRun "<<fSubRun<<" Snarl "<<fSnarl<<std::endl;
   os<<"This record corresponds to event "<<fEvtNo<<" out of "<<fEvents<<std::endl;
   os<<"The analysis that made this PID decision: "<<AsString(fDecider)<<endl;
   return os;
}
