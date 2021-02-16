#include <iostream>
#include "TObject.h"
#include "NueAna/NueHeader.h"
#include "Record/RecRecordImp.h"
using namespace std;
#include "Record/RecRecordImp.cxx"

ClassImp(NueHeader)
template class RecRecordImp<NueHeader>;

NueHeader::NueHeader():
   RecHeader(),
   fSnarl(0),
   fRun(0),
   fSubRun(0),
   fEvtNo(0),
   fEvents(0),
   fTrackLength(0),
   foundSR(false),
   foundMC(false),
   foundTH(false),
   foundMR(false)
{}

NueHeader::NueHeader(const VldContext &vld):
   RecHeader(vld),
   fSnarl(0),
   fRun(0),
   fSubRun(0),
   fEvtNo(0),
   fEvents(0),
   fTrackLength(0),
   foundSR(false),
   foundMC(false),
   foundTH(false),
   foundMR(false)
{}

NueHeader::~NueHeader()
{}

void NueHeader::Print(Option_t * /*option*/) const
{
   NueHeader::Print(std::cout);
   return;
}

std::ostream &NueHeader::Print(std::ostream &os) const
{
   os<<"Run: "<<fRun<<" SubRun "<<fSubRun<<" Snarl "<<fSnarl<<std::endl;
   os<<"This record corresponds to event "<<fEvtNo<<" out of "<<fEvents<<std::endl;
   os<<"Track has "<<fTrackLength<<" planes"<<std::endl;
   if(foundSR){
      os<<"Found SR"<<std::endl;
   }
   if(foundMC){
      os<<"Found MC"<<std::endl;
   }
   if(foundTH){
      os<<"Found TH"<<std::endl;
   }
   if(foundMR){
      os<<"Found MR"<<std::endl;
   }
   return os;
}
