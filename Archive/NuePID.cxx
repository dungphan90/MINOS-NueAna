#include "NueAna/NuePID.h"

ClassImp(NuePID)

NuePID::NuePID():
   RecRecordImp<NuePIDHeader>(),
   IsNue(0),
   likelihood(0.)
{}

NuePID::NuePID(const NuePIDHeader &pidhead):
   RecRecordImp<NuePIDHeader>(pidhead),
   IsNue(0),
   likelihood(0.)
{}

NuePID::~NuePID()
{}
