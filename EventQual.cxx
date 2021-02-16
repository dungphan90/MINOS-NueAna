#include <iostream>
#include "TClonesArray.h"
#include "NueAna/EventQual.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(EventQual)

EventQual::EventQual() 
{
  Reset(); 
}

EventQual::EventQual(const EventQual *eq):
   liTime(eq->liTime),
   passLI(eq->passLI),
   passCosmicCut(eq->passCosmicCut),
   rcBoundary(eq->rcBoundary),
   isLargestEvent(eq->isLargestEvent),
   dmxstatus(eq->dmxstatus),                                                                                
   passFarDetQuality(eq->passFarDetQuality),
   passNearDetQuality(eq->passNearDetQuality),
   triggerSource(eq->triggerSource),
   triggerTime(eq->triggerTime),
   spillType(eq->spillType),                                                                                
   coilStatus(eq->coilStatus),
   coilQuality(eq->coilQuality),
   coilDirection(eq->coilDirection),                                                              
   edgeActivityStrips(eq->edgeActivityStrips),
   edgeActivityPH(eq->edgeActivityPH),
   oppEdgeStrips(eq->oppEdgeStrips),
   oppEdgePH(eq->oppEdgePH),   
   minTimeSeparation(eq->minTimeSeparation),
   closeTimeDeltaZ(eq->closeTimeDeltaZ)
{
}

EventQual::~EventQual()
{
  //  std::cout<<"in EventQual destructor"<<std::endl;
}


void EventQual::Reset()
{
  liTime = ANtpDefVal::kFloat;

  passLI = ANtpDefVal::kInt;
  passCosmicCut = ANtpDefVal::kInt;
  rcBoundary = ANtpDefVal::kInt;
  isLargestEvent = ANtpDefVal::kInt;
  dmxstatus = ANtpDefVal::kInt;
  passFarDetQuality = ANtpDefVal::kInt;
  passNearDetQuality = ANtpDefVal::kInt;
  triggerSource = ANtpDefVal::kInt;
  triggerTime = ANtpDefVal::kDouble;
  spillType = ANtpDefVal::kInt;
  coilStatus = ANtpDefVal::kInt;
  coilQuality = ANtpDefVal::kInt;
  coilDirection = ANtpDefVal::kInt;

  edgeActivityStrips = ANtpDefVal::kInt;
  edgeActivityPH = ANtpDefVal::kFloat;
  oppEdgeStrips = ANtpDefVal::kInt;
  oppEdgePH = ANtpDefVal::kFloat;                                                                     
  minTimeSeparation = ANtpDefVal::kDouble;
  closeTimeDeltaZ = ANtpDefVal::kDouble;
}


