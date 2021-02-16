#ifndef EVENTQUAL
#define EVENTQUAL

#include "TObject.h"

class NueRecord;

class EventQual : public TObject{ 

public:
   EventQual();
   EventQual(const EventQual *eq);
   virtual ~EventQual();

   void Reset();

   //Event based cuts   
   Float_t liTime;         
   Int_t passLI;             // LISeive
   Int_t passCosmicCut;      // Steve's cosmic cuts
   Int_t rcBoundary;
   Int_t isLargestEvent;     //Only meant for the far detector
   Int_t dmxstatus;   

   //snarl or time based cuts
   Int_t passFarDetQuality;  // derived from the IsGoodData Function
   Int_t passNearDetQuality;  // derived from the IsGoodData Function

   Int_t triggerSource;
   Double_t triggerTime;
   Int_t spillType;
   Int_t coilStatus;
   Int_t coilQuality;
   Int_t coilDirection;

   Int_t edgeActivityStrips;
   Float_t edgeActivityPH;
   Int_t oppEdgeStrips;
   Float_t oppEdgePH;
   
   Double_t minTimeSeparation;
   Double_t closeTimeDeltaZ;
   

private:

   ClassDef(EventQual,2)
   
};

#endif

