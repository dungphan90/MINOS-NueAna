#include <algorithm>

#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "MessageService/MsgService.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/EventQualAna.h"

#include "DatabaseInterface/DbiResultPtr.h"
#include "DcsUser/Dcs_Mag_Near.h"
#include "Mad/fddataquality.h"

#include "DcsUser/CoilTools.h"
#include "Conventions/Munits.h"
#include "Conventions/Detector.h"
#include "DataUtil/PlaneOutline.h"
#include "DataUtil/DataQualDB.h"
                                                                                
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NtupleUtils/LISieve.h"
#include "NueAna/NueStandard.h"

CVSID("$Id: EventQualAna.cxx,v 1.7 2009/07/30 16:30:17 pawloski Exp $");

EventQualAna::EventQualAna(NueRecord &nr):
    nr(nr)
{}

EventQualAna::~EventQualAna()
{
  
}

void EventQualAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
    NtpStRecord *st = dynamic_cast<NtpStRecord*>(srobj);

//    cout<<"Starting EventQualAna - "<<evtn<<endl;
    if(st == 0) return;
 
    NtpSREventSummary *eventSummary = &(st->evthdr);
    NtpSRDmxStatus *dmxSummary = &(st->dmxstatus);
    NtpSRDetStatus *detStatus = &(st->detstatus);
    const RecCandHeader *header = &(st->GetHeader());
                                                                                                             
    int det = header->GetVldContext().GetDetector();
    SimFlag::SimFlag_t s = header->GetVldContext().GetSimFlag();
    VldContext vld = header->GetVldContext();

    
    if(s == SimFlag::kMC && evtn == -10) return; //if "empty" and mc don't need this all

    nr.eventq.triggerSource = header->GetTrigSrc();
    nr.eventq.spillType   = header->GetRemoteSpillType();
//    nr.eventq.triggerTime = eventSummary->trigtime;
    nr.eventq.coilStatus  = detStatus->coilstatus;

    nr.eventq.liTime = eventSummary->litime;
    int isLI = 0;
    if(nr.eventq.liTime != -1) isLI = 1;
    else{ isLI = (int) LISieve::IsLI((*st)); }
    
    nr.eventq.passLI = 1 - isLI;  // 1 when its a good event

//    cout<<"LITime - "<<isLI<<endl;                                                                                                             
    nr.eventq.coilQuality = 0;
    nr.eventq.coilDirection = 0;

    if(s != SimFlag::kMC){
       nr.eventq.coilQuality = CoilTools::IsOK(vld);
       if(CoilTools::IsReverse(vld))  nr.eventq.coilDirection = -1;
       else                           nr.eventq.coilDirection = 1;
    }

//    cout<<nr.eventq.coilQuality<<"  "<<nr.eventq.coilDirection<<endl;
                                                                                                             
    int dmxStatus = 0;
    dmxStatus += (dmxSummary->nonphysicalfail);
    dmxStatus += ((dmxSummary->validplanesfail<<1));
    dmxStatus += ((dmxSummary->vertexplanefail<<2));
    nr.eventq.dmxstatus = dmxStatus;
    
    nr.eventq.rcBoundary = 0;
    if(det == Detector::kFar)
    {
      nr.eventq.rcBoundary = FDRCBoundary(eventSummary);
    }

    nr.eventq.passFarDetQuality = 1;    
    if(s != SimFlag::kMC) nr.eventq.passFarDetQuality = DataUtil::IsGoodData(st);
    nr.eventq.passNearDetQuality = nr.eventq.passFarDetQuality;
    //Greg July 29, 2009: I added the variable passNearDetQuality 
    //passNearDetQuality is the same as passFarDetQuality 
    //I would have preferred to have renamed passFarDetQuality->passDetQuality 
    //but I worried that it might have caused problems with older nueana files   
    //and I didn't want to store ND quailty in something that is called FD
    //so the two identical variables, maybe they will diverge in the future 

    if(evtn == -10) return;  //exit if this is an "empty" snarl

    nr.eventq.passCosmicCut = NueStandard::PassesCosmicCutFunction(&nr);

    int lgst_evt_idx=-1;
    float big_ph=0.0;
    for(int i=0;i<eventSummary->nevent;i++){ 
      NtpSREvent *evtTemp = SntpHelpers::GetEvent(i,srobj);
      if(evtTemp == 0) continue; //no event found
      if(evtTemp->ph.mip > big_ph){ 
	big_ph=evtTemp->ph.mip; 
	lgst_evt_idx=i;
      }
    }

    nr.eventq.isLargestEvent = 0;    
    if(lgst_evt_idx == evtn)  nr.eventq.isLargestEvent = 1;

    NtpSREvent *event = SntpHelpers::GetEvent(evtn,st);

    //Loop over all the snarls
    Int_t edgeActivityStrips = 0;
    Float_t edgeActivityPH = 0;
    Int_t oppEdgeStrips = 0;
    Float_t oppEdgePH = 0;
                                                                                                             
    Double_t thisEvtTime = 99999.;
                                                                                
    Int_t thisVtxPlane = 0;
    Float_t thisVtxZ = 0.0;                                                                                
    thisVtxPlane = event->vtx.plane;
    thisVtxZ = event->vtx.z;
                                                                                
    vector<Double_t> stpTimes; // container needed to calculate median
                                                                                
    // loop over strips in current event
    for(int i = 0; i < event->nstrip; i++){
       Int_t index = SntpHelpers::GetStripIndex(i,event);
       NtpSRStrip *strip = SntpHelpers::GetStrip(index,srobj);
       if(!strip)  continue;
     
       if (strip->plane >= thisVtxPlane 
             && strip->plane  < (thisVtxPlane+5)) { 
          stpTimes.push_back(strip->time1);
       }
    }
    // calculate strip time median for this event
    sort(stpTimes.begin(), stpTimes.end());
 
    if(stpTimes.size() > 0){

    if (stpTimes.size()%2)
      thisEvtTime = *(stpTimes.begin()+stpTimes.size()/2);
    else
      thisEvtTime = ( *(stpTimes.begin()+stpTimes.size()/2)
                    + *(stpTimes.begin()+stpTimes.size()/2-1) )/2.;
                                                                              
    Double_t timeDiff = 99999.;
    int nevent = eventSummary->nevent;
    for(int i = 0; i < nevent; ++i){
      NtpSREvent *evtTemp = SntpHelpers::GetEvent(i,srobj);
      if(evtTemp == 0) continue; //no event found
      // don't want to compare time to itself
      if (i == evtn) continue;
                                                                                
      int vtxPlane = evtTemp->vtx.plane;
      float vtxZ = evtTemp->vtx.z;
   
      Double_t evtTime = 99999.;
      stpTimes.clear();
                                                                                
      // loop over strips in other event
      for(Int_t j = 0; j < evtTemp->nstrip; ++j){
         Int_t index = SntpHelpers::GetStripIndex(j,evtTemp);
         NtpSRStrip *strip = SntpHelpers::GetStrip(index,st);
         if(!strip)  continue;
         if (strip->plane >= vtxPlane && strip->plane < (vtxPlane+5)) {
           stpTimes.push_back(strip->time1);
         }
      }
      // calculate strip time median for this event
      sort(stpTimes.begin(), stpTimes.end());
      if (stpTimes.size()%2)
        evtTime = *(stpTimes.begin()+stpTimes.size()/2);
      else
        evtTime= ( *(stpTimes.begin()+stpTimes.size()/2)
               + *(stpTimes.begin()+stpTimes.size()/2-1) )/2.;
                                                                                
      Double_t deltaT = evtTime - thisEvtTime;
      if( TMath::Abs(deltaT) < TMath::Abs(timeDiff)) {
        timeDiff = deltaT;
        nr.eventq.closeTimeDeltaZ = thisVtxZ - vtxZ;
      }
  }
  //end loop over events to find min time separation between events
  nr.eventq.minTimeSeparation = timeDiff;
                                                                                
  // calculate edge activity variables
  Double_t activityTimeStart = thisEvtTime - 4e-8; // +/- 40ns window
  Double_t activityTimeStop  = thisEvtTime + 4e-8;

  // loop over strips in snarl
  for(unsigned int i = 0; i < eventSummary->nstrip; ++i){
     NtpSRStrip *strip = SntpHelpers::GetStrip(i,st);
     if(!strip)  continue;

     // consider only strips in time window
     if (strip->time1 > activityTimeStart
        && strip->time1 < activityTimeStop) {
                                                                                
      // look at U planes in calorimeter
      if (strip->plane%2==1 && (strip->tpos<-0.24)
          && strip->plane<121) {
        edgeActivityStrips++;
        edgeActivityPH += strip->ph1.sigcor;
      }
      // look at opposite edge in U (use 3 strips of fully instr.)
      if (strip->plane%2==1 && (strip->tpos>2.27)
          && strip->plane<121) {
        oppEdgeStrips++;
        oppEdgePH += strip->ph1.sigcor;
      }
      // look at V planes in calorimeter
      if (strip->plane%2==0 && (strip->tpos>0.24)
          && strip->plane<121) {
        edgeActivityStrips++;
        edgeActivityPH += strip->ph1.sigcor;
      }
      // look at opposite edge in V (use 3 strips of fully instr.)
      if (strip->plane%2==0 && (strip->tpos<-2.27)
          && strip->plane<121) {
        oppEdgeStrips++;
        oppEdgePH += strip->ph1.sigcor;
      }
                                                                                
    } // time window
  }
  }
  nr.eventq.edgeActivityPH = edgeActivityPH;
  nr.eventq.edgeActivityStrips = edgeActivityStrips;
  nr.eventq.oppEdgePH = oppEdgePH;
  nr.eventq.oppEdgeStrips = oppEdgeStrips;
                                                                                
  //....................................................................



}
  
  
Int_t EventQualAna::FDRCBoundary(NtpSREventSummary *eventSummary){
  Int_t litag=0;
  Int_t numshower=eventSummary->nshower;
  Float_t allph=eventSummary->ph.raw;
  Int_t plbeg=eventSummary->plane.beg;
  Int_t plend=eventSummary->plane.end;

  if (numshower) {
    if (allph>1e6) litag+=10;

    if (((plbeg==1 || plbeg==2) && (plend==63 || plend==64)) ||
        ((plbeg==65 || plbeg==66) && (plend==127 || plend==128)) ||
        ((plbeg==129 || plbeg==130) && (plend==191 || plend==192)) ||
        ((plbeg==193 || plbeg==194) && (plend==247 || plend==248))) litag++;
    if (((plbeg==250 || plbeg==251) && (plend==312 || plend==313)) ||
        ((plbeg==314 || plbeg==315) && (plend==376 || plend==377)) ||
        ((plbeg==378 || plbeg==379) && (plend==440 || plend==441)) ||
        ((plbeg==442 || plbeg==443) && (plend==484 || plend==485))) litag++;
  }
  return litag;
}


