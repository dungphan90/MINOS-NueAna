#include "StandardNtuple/NtpStRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "MessageService/MsgService.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/ANtpEventInfoAna.h"
#include "AnalysisNtuples/Module/ANtpRecoNtpManipulator.h"

#include "DatabaseInterface/DbiResultPtr.h"
#include "DcsUser/Dcs_Mag_Near.h"
#include "Mad/fddataquality.h"

#include "DcsUser/CoilTools.h"
#include "Conventions/Munits.h"
#include "Conventions/Detector.h"
#include "DataUtil/PlaneOutline.h"
#include "DataUtil/DataQualDB.h"
                                                                                
#include "NueAna/NueAnaTools/NueConvention.h"
#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"

CVSID("$Id: ANtpEventInfoAna.cxx,v 1.42 2008/11/19 18:22:51 rhatcher Exp $");

ANtpEventInfoAna::ANtpEventInfoAna(ANtpEventInfoNue &anei):
    fANtpEventInfo(anei),
    fDetectorType(Detector::kUnknown)
{}

ANtpEventInfoAna::~ANtpEventInfoAna()
{
  
}

//void ANtpEventInfoAna::Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord * /*mcobj*/, NtpTHRecord * /*thobj*/)
void ANtpEventInfoAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{

    fANtpEventInfo.index = evtn;

//    const RecCandHeader *ntpHeader = &(srobj->GetHeader());

    ANtpRecoNtpManipulator *ntpManipulator = 0;
    NtpStRecord *st = dynamic_cast<NtpStRecord *>(srobj);
    NtpSRRecord *sr = 0;   
    NtpSREventSummary *eventSummary;
    NtpSRDmxStatus *dmxSummary;
    NtpSRDetStatus *detStatus;
    const RecCandHeader *header;
                          
    if(st != 0){
      ntpManipulator = new ANtpRecoNtpManipulator(st);
      header = &(st->GetHeader());
      eventSummary = &(st->evthdr);
      detStatus = &(st->detstatus);
      dmxSummary = &(st->dmxstatus);
    }
    else{
      sr = dynamic_cast<NtpSRRecord *>(srobj);
      ntpManipulator = new ANtpRecoNtpManipulator(sr);
      header = &(sr->GetHeader());
      eventSummary = &(sr->evthdr);
      detStatus = &(sr->detstatus);
      dmxSummary = &(sr->dmxstatus);
    }

    fDetectorType = header->GetVldContext().GetDetector();
    SimFlag::SimFlag_t s = header->GetVldContext().GetSimFlag();
    if(s == SimFlag::kMC && evtn == -10) return;
                                                                               

    MSG("ANtpEventInfoAna",Msg::kDebug)<<"Filled event information specific to nue" << endl; 

    fANtpEventInfo.triggerSource = header->GetTrigSrc();
    fANtpEventInfo.spillType = header->GetRemoteSpillType();
    fANtpEventInfo.liTime = eventSummary->litime;    

//    fANtpEventInfo.triggerTime = eventSummary->trigtime;
    fANtpEventInfo.coilStatus  = detStatus->coilstatus;
                                                                                
    fANtpEventInfo.coilCurrent = 0;
    fANtpEventInfo.coilQuality = 0;
    VldContext vld = header->GetVldContext();

    if(s != SimFlag::kMC){
       std::pair<float, float> temp = CoilTools::CoilCurrent(vld);

       fANtpEventInfo.coilCurrent = temp.first;
       if(vld.GetDetector() == Detector::kFar)
          fANtpEventInfo.coilCurrentSM2 = temp.second;
       fANtpEventInfo.coilQuality = CoilTools::IsOK(vld);
       if(CoilTools::IsReverse(vld))  fANtpEventInfo.coilDirection = -1; 
       else                           fANtpEventInfo.coilDirection = 1;
    }

    fANtpEventInfo.dmxStatus = 0;
    fANtpEventInfo.dmxStatus += (dmxSummary->nonphysicalfail);
    fANtpEventInfo.dmxStatus += ((dmxSummary->validplanesfail<<1));
    fANtpEventInfo.dmxStatus += ((dmxSummary->vertexplanefail<<2));


    fANtpEventInfo.rcBoundary = FDRCBoundary(eventSummary);
    fANtpEventInfo.daveFDDataQuality = DataUtil::IsGoodFDData(st);

    if(evtn == -10) return;


    NtpSREvent *event = 0;
    event = SntpHelpers::GetEvent(evtn,srobj);
    if(!event){
        MSG("ANtpEventInfoAna",Msg::kError)<<"Couldn't get event "
               <<evtn<<" from Snarl "<<srobj->GetHeader().GetSnarl()<<endl;
        return;
    }


    fInfoFiller = new ANtpInfoObjectFiller();
        MSG("ANtpEventInfoAna",Msg::kDebug)<<"Created manipulator and filler "
                                           << ntpManipulator << " "
                                           << fInfoFiller <<endl;

//  Fill information from the base ANtpEventINfo
//    TClonesArray* temp = ntpManipulator->GetStripArray();
//    MSG("ANtpEventInfoAna",Msg::kDebug)<<"Filled array " << temp->GetEntries()
//                                     <<endl;
    fInfoFiller->SetStripArray(ntpManipulator->GetStripArray());
    MSG("ANtpEventInfoAna",Msg::kDebug)<<"SetStripArray " << fInfoFiller;
    fInfoFiller->FillEventInformation(ntpManipulator, event, &fANtpEventInfo);
    MSG("ANtpEventInfoAna",Msg::kDebug)<<"Filled event information "
                                       << fInfoFiller << endl;

    if(ReleaseType::IsCedar(release))
    {
      Float_t vtxx, vtxy, vtxz;

      NtpVtxFinder vtxf;
      vtxf.SetTargetEvent(evtn, st);
      if(vtxf.FindVertex() > 0){
        fANtpEventInfo.vtxX = vtxx = vtxf.VtxX();
        fANtpEventInfo.vtxY = vtxy = vtxf.VtxY();
        fANtpEventInfo.vtxZ = vtxz = vtxf.VtxZ();
        fANtpEventInfo.vertexTime = vtxf.VtxT();

        Detector::Detector_t det = ( Detector::Detector_t) fDetectorType;
        fANtpEventInfo.vtxMetersToBeam =
            fInfoFiller->MetersToBeam(det, vtxx, vtxy, vtxz);
        fANtpEventInfo.vtxMetersToCoil =
             fInfoFiller->MetersToCoil(det, vtxx, vtxy);
        fANtpEventInfo.vtxMetersToCloseEdge =
             fInfoFiller->MetersToCloseEdge(det, vtxx,vtxy, 0);
     }
   }

    FillNueEventInformation(event, &fANtpEventInfo);


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

    fANtpEventInfo.largestEvent = 0;    
    if(lgst_evt_idx == evtn)
       fANtpEventInfo.largestEvent = 1;

    FillStripVariables(event, srobj);
     
  if(fInfoFiller){
    delete fInfoFiller;
    fInfoFiller=0;
  }
  if(ntpManipulator){
    delete ntpManipulator;
    ntpManipulator=0;
  }

}
    
//----------------------------------------------------------------------
void ANtpEventInfoAna::FillNueEventInformation(NtpSREvent *ntpEvent, ANtpEventInfoNue *eventInfoNue)
{

    eventInfoNue->timeLength=ntpEvent->end.t-ntpEvent->vtx.t; // event time lenght
    eventInfoNue->phMeu =  ntpEvent->ph.mip;
    eventInfoNue->phMip =  ntpEvent->ph.mip;
    eventInfoNue->phNueGeV = ntpEvent->ph.mip/MeuPerGeV;
   
    return;
}

Float_t ANtpEventInfoAna::GetNDCoilCurrent(const VldContext& vc)
{
   const Dcs_Mag_Near* magnear =
         CoilTools::Instance().GetMagNear(vc);  // NearDet only
   if ( magnear ) return magnear->GetCurrent();
   else           return 0;
}

Int_t ANtpEventInfoAna::FDRCBoundary(NtpSREventSummary *eventSummary){
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

void ANtpEventInfoAna::FillStripVariables(NtpSREvent *ntpEvent, RecRecordImp<RecCandHeader> *srobj)
{
 
  Int_t  HotChannel = 0;
  //cut and pasted, with minor changes from MadDpAnalysis.cxx

  double trgtime= srobj->GetHeader().GetVldContext().GetTimeStamp().GetNanoSec()/1.e9;                                                                                
  double timemax=0.;
  double timemin=1.e10;

  Float_t sigfull,sigpart;
  sigfull=sigpart=0;

  int triggerPlanes = 4; //4/5 plane trigger
  int nDetPlanes = 0;
  if(fDetectorType == Detector::kNear)   nDetPlanes = 122;
  if(fDetectorType == Detector::kFar)    nDetPlanes = 485;
                                                                                     
  vector<int> dataInPlane;
  dataInPlane.assign(nDetPlanes, 0);

  for(int i=0;i<ntpEvent->nstrip;i++){
    Int_t index = SntpHelpers::GetStripIndex(i,ntpEvent);
    NtpSRStrip *ntpStrip = SntpHelpers::GetStrip(index,srobj);
    if(!ntpStrip){
      MSG("ANtpEventInfoAna",Msg::kError)<<"Couldn't get strip "
	      <<index<<" from "<<" snarl "
	      <<srobj->GetHeader().GetSnarl()
	      <<" something has gone horribly wrong, I'll just go"
	      <<" on to the next strip"<<endl;
      continue;
    }
    
    if(fDetectorType==Detector::kNear)
       if(ntpStrip->ph0.raw+ntpStrip->ph1.raw>60000) HotChannel = 1;

    //Strip Timing Information
    Double_t striptime = 0;
    if(fDetectorType==Detector::kNear){
      striptime=ntpStrip->time1-trgtime;
    }
    if(fDetectorType==Detector::kFar){
      Double_t striptime1=ntpStrip->time1-trgtime;
      Double_t striptime0=ntpStrip->time0-trgtime;
      if(striptime1>0 && striptime0<0)  striptime=striptime1;
      if(striptime0>0 && striptime1<0)  striptime=striptime0;
      if(striptime0>0 && striptime1>0)  striptime=(striptime0+striptime1)/2.; 
    }

    if(striptime<=timemin) timemin=striptime;
    if(striptime>=timemax) timemax=striptime;
    //End of Strip timing information

    // compute the amount of signal in the partially instrumented region
    // and the amount in the fully instrumented region
    //
    // this region is defined as:
    // v planes: (strip<=4 || strip>=67)
    // partial u: (strip==0 || strip=63)
    // full u: (strip<=26 || strip>=88) 

    if(fDetectorType == Detector::kNear){
      Int_t pr = NueConvention::InPartialRegion(ntpStrip->plane, ntpStrip->strip);
      if(!evtstp0mip){
       MSG("ANtpEventInfoAna",Msg::kError)<<"No mip strip information"<<endl;
       continue;
      }

      float charge = evtstp0mip[index] + evtstp1mip[index];
      if(pr==1)      {  sigpart += charge; }
      else if(pr==-1){  sigfull += charge; }
    }
                                                                                     
    if(ntpStrip->plane <  nDetPlanes) dataInPlane[ntpStrip->plane] = 1;
  }
  
  bool triggerVerified = false;
  Int_t group = 0;
                                                                       
  for(Int_t loop = 0; loop < nDetPlanes && !triggerVerified; loop++) {
    if(dataInPlane[loop] > 0) {
       group = 1;
       for(Int_t l = 1; l < triggerPlanes+1; l++) {
         if(dataInPlane[loop+l] > 0) group++;
       }
       if(group >= triggerPlanes) triggerVerified = true;
     }
   }
                                
  dataInPlane.clear();                                                     

  if(fDetectorType == Detector::kFar){
    sigfull = ntpEvent->ph.mip;
    sigpart = 0;
  }
  
  fANtpEventInfo.eventSignalFull = sigfull;
  fANtpEventInfo.eventSignalPartial = sigpart;

  fANtpEventInfo.triggerPass = (int) triggerVerified;
    
  fANtpEventInfo.eventTimeMax=timemax;
  fANtpEventInfo.eventTimeMin=timemin;
  fANtpEventInfo.hotch = HotChannel;


//    evttimeln=timemax-timemin;

}

