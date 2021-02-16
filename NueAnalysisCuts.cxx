#include "NueAna/NueAnalysisCuts.h"
#include "Calibrator/CalMIPCalibration.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRTrack.h"
#include "CandNtupleSR/NtpSRShower.h"
#include "StandardNtuple/NtpStRecord.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "MessageService/MsgService.h"
#include <fstream>
#include "TList.h"
#include "TClass.h"
#include "TDataType.h"
#include "TDataMember.h"
#include "TRealData.h"                  
#include "DataUtil/infid.h"

                                                              
ClassImp(NueAnalysisCuts)
CVSID("$Id: NueAnalysisCuts.cxx,v 1.18 2010/07/30 01:42:31 rhatcher Exp $");

#include "DatabaseInterface/DbiResultPtr.tpl"

NueAnalysisCuts::NueAnalysisCuts()
{
   kMEUPerGeV =25.66;
   kSigCorrPerMEU = 1.0;
   Reset();
}

void NueAnalysisCuts::Reset()
{
   LowCuts.Reset();
   HighCuts.Reset();
                        
   fDetectorType = Detector::kUnknown;
   fSimFlag = SimFlag::kUnknown;
                                                        
   fCheckHotChannel = 0;
   fFiducialCut = 0;
   fFidVtxChoice = FiducialAlg::kUseEvtVtx;
   fFidVolumeChoice = FiducialAlg::kExtendedVolume;
   fContainmentCut = 0;
   fBeamCut = 0;
   fTargetCurrent = -1000;
   fPhProngCut = -1;
   fFileTrim = 0;
   fEventNum = 0;
}


void NueAnalysisCuts::SetInfoObject(NueRecord* nr)
{
    nrInfo = nr;

    stInfo = NULL;
    srInfo = NULL;
    eventInfo = NULL;
    trackInfo = NULL;
    mcInfo = NULL;
    showerInfo = NULL;

    fDetectorType = nr->GetHeader().GetVldContext().GetDetector();
    fSimFlag = nrInfo->GetHeader().GetVldContext().GetSimFlag();
}
                                                                                
void NueAnalysisCuts::SetInfoObject(int evtn, NtpStRecord* st)
{
    if(!st) return;

    stInfo = st;

    VldContext vc;
    vc=st->GetHeader().GetVldContext();
 
    int run = st->GetHeader().GetRun();
    ntpManipulator.SetRecord(st);
    SetNtpInfoObjects(evtn,run, vc);
    
    nrInfo = NULL;
}

void NueAnalysisCuts::SetInfoObject(int evtn, NtpSRRecord* st, NtpMCRecord * mcobj, NtpTHRecord * thobj)
{
    if(!st) return;                                                                                                       
    srInfo = st;                                                                   
    VldContext vc;
    vc=st->GetHeader().GetVldContext();

    int run = st->GetHeader().GetRun();
    ntpManipulator.SetRecord(st, mcobj, thobj);
    SetNtpInfoObjects(evtn,run,vc);                                                                                       
    nrInfo = NULL;
}

void NueAnalysisCuts::SetNtpInfoObjects(int evtn, int run, VldContext &vc)
{
    fDetectorType = vc.GetDetector();
    fSimFlag = vc.GetSimFlag();

    fEventNum = evtn;
    ntpManipulator.SetPrimaryShowerCriteria(0,1); // nplanes, ph
    ntpManipulator.SetPrimaryTrackCriteria(0,1,0); //nplanes, length, ph

    ANtpEventManipulator * ntpEventManipulator =
                              ntpManipulator.GetEventManipulator();
    ntpEventManipulator->SetEventInSnarl(evtn);
 
    eventInfo=ntpEventManipulator->GetEvent();
    showerInfo = ntpEventManipulator->GetPrimaryShower();
    trackInfo = SntpHelpers::GetPrimaryTrack(evtn, stInfo);
    thEventInfo = ntpManipulator.GetMCManipulator()->GetNtpTHEvent(evtn);

    if(thEventInfo){
      mcInfo = ntpManipulator.GetMCManipulator()->GetNtpMCTruth(thEventInfo->neumc);
//      mcstdhep = ntpManipulator.GetMCManipulator()->GetNtpMCStdHep(thevent->neustdhep);
    }
               
    static int currentRun = 0;
    if(run != currentRun)
    {
        DbiResultPtr<CalMIPCalibration> dbp(vc);
        if(dbp.GetNumRows()>0){
           const CalMIPCalibration *m = dbp.GetRow(0);
           float mip=m->GetMIP(1.);
           if(mip>0){
             kSigCorrPerMEU=1./mip;
           }
        }
        currentRun = run;
    }
}

// End of Filling Information Variables
const Registry& NueAnalysisCuts::DefaultConfig() const
{
  static Registry r;

  r.UnLockValues();

  r.Set("HiPlaneTrackCut",-1);
  r.Set("LoPlaneEventCut",-1);
  r.Set("HiPlaneEventCut",-1);
  r.Set("HiTrackLikeCut",-1);
  r.Set("CheckHotChannel", 0);
  r.Set("FiducialCut", 0); 
  r.Set("FiducialVtx", 1); 
  r.Set("FiducialVolumeAlg", 2);
  r.Set("HiEnergyCut",-1);
  r.Set("LoEnergyCut",-1);
  r.Set("HiShowerEnergyCut",-1);
  r.Set("LoShowerEnergyCut",-1);

  // new contPlaneCount - Minerba
  //  r.Set("LoContNPlaneCut",-1);

  r.Set("CutOnClasses", 0);
  r.Set("BeamQualityCut", 0);
  r.Set("LoPhNStripCut",-1);
  r.Set("LoPhNPlaneCut",-1);
  r.Set("HiEnergyCut",-1);
  r.Set("LoEnergyCut",-1);
  r.Set("HiEnergyShowerCut",-1);
  r.Set("LoEnergyShowerCut",-1);

  // ?? new contPlaneCount - Minerba
  //  r.Set("LoContPhPlaneCut",-1.);

  r.Set("LoCurrentCut", 0.1);
  r.Set("LoHorBeamWidth", 0.0);
  r.Set("HiHorBeamWidth", 2.9);
  r.Set("LoVertBeamWidth", 0.0);
  r.Set("HiVertBeamWidth", 2.9);
  r.Set("LoNuTarZ", -1);
  r.Set("HiNuTarZ", 1000);
  r.Set("Oscillate", 0);
  r.Set("OutputFile", "HistManInfo.root");

  r.Set("FileCut", 0); 
  r.Set("TrimFile", "Default.txt");

//  r.Set("SetBackground", 0);
//  r.Set("AddBackground", 0);
  r.Set("AddSignal", 2);
//  r.Set("SetSignal", 2);
//  r.Set("PhProngCut",0);
//  r.Set("HiSigEmFracCut", -1);
//  r.Set("LoSigEmFracCut", -1); 
//  r.Set("HiBgEmFracCut", -1);
//  r.Set("LoBgEmFracCut", -1);
  r.Set("TargetHornCurrent", -185); 
  r.Set("CutonResCode", 0);


  r.LockValues();

  return r;                                                                                
}

void NueAnalysisCuts::Config(const Registry& r)
{
  int imps;
  if(r.Get("HiPlaneTrackCut",imps)) { HighCuts.srtrack.planes=imps;}
  if(r.Get("HiTrackLikeCut",imps))  { HighCuts.srtrack.trklikePlanes=imps;}
  if(r.Get("LoPlaneEventCut",imps)) { LowCuts.srevent.planes=imps;}
  if(r.Get("HiPlaneEventCut",imps)) { HighCuts.srevent.planes=imps;}
  if(r.Get("CheckHotChannel",imps)) { fCheckHotChannel=imps;}
  if(r.Get("FiducialCut", imps))       { fFiducialCut = imps;}
  if(r.Get("FiducialVtx", imps))       { fFidVtxChoice = imps;}
  if(r.Get("FiducialVolumeAlg", imps)) { fFidVolumeChoice = imps;}
  if(r.Get("ContainmentCut", imps)) {fContainmentCut = imps;}

  if(r.Get("CutOnClasses", imps)) { fCutClasses = imps; }
  if(r.Get("SetSignal", imps)) {   
       signalClasses.clear(); signalClasses.push_back(imps);
  }
  if(r.Get("AddSignal", imps)) { signalClasses.push_back(imps);  }
  if(r.Get("SetBackground", imps)) {
       bgClasses.clear(); bgClasses.push_back(imps);
  }
  if(r.Get("AddBackground", imps)) { bgClasses.push_back(imps); }

  if(r.Get("SetSignalResCode",imps)) {
     sigResonanceCodes.clear(); sigResonanceCodes.push_back(imps); }
  if(r.Get("SetBackgroundResCode",imps)) {
     bgResonanceCodes.clear(); bgResonanceCodes.push_back(imps); }
                                                                             
  if(r.Get("AddSignalResCode",imps)) {sigResonanceCodes.push_back(imps); }
  if(r.Get("AddBackgroundResCode",imps)) {bgResonanceCodes.push_back(imps); }
  if(r.Get("CutonResCode", imps)) {fCutResCode = imps;}

  if(r.Get("BeamQualityCut", imps)) {fBeamCut = imps;}
  if(r.Get("LoPhNStripCut",imps)) { LowCuts.shwfit.hiPhStripCount=imps;}
  if(r.Get("LoPhNPlaneCut",imps)) { LowCuts.shwfit.hiPhPlaneCountM2=imps;}

  // new contPlaneCount - Minerba
  //if(r.Get("LoContNPlaneCut",imps)) { LowCuts.shwfit.contPlaneCount=imps;}

  if(r.Get("FileCut", imps)) {fFileTrim = imps;}  
//  if(r.Get("LoNuTarZ", imps)) {kLoNuTarZ = imps;}
//  if(r.Get("HiNuTarZ", imps)) {kHiNuTarZ = imps;}
                   
  double fmps;
  if(r.Get("HiShowerEnergyCut",fmps)) { HighCuts.srshower.phNueGeV =fmps;}
  if(r.Get("LoShowerEnergyCut",fmps)) { LowCuts.srshower.phNueGeV =fmps;}
  if(r.Get("HiEnergyCut",fmps)) { HighCuts.srevent.phNueGeV =fmps;}
  if(r.Get("LoEnergyCut",fmps)) { LowCuts.srevent.phNueGeV  =fmps;}
  if(r.Get("PhProngCut",fmps))  { fPhProngCut=fmps; }
  if(r.Get("HiSigEmFracCut", fmps)) {
		      HighSigCuts.emShowerFraction = fmps;}
  if(r.Get("LoSigEmFracCut", fmps)) {
		      LowSigCuts.emShowerFraction = fmps;}
  if(r.Get("HiBgEmFracCut", fmps)) {
		      HighBgCuts.emShowerFraction = fmps;}
  if(r.Get("LoBgEmFracCut", fmps)) {
	    	      LowBgCuts.emShowerFraction = fmps; }
  
  if(r.Get("TargetHornCurrent", fmps)) {fTargetCurrent = fmps;}

  // ?? new contPlaneCount - Minerba
  //  if(r.Get("LoContPhPlaneCut",fmps)) { contPlaneCount=fmps;}

  const char* tmps;

  if(r.Get("SetBackground", tmps)) {
    if( std::string(tmps) == "All"){
       bgClasses.clear(); 
       for(int i = 0; i < 5; i++) bgClasses.push_back(i);
    }
  }

  if(r.Get("SetSignal", tmps)) {
    if( std::string(tmps) == "All"){
       signalClasses.clear();
       for(int i = 0; i < 5; i++) signalClasses.push_back(i);
    }
  }

  if(r.Get("TrimFile", tmps)) {fFileCutName = tmps;}  
  if(r.Get("OutputFile", tmps)) {kOutputFile = tmps;} 

}

// ******************************************************************
//      Simplifying Functions
// ******************************************************************

bool NueAnalysisCuts::IsValid()
{
   bool retval = false;
   if(nrInfo || stInfo || srInfo) retval = true; 
   return retval;
}

bool NueAnalysisCuts::IsMC()
{
   if(!IsValid()) return false;
   return (fSimFlag == SimFlag::kMC);
}


// ******************************************************************
//      Start of Actual Cuts
// ******************************************************************

bool NueAnalysisCuts::PassesAllCuts()
{
     bool goodEvent = true;
     goodEvent = goodEvent && PassesTrackPlaneCut();
     goodEvent = goodEvent && PassesTrackLikePlaneCut();
     goodEvent = goodEvent && PassesHighShowerEnergyCut();
     goodEvent = goodEvent && PassesLowShowerEnergyCut();
     goodEvent = goodEvent && PassesHighEnergyCut();
     goodEvent = goodEvent && PassesLowEnergyCut();
     goodEvent = goodEvent && PassesEventPlaneCut();
     goodEvent = goodEvent && PassesHotChannel();
     goodEvent = goodEvent && PassesFiducialVolume();
     goodEvent = goodEvent && PassesFullContainment();
     goodEvent = goodEvent && PassesPhProngCut();
     goodEvent = goodEvent && PassesCutOnClasses();
                                                                                
     goodEvent = goodEvent && PassesLowPhNPlaneCut();
     goodEvent = goodEvent && PassesLowPhNStripCut();
     goodEvent = goodEvent && PassesEMFraction();
     goodEvent = goodEvent && PassesResonanceCode();
     goodEvent = goodEvent && PassesFileCut();
                                                                                
     return goodEvent;
}

// ******************************************************************
//      Plane Cuts
// ******************************************************************

bool NueAnalysisCuts::PassesTrackPlaneCut()
{
   if(!IsValid()) return false;
   if(HighCuts.srtrack.planes < 0) return true;

   int planes = 0;
   if(nrInfo) planes = nrInfo->srtrack.planes;
   if(stInfo && trackInfo) planes = trackInfo->plane.n;

   bool goodEvent = !(planes > HighCuts.srtrack.planes);
   
   return goodEvent;
}

bool NueAnalysisCuts::PassesEventPlaneCut()
{
   if(!IsValid()) return false;
   
   int planes = 0;
   if(nrInfo) planes = nrInfo->srevent.planes; 
   if(eventInfo) planes = eventInfo->plane.n;
                                                                             
   bool goodEvent = true;
   if(HighCuts.srevent.planes > 0){
     goodEvent = !(planes > HighCuts.srevent.planes);
   }
   if(LowCuts.srevent.planes > 0){
     goodEvent = goodEvent && (planes >= LowCuts.srevent.planes);
   }

   return goodEvent;
}

bool NueAnalysisCuts::PassesTrackLikePlaneCut()
{
   if(!IsValid()) return false;
   if(HighCuts.srtrack.trklikePlanes < 0) return true;
                                                                                
   int planes = 0;
   if(nrInfo) planes = nrInfo->srtrack.trklikePlanes;
   if(stInfo && trackInfo) planes = trackInfo->plane.ntrklike;
                                                                                
   bool goodEvent = planes <= HighCuts.srtrack.trklikePlanes;
                                                                                
   return goodEvent;
}

// ******************************************************************
//      Energy Cuts
// ******************************************************************

bool NueAnalysisCuts::PassesHighShowerEnergyCut()
{
   if(!IsValid()) return false;
                                                                                
   bool goodEvent = true;
   if(HighCuts.srshower.phNueGeV < 0) return goodEvent;
                                                                                
   if(nrInfo){
      if(nrInfo->srshower.phNueGeV > HighCuts.srshower.phNueGeV)
        goodEvent = false;
   }else
   if(stInfo && showerInfo){
      float shwEnergyInNueGeV = showerInfo->ph.mip;
      shwEnergyInNueGeV = shwEnergyInNueGeV/kMEUPerGeV;

      if(shwEnergyInNueGeV > HighCuts.srshower.phNueGeV)
        goodEvent = false;
   }
                                                                                
   return goodEvent;
}

bool NueAnalysisCuts::PassesLowShowerEnergyCut()
{
   if(!IsValid()) return false;
   
   bool goodEvent = true;
   if(LowCuts.srshower.phNueGeV < 0) return goodEvent;
   
   if(nrInfo){
      int nShw = nrInfo->srevent.showers;                                                                                
      if(nShw == 0) goodEvent = false;
      if(nShw == 1 && nrInfo->srshower.phNueGeV < LowCuts.srshower.phNueGeV)
        goodEvent = false;
      if(nShw > 1 && nrInfo->srshower.phNueGeV < LowCuts.srshower.phNueGeV/2)
        goodEvent = false;

   }else
   if(stInfo && showerInfo){
      float EnergyInNueGeV = showerInfo->ph.mip;
      EnergyInNueGeV = EnergyInNueGeV/kMEUPerGeV;
      int nShw = eventInfo->nshower;

      if(nShw == 0) goodEvent = false;
      if(nShw == 1 && EnergyInNueGeV < LowCuts.srshower.phNueGeV)
        goodEvent = false;
      if(nShw > 1 && EnergyInNueGeV < LowCuts.srshower.phNueGeV/2)
        goodEvent = false;
   }
          
   return goodEvent;
}

bool NueAnalysisCuts::PassesHighEnergyCut()
{
   if(!IsValid()) return false;
   
   bool goodEvent = true;
   if(HighCuts.srevent.phNueGeV < 0) return goodEvent;
   
   if(nrInfo){
      if(nrInfo->srevent.phNueGeV > HighCuts.srevent.phNueGeV)
        goodEvent = false;
   }else
   if(stInfo && eventInfo){
      float EnergyInNueGeV = eventInfo->ph.mip;
      EnergyInNueGeV = EnergyInNueGeV/kMEUPerGeV;
                                        
      if(EnergyInNueGeV > HighCuts.srevent.phNueGeV)
        goodEvent = false;
   }
   
   return goodEvent;
}

bool NueAnalysisCuts::PassesLowEnergyCut()
{
   if(!IsValid()) return false;
   bool goodEvent = true;
   if(LowCuts.srevent.phNueGeV < 0) return goodEvent;
                                 
   if(nrInfo){
      if(nrInfo->srevent.phNueGeV < HighCuts.srevent.phNueGeV)
        goodEvent = false;
   }else
   if(stInfo && eventInfo){
      float EnergyInNueGeV = eventInfo->ph.mip;
      EnergyInNueGeV = EnergyInNueGeV/kMEUPerGeV;
                                           
      if(EnergyInNueGeV < HighCuts.srevent.phNueGeV)
        goodEvent = false;
   }
   
   return goodEvent;
}

// ******************************************************************
//      Fiducial Volume Cuts
// ******************************************************************

bool NueAnalysisCuts::PassesFullContainment()
{
    if(!IsValid()) return false;
                                                                                
    if(fContainmentCut == 0) return true;
                                                                                
    if(nrInfo == 0){
      MAXMSG("NueAnalysisCuts",Msg::kError,10)
        <<"FullContainment may only be called on ana_nue files"
        <<" the cut will not be applied"<<endl;
      return true;
    }
                                                                                
   bool goodEvent = false;
                                                                                
   if(fDetectorType == Detector::kFar)
      goodEvent = (nrInfo->anainfo.isFullyContained == 1);
                                                                                
   if(fDetectorType == Detector::kNear)
      goodEvent = (nrInfo->anainfo.isFullyContained == 1 ||
                    nrInfo->anainfo.isFullyContained == -2);
                                                                                
   return goodEvent;
}

bool NueAnalysisCuts::PassesFiducialVolume()
{
    if(!IsValid()) return false;
    if(fFiducialCut == 0) return true;

    bool goodEvent = false;
 
    float x = -9999; float y = -9999; float z = -9999;
    FillVertexPosition(x,y,z);

    if(x + y + z < -1000) return false;

    if(fDetectorType == Detector::kFar)
	goodEvent = IsInsideFarFiducial(x,y,z);
    if(fDetectorType == Detector::kNear)
        goodEvent = IsInsideNearFiducial(x,y,z);

   return goodEvent;
}

void NueAnalysisCuts::FillVertexPosition(float &x, float &y, float &z)
{
   if(nrInfo){
     if(fFidVtxChoice == FiducialAlg::kUseEvtVtx){
	x = nrInfo->srevent.vtxX;
        y = nrInfo->srevent.vtxY;
        z = nrInfo->srevent.vtxZ;
     }
     if(fFidVtxChoice == FiducialAlg::kUseTrkVtx){
        x = nrInfo->srtrack.vtxX;
        y = nrInfo->srtrack.vtxY;
        z = nrInfo->srtrack.vtxZ;
     }
   }

   if(eventInfo || trackInfo){
     if(eventInfo && fFidVtxChoice == FiducialAlg::kUseEvtVtx){
        x = eventInfo->vtx.x;
        y = eventInfo->vtx.y;
        z = eventInfo->vtx.z;
     }
     if(trackInfo && fFidVtxChoice == FiducialAlg::kUseTrkVtx){
        x = trackInfo->vtx.x;
        y = trackInfo->vtx.y;
        z = trackInfo->vtx.z;
     }
   }

   return;
}

bool NueAnalysisCuts::IsInsideNearFiducial(float x, float y, float z)
{
   if(fFidVolumeChoice == FiducialAlg::kDefaultVolume)
       return (infid(fDetectorType, fSimFlag, x,y,z) > 0);
   if(fFidVolumeChoice == FiducialAlg::kExtendedVolume){
       return (NueConvention::IsInsideNearFiducial_Nue_Extended(x,y,z) > 0);
   }
   return false;
}
                                                                                
bool NueAnalysisCuts::IsInsideFarFiducial(float x, float y, float z)
{
   if(fFidVolumeChoice == FiducialAlg::kDefaultVolume)
      return (infid(fDetectorType, fSimFlag, x,y,z) > 0);
   if(fFidVolumeChoice == FiducialAlg::kExtendedVolume){
       return (NueConvention::IsInsideFarFiducial_Nue_Extended(x,y,z) > 0);
   }
   return false;
}

// ******************************************************************
//      Single Param Cuts
// ******************************************************************
     
bool NueAnalysisCuts::PassesHotChannel()
{
   if(!IsValid()) return false;

   bool goodEvent = true;
   if(fCheckHotChannel == 0) return goodEvent;
                       
   if(nrInfo){
      if(nrInfo->srevent.hotch == 1) goodEvent = false;
   }else
   if(stInfo){
      if(fDetectorType == Detector::kFar) return true;
      TClonesArray* Strips = ntpManipulator.GetStripArray();
                                                                                
      for(int i=0;i<eventInfo->nstrip;i++){
        NtpSRStrip *strip =
          dynamic_cast<NtpSRStrip *>(Strips->At(eventInfo->stp[i]));               
        if(!strip){
          MSG("NueAnalysisCuts",Msg::kError)<<"Couldn't get strip "
              <<endl;
          return true;
        }
        if(strip->ph0.raw+strip->ph1.raw>60000) goodEvent = false;
      }
   }
                                              
   return goodEvent;
}
   
bool NueAnalysisCuts::PassesLowPhNStripCut()
{
   if(!IsValid()) return false;
                                                                             
   bool goodEvent = true;
   if(LowCuts.shwfit.hiPhStripCount < 0) return goodEvent;
                                                                             
   if(nrInfo){
      if(nrInfo->shwfit.hiPhStripCount < LowCuts.shwfit.hiPhStripCount)
        goodEvent = false;
   }else
   if(stInfo){
      MAXMSG("NueAnalysisCuts",Msg::kError,10)
        <<"PhNStripCut may only be called on ana_nue files"
        <<" the cut will not be applied"<<endl;
   }
                                                                             
   return goodEvent;
}

bool NueAnalysisCuts::PassesLowPhNPlaneCut()
{
   if(!IsValid()) return false;
                                                                             
   bool goodEvent = true;
   if(LowCuts.shwfit.hiPhPlaneCountM2 < 0) return goodEvent;
                                                                             
   if(nrInfo){
      if(nrInfo->shwfit.hiPhPlaneCountM2 < LowCuts.shwfit.hiPhPlaneCountM2)
        goodEvent = false;
   }else
   if(stInfo){
      MAXMSG("NueAnalysisCuts",Msg::kError,10)
        <<"PhNPlaneCutM2 may only be called on ana_nue files"
        <<" the cut will not be applied"<<endl;
   }
                                                                             
   return goodEvent;
}

bool NueAnalysisCuts::PassesPhProngCut()
{
   if(!IsValid()) return false;
                                                                             
   bool goodEvent = true;
   if(fPhProngCut < 0) return goodEvent;
                                                                             
   if(nrInfo){
      float prongPh = TMath::Max(nrInfo->srtrack.pulseHeight,
                          nrInfo->srshower.pulseHeight);
      if(prongPh< fPhProngCut) 
          goodEvent=false;
   }else
   if(stInfo){
     if(showerInfo || trackInfo){
	float shwEnergy = 0;
        float trkEnergy = 0;
	if(showerInfo) shwEnergy = showerInfo->ph.mip;
	if(trackInfo)  trkEnergy = trackInfo->ph.mip;        
        float prongPh = TMath::Max(shwEnergy, trkEnergy);

        if(prongPh< fPhProngCut)  goodEvent=false;
     }else{
	goodEvent = false;
     }                                          
   }
  
  return goodEvent;
}

// ******************************************************************
//      Truth Value Cuts
// ******************************************************************
                                                                                
bool NueAnalysisCuts::IsSignal()
{
   if(!IsValid()) return false;
   if(!IsMC()) return true;
   if(fCutClasses == 0) return true;
                                                                                
   bool goodEvent = false;
   Int_t cType = -10;
                                                                                
   if(nrInfo) cType = nrInfo->mctrue.fNueClass;
   if(mcInfo){
     int inu = mcInfo->inu;
     int inunoosc = mcInfo->inunoosc;
     int iaction = mcInfo->iaction;
     cType = ClassType::DetermineClassType(inu, inunoosc, iaction);
   }
                                                                                
   for(unsigned int i = 0; i < signalClasses.size(); i++){
      if(signalClasses[i] == cType)
          goodEvent = true;
   }
                                                                                
   return goodEvent;
}

bool NueAnalysisCuts::IsBackground()
{
   if(!IsValid()) return false;
   if(!IsMC()) return true;
   if(fCutClasses == 0) return true;
                                                                                
   bool goodEvent = false;
   Int_t cType = -10;
                                                                                
   if(nrInfo) cType = nrInfo->mctrue.fNueClass;
   if(mcInfo){
     int inu = mcInfo->inu;
     int inunoosc = mcInfo->inunoosc;
     int iaction = mcInfo->iaction;
     cType = ClassType::DetermineClassType(inu, inunoosc, iaction);
   }
                                                                                
   for(unsigned int i = 0; i < bgClasses.size(); i++){
      if(bgClasses[i] == cType)
          goodEvent = true;  
   }
                                                                                
   return goodEvent;
}
                                                                                
bool NueAnalysisCuts::PassesCutOnClasses()
{
   if(!IsValid()) return false;
   if(!IsMC()) return true;
   if(fCutClasses == 0) return true;
                                                                                
   bool goodEvent = false;
   if(IsSignal() || IsBackground())
      goodEvent = true;
                                                                                
  return goodEvent;
}

// This function only allows ResCodes in the list through
bool NueAnalysisCuts::PassesResonanceCode()
{
  if(!IsValid()) return false;
  if(!IsMC()) return true;
  if(fCutResCode != 1) return true;

  bool goodEvent = true;

  int code = 0;
  if(nrInfo) code = nrInfo->mctrue.resonanceCode;
  if(mcInfo) code = mcInfo->iresonance;

  if(IsSignal() && sigResonanceCodes.size() > 0)
  {
     goodEvent = false;
     for(unsigned int i = 0; i < sigResonanceCodes.size(); i++) 
     {
	if(code == sigResonanceCodes[i])  goodEvent = true;
     }
  }                                                  
  if(IsBackground() && bgResonanceCodes.size() > 0)
  {
     goodEvent = false;
     for(unsigned int i = 0; i < bgResonanceCodes.size(); i++)
     {
        if(code == bgResonanceCodes[i])  goodEvent = true;
     }
  }

  return goodEvent;
}

bool NueAnalysisCuts::PassesEMFraction()
{
  if(!IsValid()) return false;
  if(!IsMC()) return true;
  bool goodEvent = true;

  goodEvent = PassesHighEMFraction();
  goodEvent = goodEvent && PassesLowEMFraction();

  return goodEvent;
}

bool NueAnalysisCuts::PassesHighEMFraction()
{
  if(!IsValid()) return false;
  if(!IsMC()) return true;
 
  float emfrac = 0.0;
  if(nrInfo) emfrac = nrInfo->mctrue.emShowerFraction;
  if(mcInfo) emfrac = mcInfo->emfrac;

  bool goodEvent = true;
  if(IsSignal())
  {
     if(HighSigCuts.emShowerFraction < 0) return true;
     if(emfrac > HighSigCuts.emShowerFraction)
        goodEvent = false; 
  }
  if(IsBackground())
  {
     if(HighBgCuts.emShowerFraction < 0) return true;
     if(emfrac > HighBgCuts.emShowerFraction)
        goodEvent = false;
  }

  return goodEvent;
}

bool NueAnalysisCuts::PassesLowEMFraction()
{
  if(!IsValid()) return false;
  if(!IsMC()) return true;

  float emfrac = 2.0;
  if(nrInfo) emfrac = nrInfo->mctrue.emShowerFraction;
  if(mcInfo) emfrac = mcInfo->emfrac;
                                                                        
  bool goodEvent = true;
  if(IsSignal())
  {
     if(LowSigCuts.emShowerFraction < 0) return true;
     if(emfrac < LowSigCuts.emShowerFraction)
        goodEvent = false;
  }
  if(IsBackground())
  {
     if(LowBgCuts.emShowerFraction < 0) return true;
     if(emfrac < LowBgCuts.emShowerFraction)
        goodEvent = false;
  }
                                                                             
  return goodEvent;
}

// ******************************************************************
//      Trim To File List 
// ******************************************************************

bool NueAnalysisCuts::PassesFileCut()
{
   if(!IsValid()) return false;
   if(fFileTrim == 0) return true;

   static vector<FilePosition> eventList;
   static unsigned int pos = 0;

   // If this is the first call fill the file list
   if(eventList.size() == 0) {
    std::ifstream ins;
    ins.open(fFileCutName.c_str());
    if(!ins.is_open())
      MSG("NueAnalysisCuts", Msg::kError)<<"Unable to open cut file "
                            <<fFileCutName<<endl;

    while(!ins.eof()){
      int Snarl, Run, SubRun, Event;
      ins>>Run>>SubRun>>Snarl>>Event;
      FilePosition temp;
      if(!ins.eof()){
          temp.Run = Run; temp.SubRun = SubRun; 
          temp.Snarl = Snarl; temp.Event = Event;
          eventList.push_back(temp);
      }
    }
    if(eventList.size() == 0)  return false;
   }
  
   
   FilePosition current;

   if(nrInfo){
     current.Snarl = nrInfo->GetHeader().GetSnarl();
     current.Event = nrInfo->GetHeader().GetEventNo();
     current.Run = nrInfo->GetHeader().GetRun();
     current.SubRun = nrInfo->GetHeader().GetSubRun();
     if(IsMC()) current.SubRun = 0;
   }
  
   if(stInfo){
     current.Snarl = stInfo->GetHeader().GetSnarl();
     current.Event = fEventNum;
     current.Run = stInfo->GetHeader().GetRun();
     current.SubRun = stInfo->GetHeader().GetSubRun();
     if(IsMC()) current.SubRun = 0;
   }

   if(current < eventList[0]) return false;

   if(current > eventList[eventList.size() - 1]) return false;

   if(pos >= eventList.size())
   {
      if(current < eventList[eventList.size() - 1]){
        MSG("NueAnalysisCuts", Msg::kWarning)<<"Problem with file list "
            <<"the file list should always be >= to the position of "
            <<"the files being scanned\n"<<endl;
        return false;
      }
   }

   if(current > eventList[pos] )
   {
      while(current > eventList[pos] && pos < eventList.size())
          pos++;
   }

   if(current == eventList[pos])
   {
       pos++;
       return true;
   }
                                                                                
   if(pos > 0 && current < eventList[pos-1])
   {
      MSG("NueAnalysisCuts", Msg::kWarning)<<"Problem with file list "
            <<" appears to be out of order at position "<<pos<<endl;
        return false;
   }
                                                    
   return false;
}

// ******************************************************************
//      Beam Cuts
// ******************************************************************
                                                                            
bool NueAnalysisCuts::IsGoodBeam()
{
   if(!IsValid()) return false;
   if(IsMC()) return true;
   if(fBeamCut != 1) return true;

   if(nrInfo) 
     return (nrInfo->bmon.goodBeamMon == 1);

   if(eventInfo)
     MSG("NueAnalysisCuts", Msg::kWarning)<<"Beam Cut cannot be used on NtpSt data"<<endl;

   return false;
}

bool NueAnalysisCuts::PassesHornCurrent()
{
   if(!IsValid()) return false;
   if(IsMC()) return true;
   if(fTargetCurrent < 0) return true;
                                                                                
   if(nrInfo){
     Float_t diffCurrent = TMath::Abs(nrInfo->bmon.hornI) - fTargetCurrent;
     bool goodEvent = TMath::Abs(diffCurrent) < 5;
     return goodEvent;
   }
                                                                                
   if(eventInfo)
     MSG("NueAnalysisCuts", Msg::kWarning)<<"Beam Cut cannot be used on NtpSt data"<<endl;
                                                                                
   return false;
}

bool NueAnalysisCuts::PassesDataQuality()
{
  if(!IsValid()) return false;
  if(IsMC()) return true;

  if(nrInfo == 0)
     MSG("NueAnalysisCuts", Msg::kWarning)
        <<"Data Quality Cut cannot be used on NtpSt data"<<endl;
 
  bool goodEvent = true; 
  if(fDetectorType ==  Detector::kFar)
  {
     goodEvent = (nrInfo->srevent.liTime == -1);
     goodEvent = goodEvent && (nrInfo->srevent.rcBoundary == 0);
     goodEvent = goodEvent && (nrInfo->srevent.daveFDDataQuality == 1);
  }
  if(fDetectorType ==  Detector::kNear)
  {
     goodEvent = (nrInfo->srevent.coilCurrent < -1000);
  }

  return goodEvent;
}

//******************************************************************
//   The Reporting Code

void NueAnalysisCuts::ReportOnRecord(NueRecord *nr, string ID)
{
    TClass *cl;
    TRealData *rd;
    TDataMember *member;
    TDataType *membertype;
    Float_t value = 0.0;
    string vName;
                                                                           
    cl=nr->IsA();
    TIter  next(cl->GetListOfRealData());
     
    while ((rd =dynamic_cast<TRealData*>(next()))) {
      member = rd->GetDataMember();
      membertype = member->GetDataType();
      vName=rd->GetName();
      
      Int_t offset = rd->GetThisOffset();
      char *pointer = (char*)nr  + offset;
      value = -9999;
      if(!NeedsSpecialAttention(vName)){
           value=atof(membertype->AsString(pointer));
           if(!ANtpDefVal::IsDefault(value) &&
                !ANtpDefVal::IsDefault(static_cast<Double_t> (value)) &&
                !ANtpDefVal::IsDefault(static_cast<Int_t> (value))){
		     MSG("NueCutReport",Msg::kInfo)<<ID<<"Cut applied on variable "<<vName<<" with value "<<value<<endl;
           } 
      }
   }
}

void NueAnalysisCuts::Report()
{
   MSG("NueCutReport",Msg::kInfo)
       <<"-------------------------- Summary of Cuts ------------------------\n";
   ReportOnRecord(&HighCuts, "High End");
   ReportOnRecord(&LowCuts, "Low End");

   //Then go and take care of the special casses
   MSG("NueCutReport",Msg::kInfo)
       <<"-------------------------- End of Summary -------------------------\n";


}


bool NueAnalysisCuts::NeedsSpecialAttention(TString name)
{                                                                   
   //All the fHeaders and four of hte MST vars require special effort
     if(name == "fHeader.fSnarl"
        ||  name == "fHeader.fRun"
        || name == "fHeader.fSubRun"
        || name == "fHeader.fEvtNo"
        || name == "fHeader.fEvents"
        || name == "fHeader.fTrackLength"
        || name == "mstvars.eallw1"
        || name == "mstvars.oallw1"
        || name == "mstvars.eallm1"
        || name == "mstvars.oallm1")
          return true;

    return false;
}


