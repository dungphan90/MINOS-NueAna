/**
 *
 * $Id: AnalysisInfoAna.cxx,v 1.24 2009/10/13 16:06:03 scavan Exp $
 *
 * \class AnalysisInfoAna
 *
 * \package NueAna
 *
 * \brief
 *
 * Contact: J. Boehm
 *
 * Created on: Mon Apr 18 17:17:36 2005
 *
 */
#include <cassert>

#include "TSystem.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "CandNtupleSR/NtpSREvent.h"
#include "MessageService/MsgService.h"
#include "StandardNtuple/NtpStRecord.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "NueAna/AnalysisInfoAna.h"
#include "NueAna/NueAnalysisCuts.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "Mad/MadNsID.h"
#include "Mad/MadDpID.h"
#include "Mad/MadAbID.h"
#include "AnalysisNtuples/Module/ANtpRecoNtpManipulator.h"
#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"
#include "DataUtil/infid.h"

#include "NueAna/Module/SetKNNModule.h"  
#include "PhysicsNtuple/Handle.h"
#include "PhysicsNtuple/Factory.h"
  

                                                                                
MadNsID AnalysisInfoAna::nsid;
MadDpID AnalysisInfoAna::dpid;
MadAbID AnalysisInfoAna::abid;

bool AnalysisInfoAna::readabidfile=false;
                                                                                
CVSID("$Id: AnalysisInfoAna.cxx,v 1.24 2009/10/13 16:06:03 scavan Exp $");

AnalysisInfoAna::AnalysisInfoAna(AnalysisInfoNue &anai):
   fAnalysisInfo(anai)
{
    fDetectorType =  Detector::kUnknown;
}
                                                                                
AnalysisInfoAna::~AnalysisInfoAna()
{}
                                                                                
void AnalysisInfoAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
  MSG("AnalysisInfoAna",Msg::kDebug)
      <<"Entering AnalysisInfoAna::Analyze"<<endl;
                                                                                
  if(srobj==0){
      return;
  }  if(((dynamic_cast<NtpStRecord *>(srobj))==0)&&
        ((dynamic_cast<NtpSRRecord *>(srobj))==0)){
     return;
  }
                                                                                
  const RecCandHeader *ntpHeader = &(srobj->GetHeader());
  fDetectorType = ntpHeader->GetVldContext().GetDetector();
  
  //Not applicable to CalDet:
  if(fDetectorType==Detector::kCalDet) return;
  
  bool foundst=false;

  ANtpRecoNtpManipulator *ntpManipulator = 0;
                                                                                
  NtpSRRecord *sr = 0;
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(srobj);

  NtpSREvent *event = 0;
  NtpSRTrack *track = 0;
  NtpSRShower *shower = 0;
  NtpSREventSummary *eventSummary;
  SimFlag::SimFlag_t sim;

  if(st != 0)
  {
      foundst = true;
      ntpManipulator =  new ANtpRecoNtpManipulator(st);
      eventSummary = &(st->evthdr);
      sim = st->GetHeader().GetVldContext().GetSimFlag();
  }
  else{
      sr=dynamic_cast<NtpSRRecord *>(srobj);
      ntpManipulator =  new ANtpRecoNtpManipulator(sr);
      eventSummary = &(sr->evthdr);
      sim = sr->GetHeader().GetVldContext().GetSimFlag();
  }

  ntpManipulator->SetPrimaryShowerCriteria(0,1); //nplanes, total ph
  ntpManipulator->SetPrimaryTrackCriteria(0,1,0); //nplanes, length, total ph
  ANtpEventManipulator * ntpEventManipulator =
                            ntpManipulator->GetEventManipulator();
                                                                                
  ntpEventManipulator->SetEventInSnarl(evtn);
  event  = ntpEventManipulator->GetEvent();
  shower = SntpHelpers::GetPrimaryShower(evtn, st);
  track  = SntpHelpers::GetPrimaryTrack(evtn,st);

  Float_t vertexX = event->vtx.x; //Munits::meters
  Float_t vertexY = event->vtx.y; //Munits::meters
  Float_t vertexZ = event->vtx.z; //Munits::meters
                                                                               
  if(ReleaseType::IsCedar(release)){
    NtpVtxFinder vtxf;
    vtxf.SetTargetEvent(st, event, track);
    if(vtxf.FindVertex() > 0){
       vertexX = vtxf.VtxX();
       vertexY = vtxf.VtxY();
       vertexZ = vtxf.VtxZ();
    }
  }

                    
  fAnalysisInfo.isFullyContained =  
     IsFidAll(vertexX, vertexY, vertexZ, event);

  Int_t test = 0;
  bool isMC = (sim == SimFlag::kMC);
  if(fDetectorType == Detector::kNear)
      test = NueConvention::IsInsideNearFiducial_Nue_Standard(vertexX, vertexY, vertexZ, isMC);
  if(fDetectorType == Detector::kFar)
      test = NueConvention::IsInsideFarFiducial_Nue_Standard(vertexX, vertexY, vertexZ, isMC);
  
  fAnalysisInfo.inFiducialVolume = test;


  if(fDetectorType==Detector::kFar ){
       dpid.SetPHCorrection(1.018);
  }

  BeamType::BeamType_t current_beam = beam;
  //  BeamType::BeamType_t current_beam = BeamType::kL250z200i;

  string reco_version = ReleaseType::AsString(release); 
  string mc_version = ""; 
  if(ntpHeader->GetVldContext().GetSimFlag()==SimFlag::kMC) 
    mc_version = ReleaseType::AsString(release); 
  if(ReleaseType::IsCarrot(release)) reco_version = "birch"; 
  else if(ReleaseType::IsCedar(release)) reco_version = "cedar"; 
  else if(ReleaseType::IsDogwood(release)) reco_version = "cedar";  //<-- ! no dogwood file! 
  else reco_version = "birch"; 


  if(dpid.ChoosePDFs(fDetectorType,current_beam,
		     reco_version,mc_version))
    fAnalysisInfo.dpCCPID = 
          dpid.CalcPID(track, event, eventSummary, fDetectorType, 0);

   
  if(nsid.ChooseWeightFile(fDetectorType,current_beam)){
    if(!nsid.GetPID(event,track,shower,st,fDetectorType,fAnalysisInfo.nsCCPID))
       fAnalysisInfo.nsCCPID=ANtpDefVal::kFloat;
  }

  if(!readabidfile){
    string fname = AnalysisInfoAna::BuildABPIDFile();
    abid.ReadPDFs(fname.c_str());    
    MSG("AnalysisInfoAna",Msg::kInfo)<<"Reading AB pdfs from "<<fname<<endl;
    readabidfile=true;
  }

  //now that everything's all set up, compute andy's pid
  NtpStRecord *strecord = dynamic_cast<NtpStRecord *>(srobj);
  if(strecord==0){}
  else{    
    fAnalysisInfo.abCCPID=abid.CalcPID(event,strecord);
    fAnalysisInfo.abcostheta=abid.GetRecoCosTheta(event,strecord);
    fAnalysisInfo.abeventenergy=abid.GetRecoE(event,strecord);
    fAnalysisInfo.abtrackcharge=abid.GetTrackEMCharge(event,strecord);
    fAnalysisInfo.abtrackenergy=abid.GetTrackPlanes(event,strecord);
    fAnalysisInfo.abtracklikeplanes=abid.GetTrackLikePlanes(event,strecord);
    fAnalysisInfo.abtrackphfrac=abid.GetTrackPHfrac(event,strecord);
    fAnalysisInfo.abtrackphmean=abid.GetTrackPHmean(event,strecord);
    fAnalysisInfo.abtrackqpsigmaqp=abid.GetTrackQPsigmaQP(event,strecord);
    fAnalysisInfo.aby=abid.GetRecoY(event,strecord);
  }

  //rustem's pid
   //Loading up Rustems variables
   Anp::Handle<StorekNNData> data = Anp::Factory<StorekNNData>::Instance().Get("kNNData");
   if(!data.valid())
   {
      MAXMSG("NueModule", Msg::kError, 1)
          << "NueModule::Reco - Handle<StorekNNData> is invalid, assuming no Rustem variable to run"
          << endl;
   }else{

     float numubar, rel_ang, knn_pid, knn_01, knn_10, knn_20, knn_40; 

     data -> SetPrefix("SNTP");
     data -> Get(evtn, "numubar", numubar);
     data -> Get(evtn, "rel_ang", rel_ang);
     data -> Get(evtn, "knn_pid", knn_pid);
     data -> Get(evtn, "knn_01", knn_01);
     data -> Get(evtn, "knn_10", knn_10);
     data -> Get(evtn, "knn_20", knn_20);
     data -> Get(evtn, "knn_40", knn_40);

     fAnalysisInfo.roNScientPlanes = knn_01;
     fAnalysisInfo.roMeanTrkSig = knn_10;
     fAnalysisInfo.roMeanRatio = knn_20;
     fAnalysisInfo.roTrkMeanVsWindow = knn_40;
     fAnalysisInfo.roNuMuBar = numubar;
     fAnalysisInfo.roRelAng = rel_ang;
     fAnalysisInfo.roCCPID  = knn_pid;
  }
 

  if(ntpManipulator){
    delete ntpManipulator;
    ntpManipulator=0;
  }
    
}

void AnalysisInfoAna::Set3DHit(DeqFloat_t &x
                             , DeqFloat_t &y
                             , DeqFloat_t &z
                             , DeqFloat_t &e){
                                                                                
     if(x.size()!=0 && y.size()!=0 && z.size()!=0 && e.size()!=0){
                                                                                
           fX=x;
           fY=y;
           fZ=z;
           fE=e;
     }
}

Int_t AnalysisInfoAna::IsFidAll(Float_t vtxX, Float_t vtxY, Float_t vtxZ, NtpSREvent *event)
{
  // with 3d hits shut off this function won't do anythnig interesting
                                                                                
  Bool_t contained = true;
  Int_t test = 0;
  if(fX.size()==0 || fY.size()==0 || fZ.size()==0 || fE.size()==0){
        return ANtpDefVal::kInt;
  }
                                                                                
  for(UInt_t i = 0; i < fX.size() && contained; i++){
     test = 0;
     Float_t x = fX[i] + vtxX;
     Float_t y = fY[i] + vtxY;
     Float_t z = fZ[i] + vtxZ;
     if(fDetectorType == Detector::kNear)
      test = NueConvention::IsInsideNearFiducial_Nue_Extended(x,y,z);
     if(fDetectorType == Detector::kFar)
      test = NueConvention::IsInsideFarFiducial_Nue_Extended(x,y,z);      
     if(test <= 0) contained = false;
  }
                                                                                
  if(event)
  if(contained && event->plane.n > 66) //this is approximately 4 meteres which
// is the 3D hit cutoff
   {
      Float_t x = event->end.x;
      Float_t y = event->end.y;
      Float_t z = event->end.z;

     if(fDetectorType == Detector::kNear)
      test = NueConvention::IsInsideNearFiducial_Nue_Extended(x,y,z);
     if(fDetectorType == Detector::kFar)
      test = NueConvention::IsInsideFarFiducial_Nue_Extended(x,y,z);
     if(test <= 0) contained = false;
   }
                                                                                
  if(contained) test = 1;
  return test;
}

string AnalysisInfoAna::BuildABPIDFile()
{
    string base="";
    string pub = getenv("SRT_PUBLIC_CONTEXT");
    string priv = getenv("SRT_PRIVATE_CONTEXT");

    if(priv!=""&&priv!="."){
      // check if directory exists in SRT_PRIVATE_CONTEXT
      std::string path = priv + "/Mad/data";
      void *dir_ptr = gSystem -> OpenDirectory(path.c_str());

      if(!dir_ptr){
        base=pub; // if it doesn't exist use SRT_PUBLIC
      }
      else base = priv;
    }
    else{
      base=pub;
    }
    if(base=="") {
      MSG("AnalysisInfoAna",Msg::kFatal)<<"No SRT_PUBLIC_CONTEXT set "
                                        <<"Do not know where to look "
                                        <<"for AB pdf files "<<std::endl;
      assert(false);
    }
    base+="/Mad/data";
    string fname=base;
    string rmc="";
    if(ReleaseType::IsCedar(release)&&ReleaseType::IsDaikon(release)){
      rmc="cedar_daikon";
    }
    else if(ReleaseType::IsCedar(release)&&ReleaseType::IsData(release)){
      rmc="cedar_daikon";
    }
    else{
      MSG("AnalysisInfoAna",Msg::kWarning)<<"Dont know reco/mc versions "
                                          <<"defaulting to cedar_daikon"<<endl;
      rmc="cedar_daikon";
    }
    string sdet="";
    if(fDetectorType==Detector::kNear){
      sdet="near";
    }
    else if(fDetectorType==Detector::kFar){
      sdet="far";
    }
    else{
      MSG("AnalysisInfoAna",Msg::kWarning)<<"Dont know detector type "
                                          <<"defaulting to far"<<endl;
      sdet="far";
    }

    string sbeam="";
    switch (beam){
    case BeamType::kL010z000i:     sbeam="le0";      break;
    case BeamType::kL010z170i:     sbeam="le170";    break;
    case BeamType::kL010z185i:     sbeam="le";       break;
    case BeamType::kL010z200i:     sbeam="le200";    break;
    case BeamType::kL100z200i:     sbeam="pme";      break;
    case BeamType::kL150z200i:     sbeam="pme";      break;
    case BeamType::kL250z200i:     sbeam="phe";      break;
    case BeamType::kLE:            sbeam="le";       break;
    default:
      MSG("AnalysisInfoAna",Msg::kWarning)<<"Don't know beam type "
                                          <<" defaulting to LE"<<endl;
      sbeam="le";
      break;
    }
                                                                                
    fname+="/ab_pdf_"+sdet+"_"+sbeam+"_"+rmc+".root";
 
    return fname; 
}

string AnalysisInfoAna::BuildROPIDFile()
{
    string base="";
    base=getenv("SRT_PRIVATE_CONTEXT");
    if(base!=""&&base!="."){
      // check if directory exists in SRT_PRIVATE_CONTEXT
      std::string path = base + "/Mad/data";
      void *dir_ptr = gSystem -> OpenDirectory(path.c_str());
      if(!dir_ptr){
        base=getenv("SRT_PUBLIC_CONTEXT"); // if not there use SRT_PUBLIC_CONTEXT
      }
    }
    else{
      base=getenv("SRT_PUBLIC_CONTEXT");
    }
    if(base=="") {
      MSG("AnalysisInfoAna",Msg::kFatal)<<"No SRT_PUBLIC_CONTEXT set "
                                        <<"Do not know where to look "
                                        <<"for RO pdf files "<<std::endl;
      assert(false);
    }
    base+="/Mad/data";
    //fname will depend on reco/mc version and detector, but not on beam
    string rmc="";
    if(ReleaseType::IsCedar(release)&&ReleaseType::IsDaikon(release)){
      rmc="cedar.daikon";
    }
    else if(ReleaseType::IsCedar(release)&&ReleaseType::IsData(release)){
      rmc="cedar.daikon";
    }
    else {
      MSG("AnalysisInfoAna",Msg::kWarning)<<"Dont know reco/mc versions "
                             <<"defaulting to cedar_daikon "<<endl;
      rmc="cedar.daikon";
    }
    string sdet="";
    if(fDetectorType==Detector::kNear){
      sdet="near";
    }
    else if(fDetectorType==Detector::kFar){
      sdet="far";
    }
    else{
      MSG("AnalysisInfoAna",Msg::kWarning)<<"Dont know detector type "
                                          <<"defaulting to far"<<endl;
      sdet="far";
    }
    string fname=base;
    fname+="/knn.train."+sdet+"."+rmc+".root";
                                                                                
    return fname;
}
