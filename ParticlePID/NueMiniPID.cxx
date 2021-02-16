/////////////////////////////////////////////////////////////////////
//$Id: NueMiniPID.cxx,v 1.12 2009/09/21 14:32:26 scavan Exp $
//
//NueMiniPID
//
//J Boehm 12/2007
////////////////////////////////////////////////////////////////////
                                                                                                                                                      
#include "NueAna/ParticlePID/NueMiniPID.h"
#include "MessageService/MsgService.h"
#include "TMath.h"
#include "NueAna/NueRecord.h"

CVSID("$Id: NueMiniPID.cxx,v 1.12 2009/09/21 14:32:26 scavan Exp $");                                                                                                                                                      
ClassImp(NueMiniPID)

//-------------------------------------------------------------------
NueMiniPID::NueMiniPID() :
  fBeam(BeamType::kL010z185i),
  fDet(Detector::kFar),
//  fRelease(ReleaseType::kCedarPhyDaikon)
  fRelease(ReleaseType::kDogwood1Daikon)
{
 fPOT = 0.0;
}

NueMiniPID::NueMiniPID(BeamType::BeamType_t beam,
             Detector::Detector_t det,
             ReleaseType::Release_t rel) :
  fBeam(beam),
  fDet(det),
  fRelease(rel)
{
   evtRecoNueEnergy = 0;
   evtRecoMEUEnergy = 0; ;
   trkRecoCCEnergy = 0;
   nshower = 0;
   contPlanes = 0;
   trkPlanes = 0;
   trkEndPlane = 0;
   trkBegPlane = 0;
   trkLikePlanes = 0;
   shwRecoCCEnergy = 0;
   ann14 = 0;
   ann6 = 0;
   ann30 = 0;
   ann2pe = 0;
   ann2pe_daikon04 = 0;
   ssPID = 0;
   cutPID = 0;
   mcnnPID = 0;
   mcnnMatch = 0;                                                                             
   //truth values
   shiEpi0 = 0;
   shiEmEnergy = 0;
   nuEnergy = 0;
   nuFlavor = 0;
   nonOscNuFlavor = 0;
   nueClass = 0;
   nueOscProb = 0;
   interactionType = 0;
   nuDCosX = 0;
   nuDCosY = 0;
   nuDCosZ = 0;
   hadronicY = 0;
   hadronicFinalState = 0;
   w2 = 0;
   q2 = 0;
   bjorkenX = 0;
   targetPX = 0;
   targetPY = 0;
   targetPZ = 0;
   targetEnergy = 0;
   atomicNumber = 0;
   atomicWeight = 0;
   initialState = 0;
   resonanceCode = 0;
   skzpWeight = 0;
   MCWeight = 0;

   weight = 0;

   endPlaneU = 0;
   endPlaneV = 0;
   deltaUVVtx = 0;

   neugenStdXsec = 0;
   cosmicCut = 0;
   largestEvent = 0;
   
   mri_roCCPID = 0;                                                                             
   mri_abCCPID = 0;
   mri_trkPass = 0;
   gapPlanes = 0;

   fPOT = 0.0;
   
   infid = 0;
   contained = 0;
   ntot = 0;
   longest_s = 0.0;
   event_length = 0.0;
   event_energy = 0.0;
   particle_energy = 0.0;
   
   pidA = 0.0;
   pidB = 0.0;
   pidC = 0.0;
   pidD = 0.0;
   pidE = 0.0;
   pidF = 0.0;
   
   mrcc_s = 0.0;

   pass_var_check = 0;
   pass_nvar_check = 0;
   for(int i=0;i<14;i++)pars[i]=0.0; 




   nueVtxX = 0.0;
   nueVtxY = 0.0;
   nueVtxZ = 0.0;
   vtxU = 0.0;
   vtxV = 0.0;
   vtxZ = 0.0;



   
   run = 0;
   subrun = 0;
   event = 0;
   snarl = 0;
     
   passes_NueStandard_PassesDataQuality = 0; 
   passes_NueStandard_IsInFid = 0;     
   passes_NueStandard_PassesPOTStandards = 0;
   passes_NueStandard_PassesCosmicCut = 0;
   passes_NueStandard_PassesNonHEPreSelection = 0;
   passes_NueStandard_PassesPreSelection = 0;    
   
   passes_NueStandard_PassesMinPlaneCut = 0;
   passes_NueStandard_PassesShowerCut = 0;
   passes_NueStandard_PassesTrackPlaneCut = 0;
   passes_NueStandard_PassesTrackLikePlaneCut = 0;
   passes_NueStandard_PassesLowEnergyCut = 0;
   passes_NueStandard_PassesHighEnergyCut = 0;

   passes_NueStandard_PassesMRCCFiducial = 0;
   passes_NueStandard_PassesMRCCPreSelection = 0;     

   //time stamp
   timestamp=0;      
   
   
   //nue mrcc vars
   mri_qp=0;
   mri_orig_cc_pid=0;
   mri_SigmaQP=0;
	 
   mri_best_complete=0;
   mri_fitp=0;
}

NueMiniPID::~NueMiniPID()
{}

