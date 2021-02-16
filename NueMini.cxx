/////////////////////////////////////////////////////////////////////
//$Id: NueMini.cxx,v 1.14 2009/09/13 22:15:44 jjling Exp $
//
//NueMini
//
//J Boehm 12/2007
////////////////////////////////////////////////////////////////////
                                                                                                                                                      
#include "NueAna/NueMini.h"
#include "MessageService/MsgService.h"
#include "TMath.h"
#include "NueAna/NueRecord.h"

CVSID("$Id: NueMini.cxx,v 1.14 2009/09/13 22:15:44 jjling Exp $");                                                                                                                                                      
ClassImp(NueMini)

//-------------------------------------------------------------------
NueMini::NueMini() :
  fBeam(BeamType::kL010z185i),
  fDet(Detector::kFar),
//  fRelease(ReleaseType::kCedarPhyDaikon)
  fRelease(ReleaseType::kDogwood1Daikon)
{
 fPOT = 0.0;
}

NueMini::NueMini(BeamType::BeamType_t beam,
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
   annpid_11inp = 0;
   annpid_11inp_daikon04 = 0;
   annpid_14inp_daikon04 = 0;
   annpid = 0;
   ssPID = 0;

   //LEM variables:
   mcnnv_fracCC = 0;
   mcnnv_ymean = 0;                                                                mcnnv_meanfracQmatched = 0;
   mcnnPID = 0;
   mcnnMatch = 0;                                                                  mcnnv_var2 = 0;
   mcnnv_fracCCy = 0;
   mcnnv_qtot = 0;

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
 

   endPlaneU = 0;
   endPlaneV = 0;
   deltaUVVtx = 0;

   neugenStdXsec = 0;
   cosmicCut = 0;
   largestEvent = 0;
   
   mri_orig_roCCPID = 0;                                                                             
   mri_orig_abCCPID = 0;
   mri_trkPass = 0;
   gapPlanes = 0;

   fPOT = 0.0;
   
   
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
   
   
   //MRCC variables   
   mri_best_complete_phw = 0; 
   mri_best_purity_phw = 0; 
   mri_shwe = 0; 
   mri_pmux = 0; 
   mri_pmuy = 0; 
   mri_pmuz = 0; 
   mri_vtxx = 0; 
   mri_vtxy = 0; 
   mri_vtxz = 0; 





   //ANN 11 Variables:
   shwfit_par_a = 0; 
   shwfit_par_b = 0; 
   shwfit_uv_molrad_peak_9s_2pe_dw = 0; 
   shwfit_uv_rms_9s_2pe_dw = 0; 
   mstvars_e4w_PLUS_mstvars_o4w = 0; 
   fracvars_fract_road = 0; 
   fracvars_fract_2_planes = 0; 
   fracvars_fract_4_planes = 0; 
   fracvars_fract_6_planes = 0; 
   fracvars_fract_8_counters = 0; 
   shwfit_LongE = 0; 

  
   //some new added ANN14 Variables:
   
   fracvars_shw_max = 0;
   fracvars_shw_slp = 0;
   fracvars_dis2stp = 0;
   fracvars_fract_1_plane = 0;
   fracvars_fract_5_planes = 0;
   fracvars_fract_6_counters = 0;
   fracvars_fract_20_counters = 0; 
   event_length = 0;

   //SSPID variables:
   ssvar1 = 0; 
   ssvar2 = 0; 
   ssvar3 = 0; 
   ssvar4 = 0; 
   
   //event vertex:
   vtxX = 0;
   vtxY = 0;
   vtxZ = 0;
   
   //time stamp
   timestamp=0;   
   
   
}

NueMini::~NueMini()
{}

