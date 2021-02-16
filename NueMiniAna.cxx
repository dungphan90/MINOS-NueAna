/////////////////////////////////////////////////////////////////////
//$Id: NueMiniAna.cxx,v 1.15 2009/09/13 22:16:01 jjling Exp $
//
//NueMini
//
//J Boehm 12/2007
////////////////////////////////////////////////////////////////////
                                                                                                                                                      
#include "NueAna/NueMini.h"
#include "NueAna/NueMiniAna.h"

#include "MessageService/MsgService.h"
#include "TMath.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NueStandard.h"
#include "NueAna/NueHeader.h"


CVSID("$Id: NueMiniAna.cxx,v 1.15 2009/09/13 22:16:01 jjling Exp $");
                                                                                

//-------------------------------------------------------------------
void NueMiniAna::FillMini(NueRecord *nr, NueMini *nm)
{
   if(!nr){
     MSG("NueMiniAna", Msg::kError) << "No NueRecord Passed"<<endl;
     return;
   }

   if(nr->GetHeader().GetRelease() != nm->fRelease){
     MSG("NueMiniAna", Msg::kError) << "ReleaseType MisMatch: " 
            << nr->GetHeader().GetRelease()<<" != "<<nm->fRelease<<endl;
   }
   if(nr->GetHeader().GetBeamType() != nm->fBeam){
     MSG("NueMiniAna", Msg::kError) << "BeamType MisMatch: "
            << nr->GetHeader().GetBeamType()<<" != "<<nm->fBeam<<endl;
   }
   if(nr->GetHeader().GetVldContext().GetDetector() != nm->fDet){
     MSG("NueMiniAna", Msg::kError) << "Detector MisMatch: "
            << nr->GetHeader().GetVldContext().GetDetector() 
            <<" != "<<nm->fDet<<endl;
   }

   nm->evtRecoNueEnergy = nr->srevent.phNueGeV;
   nm->evtRecoMEUEnergy = nr->srevent.phMip;

   int tP = 0;    int tLP = 0;     float trkECC = 0;
   if(nr->srtrack.phCCGeV > -10){
     trkECC = nr->srtrack.phCCGeV;
     tP = TMath::Abs(nr->srtrack.endPlane - nr->srtrack.begPlane);
     tLP = nr->srtrack.trklikePlanes;
   }
   nm->trkRecoCCEnergy = trkECC;
   nm->trkPlanes = tP;
   nm->trkEndPlane = nr->srtrack.endPlane;
   nm->trkBegPlane = nr->srtrack.begPlane;
   nm->trkLikePlanes = tLP;
   nm->nshower = nr->srevent.showers;
   nm->contPlanes = nr->shwfit.contPlaneCount050;
   

   if(nr->srshower.phCCGeV > -10)
     nm->shwRecoCCEnergy = nr->srshower.phCCGeV;
   else nm->shwRecoCCEnergy = 0;

   nm->annpid_11inp = nr->ann.pid_11inp;
   nm->annpid_11inp_daikon04 = nr->ann.pid_11inp_daikon04;
   nm->annpid_14inp_daikon04 = nr->ann.pid_14inp_daikon04;
   nm->annpid = nr->ann.pid;


   nm->ssPID = nr->subshowervars.pid;


   //LEM variables
   nm->mcnnv_meanfracQmatched = nr->mcnnv.meanfracQmatched;
   nm->mcnnv_var2 = nr->mcnnv.mcnn_var2;
   nm->mcnnv_fracCCy = nr->mcnnv.fracCCy;
   nm->mcnnv_qtot = nr->mcnnv.qtot;
   nm->mcnnPID = nr->mcnnv.mcnn_var1;
   nm->mcnnMatch = nr->mcnnv.bestmatches;
 
   nm->mcnnv_fracCC = nr->mcnnv.fracCC;
   nm->mcnnv_ymean = nr->mcnnv.ymean;
   nm->mcnnv_meanfracQmatched = nr->mcnnv.meanfracQmatched;

   //truth values
   nm->shiEpi0 = nr->shi.epi0;
   nm->shiEmEnergy = nr->shi.emenergy;
   nm->nuEnergy = nr->mctrue.nuEnergy;
   nm->nuFlavor = nr->mctrue.nuFlavor;
   nm->nonOscNuFlavor = nr->mctrue.nonOscNuFlavor;
   nm->nueClass = nr->mctrue.fNueClass;

   nm->nueOscProb = nr->mctrue.fOscProb;
   nm->interactionType = nr->mctrue.interactionType;
   nm->nuDCosX = nr->mctrue.nuDCosX;
   nm->nuDCosY = nr->mctrue.nuDCosY;
   nm->nuDCosZ = nr->mctrue.nuDCosZ;
   nm->hadronicY = nr->mctrue.hadronicY;
   nm->hadronicFinalState = nr->mctrue.hadronicFinalState;
   nm->w2 = nr->mctrue.w2;
   nm->q2 = nr->mctrue.q2;
   nm->bjorkenX = nr->mctrue.bjorkenX;
   nm->targetPX = nr->mctrue.targetPX;
   nm->targetPY = nr->mctrue.targetPY;
   nm->targetPZ = nr->mctrue.targetPZ;
   nm->targetEnergy = nr->mctrue.targetEnergy;
   nm->atomicNumber = nr->mctrue.atomicNumber;
   nm->atomicWeight = nr->mctrue.atomicWeight;
   nm->initialState = nr->mctrue.initialState;
   nm->resonanceCode = nr->mctrue.resonanceCode;

   nm->skzpWeight = nr->fluxweights.totskzpweight;
   

   nm->ntrack =  nr->srevent.tracks;
   nm->trkPass = nr->srtrack.passedFit;
   nm->endPlaneU = nr->srtrack.endPlaneU;
   nm->endPlaneV = nr->srtrack.endPlaneV;
   nm->deltaUVVtx = nr->srtrack.deltaUVVtx;
   nm->abCCPID = nr->anainfo.abCCPID;
   nm->roCCPID = nr->anainfo.roCCPID;

   nm->mri_orig_roCCPID = nr->mri.orig_roCCPID;
   nm->mri_orig_abCCPID = nr->mri.orig_abCCPID;
   nm->mri_trkPass = nr->mri.fitp;
   nm->gapPlanes = nr->srshower.gapPlanes;
   nm->neugenStdXsec = 1;
   nm->cosmicCut = (int) nr->eventq.passCosmicCut;
   nm->largestEvent = nr->srevent.largestEvent;
   
   
   
   
   nm->run=nr->GetHeader().GetRun();
   nm->subrun=nr->GetHeader().GetSubRun();   
   nm->event=nr->GetHeader().GetEventNo();
   nm->snarl=nr->GetHeader().GetSnarl();
   
   nm->passes_NueStandard_PassesDataQuality=NueStandard::PassesDataQuality(nr); 
   nm->passes_NueStandard_IsInFid=NueStandard::IsInFid(nr);     
   nm->passes_NueStandard_PassesPOTStandards=NueStandard::PassesPOTStandards(nr);
   nm->passes_NueStandard_PassesCosmicCut=NueStandard::PassesCosmicCut(nr);
   nm->passes_NueStandard_PassesNonHEPreSelection=NueStandard::PassesNonHEPreSelection(nr);
   nm->passes_NueStandard_PassesPreSelection=NueStandard::PassesPreSelection(nr);    
   
   nm->passes_NueStandard_PassesMinPlaneCut=NueStandard::PassesMinPlaneCut(nr);
   nm->passes_NueStandard_PassesShowerCut=NueStandard::PassesShowerCut(nr);
   nm->passes_NueStandard_PassesTrackPlaneCut=NueStandard::PassesTrackPlaneCut(nr);
   nm->passes_NueStandard_PassesTrackLikePlaneCut=NueStandard::PassesTrackLikePlaneCut(nr);
   nm->passes_NueStandard_PassesLowEnergyCut=NueStandard::PassesLowEnergyCut(nr);
   nm->passes_NueStandard_PassesHighEnergyCut=NueStandard::PassesHighEnergyCut(nr);

   nm->passes_NueStandard_PassesMRCCFiducial=NueStandard::PassesMRCCFiducial(nr);
   nm->passes_NueStandard_PassesMRCCPreSelection=NueStandard::PassesMRCCPreSelection(nr);     
   
   
   
   
   

   //MRCC variables
   nm->mri_best_complete_phw = nr->mri.best_complete_phw;
   nm->mri_best_purity_phw = nr->mri.best_purity_phw;
   nm->mri_shwe = nr->mri.shwe;
   nm->mri_pmux = nr->mri.pmux;
   nm->mri_pmuy = nr->mri.pmuy;
   nm->mri_pmuz = nr->mri.pmuz;
   nm->mri_vtxx = nr->mri.vtxx;
   nm->mri_vtxy = nr->mri.vtxy;
   nm->mri_vtxz = nr->mri.vtxz;





     //ANN 11 Variables:
   nm->shwfit_par_a = nr->shwfit.par_a;
   nm->shwfit_par_b = nr->shwfit.par_b;
   nm->shwfit_uv_molrad_peak_9s_2pe_dw = nr->shwfit.uv_molrad_peak_9s_2pe_dw;
   nm->shwfit_uv_rms_9s_2pe_dw = nr->shwfit.uv_rms_9s_2pe_dw;
   nm->mstvars_e4w_PLUS_mstvars_o4w = nr->mstvars.e4w + nr->mstvars.o4w;
   nm->fracvars_fract_road = nr->fracvars.fract_road;
   nm->fracvars_fract_2_planes = nr->fracvars.fract_2_planes;
   nm->fracvars_fract_4_planes = nr->fracvars.fract_4_planes;
   nm->fracvars_fract_6_planes = nr->fracvars.fract_6_planes;
   nm->fracvars_fract_8_counters = nr->fracvars.fract_8_counters;
   nm->shwfit_LongE = nr->shwfit.LongE;


   // ANN14 Variables, new added:
   nm->fracvars_shw_max = nr->fracvars.shw_max;
   nm->fracvars_shw_slp = nr->fracvars.shw_slp;
   nm->fracvars_dis2stp = nr->fracvars.dis2stp;
   nm->fracvars_fract_1_plane = nr->fracvars.fract_1_plane;
   nm->fracvars_fract_5_planes = nr->fracvars.fract_5_planes;
   nm->fracvars_fract_6_counters = nr->fracvars.fract_6_counters;
   nm->fracvars_fract_20_counters = nr->fracvars.fract_20_counters; 
   nm->event_length = TMath::Abs(nr->srevent.endPlane - nr->srevent.begPlane);


     //SSPID variables:
   nm->ssvar1 = nr->subshowervars.PHFracRMSU + nr->subshowervars.PHFracRMSV;
   nm->ssvar2 = ((nr->subshowervars.PHAvgProbEMU * nr->subshowervars.PHAvgProbEMU)+
   				(nr->subshowervars.PHAvgProbEMV * nr->subshowervars.PHAvgProbEMV))/2.;
   nm->ssvar3 = (nr->subshowervars.PHAvgIDU + nr->subshowervars.PHAvgIDV)/2.;
   nm->ssvar4 = nr->subshowervars.PHAvgDevU + nr->subshowervars.PHAvgDevV;

   nm->vtxX = nr->srevent.vtxX;
   nm->vtxY = nr->srevent.vtxY;
   nm->vtxZ = nr->srevent.vtxZ;
   
   
   nm->timestamp = nr->GetHeader().GetVldContext().GetTimeStamp().GetSec();
   
   
   
   
}

void NueMiniAna::FillRecord(NueRecord *nr, NueMini *nm)
{  
   if(!nr || !nm){
     MSG("NueMini", Msg::kError) << "No NueRecord Passed"<<endl;
     return;
   }

   nr->srevent.phNueGeV = nm->evtRecoNueEnergy;
   nr->srevent.phMip = nm->evtRecoMEUEnergy;
   nr->srtrack.phCCGeV = nm->trkRecoCCEnergy;
   nr->srtrack.trklikePlanes = nm->trkLikePlanes;
   nr->srtrack.begPlane = nm->trkBegPlane;
   nr->srtrack.endPlane = nm->trkEndPlane;
   nr->srshower.phCCGeV = nm->shwRecoCCEnergy;

   nr->srevent.showers = nm->nshower; 
   nr->shwfit.contPlaneCount050 = nm->contPlanes;
                                     

   
   nr->ann.pid_11inp = nm->annpid_11inp;
   nr->ann.pid_11inp_daikon04 = nm->annpid_11inp_daikon04;
   nr->ann.pid = nm->annpid;   
                                             
   nr->subshowervars.pid = nm->ssPID;
   nr->mcnnv.mcnn_var1 = nm->mcnnPID;
   nr->mcnnv.bestmatches = nm->mcnnMatch;
                                                                                
   //truth values
   nr->shi.epi0 = nm->shiEpi0;
   nr->shi.emenergy = nm->shiEmEnergy;
   nr->mctrue.nuEnergy = nm->nuEnergy;
   nr->mctrue.nuFlavor = nm->nuFlavor;
   nr->mctrue.nonOscNuFlavor = nm->nonOscNuFlavor;
   nr->mctrue.fNueClass = nm->nueClass;
                                                                                
   nr->mctrue.fOscProb = nm->nueOscProb;
   nr->mctrue.interactionType = nm->interactionType;
   nr->mctrue.nuDCosX = nm->nuDCosX;
   nr->mctrue.nuDCosY = nm->nuDCosY;
   nr->mctrue.nuDCosZ = nm->nuDCosZ;
   nr->mctrue.hadronicY = nm->hadronicY;
   nr->mctrue.hadronicFinalState = nm->hadronicFinalState;
   nr->mctrue.w2 = nm->w2;
   nr->mctrue.q2 = nm->q2;
   nr->mctrue.bjorkenX = nm->bjorkenX;
   nr->mctrue.targetPX = nm->targetPX;
   nr->mctrue.targetPY = nm->targetPY;
   nr->mctrue.targetPZ = nm->targetPZ;
   nr->mctrue.targetEnergy = nm->targetEnergy;
   nr->mctrue.atomicNumber = nm->atomicNumber;
   nr->mctrue.atomicWeight = nm->atomicWeight;
   nr->mctrue.initialState = nm->initialState;
   nr->mctrue.resonanceCode = nm->resonanceCode;
                                                                                
   nr->fluxweights.totbeamweight = nm->skzpWeight;

   nr->srevent.tracks = nm->ntrack;
   nr->srtrack.passedFit = nm->trkPass;
   nr->srtrack.endPlaneU = nm->endPlaneU;
   nr->srtrack.endPlaneV = nm->endPlaneV;
   nr->srtrack.deltaUVVtx = nm->deltaUVVtx;
   nr->anainfo.abCCPID = nm->abCCPID;
   nr->anainfo.roCCPID = nm->roCCPID;

   nr->mri.orig_roCCPID = nm->mri_orig_roCCPID;
   nr->mri.orig_abCCPID = nm->mri_orig_abCCPID;
   nr->mri.fitp = nm->mri_trkPass;
   nr->srshower.gapPlanes = nm->gapPlanes;
   nr->eventq.passCosmicCut = nm->cosmicCut;
   nr->srevent.largestEvent = nm->largestEvent;
   
   
   NueHeader * nh = const_cast<NueHeader *> (&(nr->GetHeader()));
   nh->SetRun(nm->run);
   nh->SetSubRun(nm->subrun);   
   nh->SetEventNo(nm->event);
   nh->SetSnarl(nm->snarl);   
   
   
   
   //MRCC variables
   nr->mri.best_complete_phw = nm->mri_best_complete_phw;
   nr->mri.best_purity_phw = nm->mri_best_purity_phw;
   nr->mri.shwe = nm->mri_shwe;
   nr->mri.pmux = nm->mri_pmux;
   nr->mri.pmuy = nm->mri_pmuy;
   nr->mri.pmuz = nm->mri_pmuz;
   nr->mri.vtxx = nm->mri_vtxx;
   nr->mri.vtxy = nm->mri_vtxy;
   nr->mri.vtxz = nm->mri_vtxz;





   //ANN 11 Variables:
   nr->shwfit.par_a = nm->shwfit_par_a;
   nr->shwfit.par_b = nm->shwfit_par_b;
   nr->shwfit.uv_molrad_peak_9s_2pe_dw = nm->shwfit_uv_molrad_peak_9s_2pe_dw;
   nr->shwfit.uv_rms_9s_2pe_dw = nm->shwfit_uv_rms_9s_2pe_dw;
   nr->mstvars.e4w = nm->mstvars_e4w_PLUS_mstvars_o4w;
   nr->fracvars.fract_road = nm->fracvars_fract_road;
   nr->fracvars.fract_2_planes = nm->fracvars_fract_2_planes;
   nr->fracvars.fract_4_planes = nm->fracvars_fract_4_planes;
   nr->fracvars.fract_6_planes = nm->fracvars_fract_6_planes;
   nr->fracvars.fract_8_counters = nm->fracvars_fract_8_counters;
   nr->shwfit.LongE  =nm->shwfit_LongE;


   //LEM variables:
   nr->mcnnv.fracCC = nm->mcnnv_fracCC;
   nr->mcnnv.ymean = nm->mcnnv_ymean;
   nr->mcnnv.meanfracQmatched = nm->mcnnv_meanfracQmatched; 

     //SSPID variables:
	//We can't fill these back in...




   nr->srevent.vtxX = nm->vtxX;
   nr->srevent.vtxY = nm->vtxY;
   nr->srevent.vtxZ = nm->vtxZ;


   //time stamp
   //can't reload the time stamp...
  

   
}

