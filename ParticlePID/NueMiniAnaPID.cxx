/////////////////////////////////////////////////////////////////////
//$Id: NueMiniAnaPID.cxx,v 1.10 2009/09/21 14:32:26 scavan Exp $
//
//NueMiniPID
//
//J Boehm 12/2007
////////////////////////////////////////////////////////////////////
                                                                                                                                                      
#include "NueAna/ParticlePID/NueMiniPID.h"
#include "NueAna/ParticlePID/NueMiniAnaPID.h"

#include "MessageService/MsgService.h"
#include "TMath.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NueStandard.h"
#include "NueAna/ParticlePID/ParticleAna/PRecord.h"
#include "NueAna/NueHeader.h"


CVSID("$Id: NueMiniAnaPID.cxx,v 1.10 2009/09/21 14:32:26 scavan Exp $");
                                                                                

//-------------------------------------------------------------------
void NueMiniAnaPID::FillMini(NueRecord *nr, NueMiniPID *nm, int domrcc)
{
   if(!nr){
     MSG("NueMiniAnaPID", Msg::kError) << "No NueRecord Passed"<<endl;
     return;
   }

   if(nr->GetHeader().GetRelease() != nm->fRelease){
     MSG("NueMiniAnaPID", Msg::kError) << "ReleaseType MisMatch: " 
            << nr->GetHeader().GetRelease()<<" != "<<nm->fRelease<<endl;
   }
   if(nr->GetHeader().GetBeamType() != nm->fBeam){
     MSG("NueMiniAnaPID", Msg::kError) << "BeamType MisMatch: "
            << nr->GetHeader().GetBeamType()<<" != "<<nm->fBeam<<endl;
   }
   if(nr->GetHeader().GetVldContext().GetDetector() != nm->fDet){
     MSG("NueMiniAnaPID", Msg::kError) << "Detector MisMatch: "
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

   nm->ann14 = nr->ann.pid;
   nm->ann30 = nr->ann.pid_30inp;
   nm->ann6 = nr->ann.pid_6inp;
   nm->ann2pe = nr->ann.pid_11inp;
   nm->ann2pe_daikon04 = nr->ann.pid_11inp_daikon04;

   nm->ssPID = nr->subshowervars.pid;
   nm->cutPID = nr->treepid.fCutPID;
   nm->mcnnPID = nr->mcnnv.mcnn_var1;
   nm->mcnnMatch = nr->mcnnv.bestmatches;

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

   nm->skzpWeight = nr->fluxweights.totbeamweight;
   nm->MCWeight = NueStandard::GetMCWeights(nr); 

   nm->ntrack =  nr->srevent.tracks;
   nm->trkPass = nr->srtrack.passedFit;
   nm->endPlaneU = nr->srtrack.endPlaneU;
   nm->endPlaneV = nr->srtrack.endPlaneV;
   nm->deltaUVVtx = nr->srtrack.deltaUVVtx;
   nm->abCCPID = nr->anainfo.abCCPID;
   nm->roCCPID = nr->anainfo.roCCPID;

   nm->mri_roCCPID = nr->mri.orig_roCCPID;
   nm->mri_abCCPID = nr->mri.orig_abCCPID;
   nm->mri_trkPass = nr->mri.fitp;
   nm->gapPlanes = nr->srshower.gapPlanes;
   nm->neugenStdXsec = 1;
   nm->cosmicCut = (int) nr->eventq.passCosmicCut;
   nm->largestEvent = nr->srevent.largestEvent;
   
   PRecord *pr = & nr->precord;
   if(domrcc)
   	pr=& nr->precordMRCC;
   	
   
   nm->longest_s=pr->particles.longest_s_particle_s;
   nm->event_length=pr->event.max_z-pr->event.min_z;
   nm->ntot=pr->particles.ntot;
   nm->infid=pr->event.inFiducial == 1;
   nm->contained=pr->event.contained == 1;
   nm->event_energy = pr->event.visenergy;
   nm->particle_energy = pr->particles.totvise;
   
   
   nm->pidA=pr->event.pidA;
   nm->pidB=pr->event.pidB;
   nm->pidC=pr->event.pidC;
   nm->pidD=pr->event.pidD;
   nm->pidE=pr->event.pidE;
   nm->pidF=pr->event.pidF;
   
   if(domrcc)
   		nm->mrcc_s = pr->mrccinfo.hasMRCC ? pr->mrccinfo.particle_s : 0;
  




   bool pass=1;

   if(pr->particles.longest_s_particle_s<0 || pr->particles.longest_s_particle_s>6)pass=0;
   if(pr->particles.longest_z<0 || pr->particles.longest_z>6)pass=0;
   if(pr->particles.ntot<0 || pr->particles.ntot>50)pass=0;
   if(pr->particles.rms_r<0 || pr->particles.rms_r>100)pass=0;
   if(pr->particles.prim_par_e0<0 || pr->particles.prim_par_e0>40e3)pass=0;
   if(pr->particles.prim_par_chisq<0 || pr->particles.prim_par_chisq>1000)pass=0;
   if(pr->particles.largest_particle_peakdiff<-200 || pr->particles.largest_particle_peakdiff>200)pass=0;
   if(nm->event_length<0 || nm->event_length>6)pass=0;

   nm->pass_var_check=pass;
    
   pass=1;
   double mstvar_combine = nr->mstvars.e4w+nr->mstvars.o4w;
   if(nr->shwfit.par_b<0 || nr->shwfit.par_b>6)pass=0;
   if(mstvar_combine<-1000 || mstvar_combine>1000)pass=0;
   if(nr->shwfit.LongE<0 || nr->shwfit.LongE>1200)pass=0;
   nm->pass_nvar_check=pass;


   double largest_frac=pr->particles.totvise ?
        pr->particles.largest_particle_e/pr->particles.totvise : 1;
   double prim_ae0=pr->particles.prim_par_e0 ?
        pr->particles.prim_par_a/pr->particles.prim_par_e0 : 0;
   double largest_cmp_chisqndf = pr->particles.largest_particle_cmp_ndf ?
        pr->particles.largest_particle_cmp_chisq /
            pr->particles.largest_particle_cmp_ndf : 0;


   int z=0;
   nm->pars[z++]=pr->particles.longest_s_particle_s;
   nm->pars[z++]=pr->particles.mol_rad_r;
   nm->pars[z++]=pr->particles.emfrac;
   nm->pars[z++]=pr->particles.ntot;
   nm->pars[z++]=pr->particles.weighted_phi;
   nm->pars[z++]=largest_frac;
   nm->pars[z++]=pr->particles.prim_par_b;
   nm->pars[z++]=pr->particles.prim_par_e0;
   nm->pars[z++]=pr->particles.prim_par_chisq;
   nm->pars[z++]=pr->particles.largest_particle_peakdiff;
   nm->pars[z++]=largest_cmp_chisqndf;
   nm->pars[z++]=pr->particles.prim_par_a;
   nm->pars[z++]=pr->event.nclusters;
   nm->pars[z++]=prim_ae0;



   nm->vtxU=pr->event.vtx_u;
   nm->vtxV=pr->event.vtx_v;
   nm->vtxZ=pr->event.vtx_z;

   nm->nueVtxX=nr->srevent.vtxX;
   nm->nueVtxY=nr->srevent.vtxY;
   nm->nueVtxZ=nr->srevent.vtxZ;

 
   
   
   
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



   nm->timestamp = nr->GetHeader().GetVldContext().GetTimeStamp().GetSec();


   nm->mri_qp=nr->mri.qp;
   nm->mri_orig_cc_pid=nr->mri.orig_cc_pid;
   nm->mri_SigmaQP=nr->mri.SigmaQP;
	 
   nm->mri_best_complete=nr->mri.best_complete;
   nm->mri_fitp=nr->mri.fitp;

 
}

void NueMiniAnaPID::FillRecord(NueRecord *nr, NueMiniPID *nm, int domrcc)
{  
   if(!nr || !nm){
     MSG("NueMiniPID", Msg::kError) << "No NueRecord Passed"<<endl;
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
       
   nr->ann.pid= nm->ann14;                              
   nr->ann.pid_11inp = nm->ann2pe;                                           
   nr->ann.pid_11inp_daikon04 = nm->ann2pe_daikon04;                                           
   nr->ann.pid_30inp = nm->ann30;
   nr->ann.pid_6inp = nm->ann6;
   nr->subshowervars.pid = nm->ssPID;
   nr->treepid.fCutPID = nm->cutPID;
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

   nr->mri.orig_roCCPID = nm->mri_roCCPID;
   nr->mri.orig_abCCPID = nm->mri_abCCPID;
   nr->mri.fitp = nm->mri_trkPass;
   nr->srshower.gapPlanes = nm->gapPlanes;
   nr->eventq.passCosmicCut = nm->cosmicCut;
   nr->srevent.largestEvent = nm->largestEvent;
   
   PRecord *pr = & nr->precord;
   if(domrcc)
   	pr=& nr->precordMRCC;
   	
   
   pr->particles.longest_s_particle_s=nm->longest_s;

   //we don't enough have information here!
   pr->event.max_z=nm->event_length;
   pr->event.min_z=0;
   
   pr->particles.ntot=nm->ntot;
   pr->event.inFiducial=nm->infid;
   pr->event.contained=nm->contained;
   pr->event.visenergy=nm->event_energy;
   pr->particles.totvise=nm->particle_energy;
   
   pr->event.pidA=nm->pidA;   
   pr->event.pidB=nm->pidB;
   pr->event.pidC=nm->pidC; 
   pr->event.pidD=nm->pidD;
   pr->event.pidE=nm->pidE;   
   pr->event.pidF=nm->pidF;
   
   if(domrcc)
   		pr->mrccinfo.particle_s = nm->mrcc_s;



   pr->event.vtx_u=nm->vtxU;
   pr->event.vtx_v=nm->vtxV;
   pr->event.vtx_z=nm->vtxZ;

   nr->srevent.vtxX=nm->nueVtxX;
   nr->srevent.vtxY=nm->nueVtxY;
   nr->srevent.vtxZ=nm->nueVtxZ;

   NueHeader * nh = const_cast<NueHeader *> (&(nr->GetHeader()));
   nh->SetRun(nm->run);
   nh->SetSubRun(nm->subrun);   
   nh->SetEventNo(nm->event);
   nh->SetSnarl(nm->snarl);   
   
   

   nr->mri.qp=nm->mri_qp;
   nr->mri.orig_cc_pid=nm->mri_orig_cc_pid;
   nr->mri.SigmaQP=nm->mri_SigmaQP;
	 
   nr->mri.best_complete=nm->mri_best_complete;
   nr->mri.fitp=nm->mri_fitp;   
   
}

