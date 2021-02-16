/////////////////////////////////////////////////////////////////////
//$Id: NueData.cxx,v 1.12 2009/08/27 15:20:01 scavan Exp $
//
//NueData
//
//J Boehm 12/2007
////////////////////////////////////////////////////////////////////
                                                                                                                                                      
#include "NueAna/Extrapolation/NueData.h"
#include "MessageService/MsgService.h"
#include "TMath.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NueMini.h"
#include "NueAna/NueMiniAna.h"

CVSID("$Id: NueData.cxx,v 1.12 2009/08/27 15:20:01 scavan Exp $");                                                                                                                                                      
//-------------------------------------------------------------------
NueData::NueData() :
  fBeam(BeamType::kL010z185i),
  fDet(Detector::kFar),
  //  fRelease(ReleaseType::kCedarPhyDaikon)
  fRelease(ReleaseType::kDogwood1Daikon)
{
 fPOT = 0.0;
 kIsMR = false;
 fNM = new NueMini(fBeam, fDet, fRelease);
}

NueData::NueData(BeamType::BeamType_t beam,
             Detector::Detector_t det,
             ReleaseType::Release_t rel, bool isNueData) :
  fBeam(beam),
  fDet(det),
  fRelease(rel),
  kIsNueData(isNueData)
{
  fPOT = 0.0;
  kIsMR = false;
  fNM = new NueMini(fBeam, fDet, fRelease);
}

void NueData::SetupNueHeader(NueHeader &nh)
{
   const VldTimeStamp vts;
   const VldContext Vld(fDet, SimFlag::kMC, vts);
   nh.SetVldContext(Vld);
   nh.SetRelease(fRelease);
   nh.SetBeamType(fBeam);
   nh.SetFoundBits(0,0,0,kIsMR);
}

void NueData::AddEvent(NueRecord *nr)
{
   if(!nr){
     MSG("NueData", Msg::kError) << "No NueRecord Passed"<<endl;
     return;
   }

   if(nr->GetHeader().GetRelease() != fRelease){
     MSG("NueData", Msg::kError) << "ReleaseType MisMatch" 
            << nr->GetHeader().GetRelease()<<" != "<<fRelease<<endl;
   }
   if(nr->GetHeader().GetBeamType() != fBeam){
     MSG("NueData", Msg::kError) << "BeamType MisMatch"
            << nr->GetHeader().GetBeamType()<<" != "<<fBeam<<endl;
   }
   if(nr->GetHeader().GetVldContext().GetDetector() != fDet){
     MSG("NueData", Msg::kError) << "Detector MisMatch"
            << nr->GetHeader().GetVldContext().GetDetector() 
            <<" != "<<fDet<<endl;
   }

   fNMA.FillMini(nr, fNM);
   AddEvent(fNM);

/* 
   evtRecoNueEnergy.push_back(nr->srevent.phNueGeV);
   evtRecoMEUEnergy.push_back(nr->srevent.phMip);

   int tP = 0;    int tLP = 0;     float trkECC = 0;
   if(nr->srtrack.phCCGeV > -10){
     trkECC = nr->srtrack.phCCGeV;
     tP = TMath::Abs(nr->srtrack.endPlane - nr->srtrack.begPlane);
     tLP = nr->srtrack.trklikePlanes;
   }
   trkRecoCCEnergy.push_back(trkECC);
   trkPlanes.push_back(tP);
   trkEndPlane.push_back(nr->srtrack.endPlane);
   trkBegPlane.push_back(nr->srtrack.begPlane);
   trkLikePlanes.push_back(tLP);

   nshower.push_back(nr->srevent.showers);
   contPlanes.push_back(nr->shwfit.contPlaneCount050);

   if(nr->srshower.phCCGeV > -10)
     shwRecoCCEnergy.push_back(nr->srshower.phCCGeV);
   else shwRecoCCEnergy.push_back(0);

   ann2pe.push_back(nr->ann.pid_11inp);
   ann30.push_back(nr->ann.pid_30inp);
   ann6.push_back(nr->ann.pid_6inp);

   cutPID.push_back(nr->treepid.fCutPID);
   ssPID.push_back(nr->subshowervars.pid);
   mcnnPID.push_back(nr->mcnnv.mcnn_var1);
   mcnnMatch.push_back(nr->mcnnv.bestmatches);

   //truth values
   shiEpi0.push_back(nr->shi.epi0);
   shiEmEnergy.push_back(nr->shi.emenergy);
   nuEnergy.push_back(nr->mctrue.nuEnergy);
   nuFlavor.push_back(nr->mctrue.nuFlavor);
   nonOscNuFlavor.push_back(nr->mctrue.nonOscNuFlavor);
   nueClass.push_back(nr->mctrue.fNueClass);

   nueOscProb.push_back(nr->mctrue.fOscProb);
   interactionType.push_back(nr->mctrue.interactionType);
   nuDCosX.push_back(nr->mctrue.nuDCosX);
   nuDCosY.push_back(nr->mctrue.nuDCosY);
   nuDCosZ.push_back(nr->mctrue.nuDCosZ);
   hadronicY.push_back(nr->mctrue.hadronicY);
   hadronicFinalState.push_back(nr->mctrue.hadronicFinalState);
   w2.push_back(nr->mctrue.w2);
   q2.push_back(nr->mctrue.q2);
   bjorkenX.push_back(nr->mctrue.bjorkenX);
   targetPX.push_back(nr->mctrue.targetPX);
   targetPY.push_back(nr->mctrue.targetPY);
   targetPZ.push_back(nr->mctrue.targetPZ);
   targetEnergy.push_back(nr->mctrue.targetEnergy);
   atomicNumber.push_back(nr->mctrue.atomicNumber);
   atomicWeight.push_back(nr->mctrue.atomicWeight);
   initialState.push_back(nr->mctrue.initialState);
   resonanceCode.push_back(nr->mctrue.resonanceCode);

   skzpWeight.push_back(nr->fluxweights.totbeamweight);


   ntrack.push_back(nr->srevent.tracks);
   trkPass.push_back(nr->srtrack.passedFit);
   endPlaneU.push_back(nr->srtrack.endPlaneU);
   endPlaneV.push_back(nr->srtrack.endPlaneV);
   deltaUVVtx.push_back(nr->srtrack.deltaUVVtx);
   abCCPID.push_back(nr->anainfo.abCCPID);
   roCCPID.push_back(nr->anainfo.roCCPID);

   mri_roCCPID.push_back(nr->mri.orig_roCCPID);
   mri_abCCPID.push_back(nr->mri.orig_abCCPID);
   mri_trkPass.push_back(nr->mri.fitp);
   gapPlanes.push_back(nr->srshower.gapPlanes);
   neugenStdXsec.push_back(1);
   cosmicCut.push_back((int) nr->dtree.bt_var1);
   largestEvent.push_back(nr->srevent.largestEvent);
*/
}

void NueData::AddEvent(NueMini *nm)
{
   if(!nm){
     MSG("NueData", Msg::kError) << "No NueRecord Passed"<<endl;
     return;
   }
                                                                                                                            
   if(nm->fRelease != fRelease){
     MSG("NueData", Msg::kError) << "ReleaseType MisMatch"
            << nm->fRelease<<" != "<<fRelease<<endl;
   }
   if(nm->fBeam != fBeam){
     MSG("NueData", Msg::kError) << "BeamType MisMatch"
            << nm->fBeam<<" != "<<fBeam<<endl;
   }
   if(nm->fDet != fDet){
     MSG("NueData", Msg::kError) << "Detector MisMatch"
            << nm->fDet <<" != "<<fDet<<endl;
   }                                                                                                                          
   evtRecoNueEnergy.push_back(nm->evtRecoNueEnergy);
   evtRecoMEUEnergy.push_back(nm->evtRecoMEUEnergy);
   trkRecoCCEnergy.push_back(nm->trkRecoCCEnergy);
   trkPlanes.push_back(nm->trkPlanes);
   trkEndPlane.push_back(nm->trkEndPlane);
   trkBegPlane.push_back(nm->trkBegPlane);
   trkLikePlanes.push_back(nm->trkLikePlanes);
   shwRecoCCEnergy.push_back(nm->shwRecoCCEnergy);

   nshower.push_back(nm->nshower);
   contPlanes.push_back(nm->contPlanes);

   //ann2pe_daikon04.push_back(nm->ann2pe_daikon04);
   //ann2pe.push_back(nm->ann2pe);
   
   printf("we've adjusted the NueMini structure (changed ann var names)... but didn't fix NueData.cxx... do that now... (ask Steve)\n");
   exit(1);
   
   //ann30.push_back(nm->ann30);
   //ann6.push_back(nm->ann6);
   ssPID.push_back(nm->ssPID);
   //cutPID.push_back(nm->cutPID);                                                                                                                            
   mcnnPID.push_back(nm->mcnnPID);
   mcnnMatch.push_back(nm->mcnnMatch);

   //truth values
   shiEpi0.push_back(nm->shiEpi0);
   shiEmEnergy.push_back(nm->shiEmEnergy);
   nuEnergy.push_back(nm->nuEnergy);
   nuFlavor.push_back(nm->nuFlavor);
   nonOscNuFlavor.push_back(nm->nonOscNuFlavor);
   nueClass.push_back(nm->nueClass);
                                                                                                                            
   nueOscProb.push_back(nm->nueOscProb);
   interactionType.push_back(nm->interactionType);
   nuDCosX.push_back(nm->nuDCosX);
   nuDCosY.push_back(nm->nuDCosY);
   nuDCosZ.push_back(nm->nuDCosZ);
   hadronicY.push_back(nm->hadronicY);
   hadronicFinalState.push_back(nm->hadronicFinalState);
   w2.push_back(nm->w2);
   q2.push_back(nm->q2);
   bjorkenX.push_back(nm->bjorkenX);
   targetPX.push_back(nm->targetPX);
   targetPY.push_back(nm->targetPY);
   targetPZ.push_back(nm->targetPZ);
   targetEnergy.push_back(nm->targetEnergy);
   atomicNumber.push_back(nm->atomicNumber);
   atomicWeight.push_back(nm->atomicWeight);
   initialState.push_back(nm->initialState);
   resonanceCode.push_back(nm->resonanceCode);
                
//the meaning of skzpWeight might have changed in nuemini!

//it seems to have been RPWeight....
// now it is actually skzpWeight .....
// RPWeight is now nm->RPWeight                                                                                                            
   //skzpWeight.push_back(nm->skzpWeight);
	exit(1);
   ntrack.push_back(nm->ntrack);
   trkPass.push_back(nm->trkPass);
   endPlaneU.push_back(nm->endPlaneU);
   endPlaneV.push_back(nm->endPlaneV);
   deltaUVVtx.push_back(nm->deltaUVVtx);
   abCCPID.push_back(nm->abCCPID);
   roCCPID.push_back(nm->roCCPID);
   mri_roCCPID.push_back(nm->mri_orig_roCCPID);
   mri_abCCPID.push_back(nm->mri_orig_abCCPID);
   mri_trkPass.push_back(nm->mri_trkPass);   
   gapPlanes.push_back(nm->gapPlanes);
   neugenStdXsec.push_back(1);
   cosmicCut.push_back((int) nm->cosmicCut);
   largestEvent.push_back(nm->largestEvent);

}



void NueData::FillRecord(NueRecord *nr, int i)
{  
   if(!nr){
     MSG("NueData", Msg::kError) << "No NueRecord Passed"<<endl;
     return;
   }

   if(nr->GetHeader().GetVldContext().GetDetector() != fDet){
     MSG("NueData", Msg::kError) << "Detector MisMatch"
            << nr->GetHeader().GetVldContext().GetDetector()
            <<" != "<<fDet<<endl;
   }
   if(nr->GetHeader().GetRelease() != fRelease){
     MSG("NueData", Msg::kError) << "ReleaseType MisMatch"
            << nr->GetHeader().GetRelease()<<" != "<<fRelease<<endl;
   }
   if(nr->GetHeader().GetBeamType() != fBeam){
     MSG("NueData", Msg::kError) << "BeamType MisMatch"
            << nr->GetHeader().GetBeamType()<<" != "<<fBeam<<endl;
   }

   nr->srevent.phNueGeV = evtRecoNueEnergy[i];
   nr->srevent.phMip = evtRecoMEUEnergy[i];
   nr->srtrack.phCCGeV = trkRecoCCEnergy[i];
   nr->srtrack.trklikePlanes = trkLikePlanes[i];
   nr->srtrack.begPlane = trkBegPlane[i];
   nr->srtrack.endPlane = trkEndPlane[i];
   nr->srshower.phCCGeV = shwRecoCCEnergy[i];
   nr->srevent.showers = nshower[i];
   nr->shwfit.contPlaneCount050 = contPlanes[i];
   nr->ann.pid_11inp_daikon04 = ann2pe_daikon04[i];
   nr->ann.pid_11inp = ann2pe[i];
                                                                                
   nr->ann.pid_30inp = ann30[i];
   nr->ann.pid_6inp = ann6[i];

   nr->subshowervars.pid = ssPID[i];
   nr->treepid.fCutPID = cutPID[i];   
   nr->mcnnv.mcnn_var1 = mcnnPID[i];
   nr->mcnnv.bestmatches = mcnnMatch[i];

                                                                             
   //truth values
   nr->shi.epi0 = shiEpi0[i];
   nr->shi.emenergy = shiEmEnergy[i];
   nr->mctrue.nuEnergy = nuEnergy[i];
   nr->mctrue.nuFlavor = nuFlavor[i];
   nr->mctrue.nonOscNuFlavor = nonOscNuFlavor[i];
   nr->mctrue.fNueClass = nueClass[i];
                                                                                
   nr->mctrue.fOscProb = nueOscProb[i];
   nr->mctrue.interactionType = interactionType[i];
   nr->mctrue.nuDCosX = nuDCosX[i];
   nr->mctrue.nuDCosY = nuDCosY[i];
   nr->mctrue.nuDCosZ = nuDCosZ[i];
   nr->mctrue.hadronicY = hadronicY[i];
   nr->mctrue.hadronicFinalState = hadronicFinalState[i];
   nr->mctrue.w2 = w2[i];
   nr->mctrue.q2 = q2[i];
   nr->mctrue.bjorkenX = bjorkenX[i];
   nr->mctrue.targetPX = targetPX[i];
   nr->mctrue.targetPY = targetPY[i];
   nr->mctrue.targetPZ = targetPZ[i];
   nr->mctrue.targetEnergy = targetEnergy[i];
   nr->mctrue.atomicNumber = atomicNumber[i];
   nr->mctrue.atomicWeight = atomicWeight[i];
   nr->mctrue.initialState = initialState[i];
   nr->mctrue.resonanceCode = resonanceCode[i];
                                                                                
   nr->fluxweights.totbeamweight = skzpWeight[i];

   nr->srevent.tracks = ntrack[i];
   nr->srtrack.passedFit = trkPass[i];
   nr->srtrack.endPlaneU = endPlaneU[i];
   nr->srtrack.endPlaneV = endPlaneV[i];
   nr->srtrack.deltaUVVtx = deltaUVVtx[i];
   nr->anainfo.abCCPID = abCCPID[i];
   nr->xsecweights.xsecweight = neugenStdXsec[i];
   nr->mri.orig_abCCPID = mri_abCCPID[i];
   nr->mri.fitp = mri_trkPass[i];
   nr->srshower.gapPlanes = gapPlanes[i];
   nr->eventq.passCosmicCut = 1; // cosmicCut[i]; this is a hack for now due to cheer problems
   nr->srevent.largestEvent = largestEvent[i];

   nr->anainfo.roCCPID = roCCPID[i];
   nr->mri.orig_roCCPID = mri_roCCPID[i];

}

double NueData::GetNeugenStdXsec(int i)
{
   return neugenStdXsec[i];
}
void   NueData::SetNeugenStdXsec(double xsec, int i){
  neugenStdXsec[i] = xsec;
}

void NueData::Clear()
{
  trkRecoCCEnergy.clear();
  shwRecoCCEnergy.clear();
  evtRecoNueEnergy.clear();
  evtRecoMEUEnergy.clear();
  skzpWeight.clear();
  weight.clear();
  trkPlanes.clear();
  trkEndPlane.clear();
  trkBegPlane.clear();
  trkLikePlanes.clear();
  shiEpi0.clear();
  shiEmEnergy.clear();
  nuEnergy.clear();
  nuFlavor.clear();
  nonOscNuFlavor.clear();
  nueOscProb.clear();
  nueClass.clear();
  ParentType.clear();
  interactionType.clear();
  nuDCosX.clear();
  nuDCosY.clear();
  nuDCosZ.clear();
  hadronicY.clear();
  hadronicFinalState.clear();
  w2.clear();
  q2.clear();
  bjorkenX.clear();
  targetPX.clear();
  targetPY.clear();
  targetPZ.clear();
  targetEnergy.clear();
  atomicNumber.clear();
  atomicWeight.clear();
  initialState.clear();
  resonanceCode.clear();
  ann30.clear();
  ann6.clear();
  ssPID.clear();
  mcnnPID.clear();
  mcnnMatch.clear();
  cutPID.clear();
  abCCPID.clear();
  ntrack.clear();
  trkPass.clear();
  endPlaneU.clear();
  endPlaneV.clear();
  deltaUVVtx.clear();
  mri_abCCPID.clear();
  mri_trkPass.clear();
  gapPlanes.clear();
  neugenStdXsec.clear();
  cosmicCut.clear();
  largestEvent.clear();

  contPlanes.clear();
  ann2pe_daikon04.clear();
  ann2pe.clear();
  nshower.clear();
  roCCPID.clear();
  mri_roCCPID.clear();


  vector<double>(trkRecoCCEnergy).swap(trkRecoCCEnergy);
  vector<double>(shwRecoCCEnergy).swap(shwRecoCCEnergy);
  vector<double>(evtRecoNueEnergy).swap(evtRecoNueEnergy);
  vector<double>(evtRecoMEUEnergy).swap(evtRecoMEUEnergy);
  vector<double>(skzpWeight).swap(skzpWeight);
  vector<double>(weight).swap(weight);
  vector<int>(trkPlanes).swap(trkPlanes);
  vector<int>(trkEndPlane).swap(trkEndPlane);
  vector<int>(trkBegPlane).swap(trkBegPlane);
  vector<int>(trkLikePlanes).swap(trkLikePlanes);      
  vector<double>(shiEpi0).swap(shiEpi0);
  vector<double>(shiEmEnergy).swap(shiEmEnergy);
  vector<double>(nuEnergy).swap(nuEnergy);
  vector<int>(nuFlavor).swap(nuFlavor);
  vector<int>(nonOscNuFlavor).swap(nonOscNuFlavor);
  vector<double>(nueOscProb).swap(nueOscProb);
  vector<int>(nueClass).swap(nueClass);
  vector<double>(ParentType).swap(ParentType);
  vector<int>(interactionType).swap(interactionType);
  vector<double>(nuDCosX).swap(nuDCosX);
  vector<double>(nuDCosY).swap(nuDCosY);
  vector<double>(nuDCosZ).swap(nuDCosZ);
  vector<double>(hadronicY).swap(hadronicY);
  vector<int>(hadronicFinalState).swap(hadronicFinalState);

  vector<double>(w2).swap(w2);
  vector<double>(q2).swap(q2);
  vector<double>(bjorkenX).swap(bjorkenX);
  vector<double>(targetPX).swap(targetPX);
  vector<double>(targetPY).swap(targetPY);
  vector<double>(targetPZ).swap(targetPZ);
  vector<double>(targetEnergy).swap(targetEnergy);
  vector<double>(atomicNumber).swap(atomicNumber);
  vector<double>(atomicWeight).swap(atomicWeight);
  vector<int>(initialState).swap(initialState);
  vector<int>(resonanceCode).swap(resonanceCode);
  vector<double>(ann30).swap(ann30);
  vector<double>(ann6).swap(ann6);
  vector<double>(ssPID).swap(ssPID);
  vector<double>(mcnnPID).swap(mcnnPID);
  vector<int>(mcnnMatch).swap(mcnnMatch);
  vector<int>(cutPID).swap(cutPID);
  vector<double>(abCCPID).swap(abCCPID);
  vector<int>(ntrack).swap(ntrack);
  vector<int>(trkPass).swap(trkPass);
  vector<int>(endPlaneU).swap(endPlaneU);
  vector<int>(endPlaneV).swap(endPlaneV);
  vector<int>(deltaUVVtx).swap(deltaUVVtx);
  vector<double>(mri_abCCPID).swap(mri_abCCPID);
  vector<int>(mri_trkPass).swap(mri_trkPass);
  vector<int>(gapPlanes).swap(gapPlanes);
  vector<double>(neugenStdXsec).swap(neugenStdXsec);
  vector<int>(cosmicCut).swap(cosmicCut);
  vector<int>(largestEvent).swap(largestEvent);

  vector<int>(contPlanes).swap(contPlanes); 
  vector<double>(ann2pe_daikon04).swap(ann2pe_daikon04); 
  vector<double>(ann2pe).swap(ann2pe); 
  vector<int>(nshower).swap(nshower); 
  vector<double>(roCCPID).swap(roCCPID); 
  vector<double>(mri_roCCPID).swap(mri_roCCPID); 
 



}

