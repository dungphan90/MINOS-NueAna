#include "NueAna/Extrapolation/NueExtrapHelper.h"
#include "NueAna/NueStandard.h"
#include "TMath.h"
#include <iostream>
#include "NueAna/NueAnaTools/NueConvention.h"

NueExtrapHelper::NueExtrapHelper() :
  fNXBins(0),fXBins(0),fNYBins(0),fYBins(0),fRecord(0),fMini(0)
{
  this->Init();
}

NueExtrapHelper::NueExtrapHelper(Int_t nx,Double_t lx,Double_t ux,
				 Int_t ny,Double_t ly,Double_t uy) :
  fNXBins(nx),fXBins(0),fNYBins(ny),fYBins(0),fRecord(0),fMini(0)
{
  fXBins = new Double_t[fNXBins+1];
  Float_t bwidth = (ux-lx)/float(fNXBins);
  for(int i=0;i<fNXBins+1;i++) fXBins[i] = lx + float(i)*bwidth;
  if(fNYBins>0) {
    fYBins = new Double_t[fNYBins+1];
    bwidth = (uy-ly)/float(fNYBins);
    for(int i=0;i<fNYBins+1;i++) fYBins[i] = ly + float(i)*bwidth;
  }
  this->Init();
}

NueExtrapHelper::NueExtrapHelper(Int_t nx,Double_t *x,
				 Int_t ny,Double_t *y) :
  fNXBins(nx),fXBins(x),fNYBins(ny),fYBins(y),fRecord(0),fMini(0)
{
  this->Init();
}

void NueExtrapHelper::Init()
{
  fNearChain = NULL;
  fFarChain = NULL;
  fNearPOT = 0.;
  fFarPOT = 0.; 
  fCurSel = Selection::kUnknown;
}

NueExtrapHelper::~NueExtrapHelper()
{
  delete [] fXBins;
  delete [] fYBins;
}

Bool_t NueExtrapHelper::PassBasicCuts()
{
  if(fRecord->GetHeader().GetEventNo()<0) return false;
  if(fRecord->anainfo.inFiducialVolume!=1) return false;
  if(fRecord->srtrack.planes>=25) return false;  
  if(fRecord->srtrack.trklikePlanes>=18) return false;
  //if(fRecord->srevent.phMeu>150) return false;
  //if(TMath::Max(fRecord->srtrack.pulseHeight,
  //	fRecord->srshower.pulseHeight)<5000) return false;  
  if(fRecord->GetHeader().GetVldContext().GetDetector()==Detector::kNear && 
     fRecord->anainfo.isFullyContained==-2) return true;
  if(fRecord->anainfo.isFullyContained!=1) return false;
  return true;
}

Bool_t NueExtrapHelper::PassCuts(NueRecord *nr, Selection::Selection_t sel)
{
  return NueStandard::PassesSelection(nr, sel);
}

Bool_t NueExtrapHelper::PassCuts(Selection::Selection_t sel)
{
  return PassCuts(fRecord, sel);
}

void NueExtrapHelper::AddNueSystematic(NueSystematic */*nueSys*/)
{
}

void NueExtrapHelper::MakeHelpers(Selection::Selection_t sel)
{
  fCurSel = sel;
}

void NueExtrapHelper::WriteFile(std::string /*fname*/)
{
}

void NueExtrapHelper::SetChains(TChain *nearChain, TChain *farChain,
				Double_t nearPOT, Double_t farPOT)
{
  fNearChain = nearChain;
  fFarChain = farChain;
  SetUpNueAnaChain(fNearChain);
  SetUpNueAnaChain(fFarChain);
  fNearPOT = nearPOT;
  fFarPOT = farPOT;
}

void NueExtrapHelper::SetUpNueMiniChain(TChain *chain)
{
  if(!fMini) fMini = new NueMini();
  chain->SetBranchAddress("NueMini",&fMini);
}

void NueExtrapHelper::SetUpNueAnaChain(TChain *chain)
{
  if(!fRecord) fRecord = new NueRecord();
  chain->SetBranchAddress("NueRecord",&fRecord);

//  chain->SetBranchStatus("dtree*", 0);
  chain->SetBranchStatus("mcnnv*", 0);
  chain->SetBranchStatus("timing*", 0);
  chain->SetBranchStatus("cdi*", 0);
  chain->SetBranchStatus("mri*", 0);
  chain->SetBranchStatus("mri*", 0);
  chain->SetBranchStatus("treepid*", 0);
  chain->SetBranchStatus("highhit*", 0);
  chain->SetBranchStatus("mstvar*", 0);
  chain->SetBranchStatus("fracvars*", 0);
  chain->SetBranchStatus("shield*", 0);
  chain->SetBranchStatus("angcluster*", 0);
  chain->SetBranchStatus("shwfit*", 0);
  chain->SetBranchStatus("hitcalc*", 0);
  chain->SetBranchStatus("anainfo*", 0);
  chain->SetBranchStatus("mda*", 0);
  chain->SetBranchStatus("srshower*", 0);
  chain->SetBranchStatus("srshower.phNueGeV", 1);
  chain->SetBranchStatus("srshower.phCCGeV",  1);
  chain->SetBranchStatus("srtrack*", 0);
  chain->SetBranchStatus("srtrack.begPlane", 1);
  chain->SetBranchStatus("srtrack.endPlane*", 1);
  chain->SetBranchStatus("srtrack.planes", 1);
  chain->SetBranchStatus("srtrack.trklikePlanes", 1);
  chain->SetBranchStatus("srtrack.phCCGeV",  1);
  chain->SetBranchStatus("srtrack.passedFit",  1);
  chain->SetBranchStatus("srtrack.deltaUVVtx",  1);


  chain->SetBranchStatus("ann*", 0);
  chain->SetBranchStatus("ann.pid_30inp", 1);
  chain->SetBranchStatus("shi*", 0);
  chain->SetBranchStatus("shi.epi0", 1);
  chain->SetBranchStatus("shi.emenergy", 1);
  chain->SetBranchStatus("bmon*", 0);
  chain->SetBranchStatus("bmon.goodBeamMon", 1);
  chain->SetBranchStatus("bmon.dt_stnd", 1);
  chain->SetBranchStatus("subshowervars*", 0);
  chain->SetBranchStatus("subshowervars.pid*", 1);
  chain->SetBranchStatus("fluxinfo*", 0);
  chain->SetBranchStatus("fluxinfo.tptype", 1);

  chain->SetBranchStatus("mctrue*", 0);
  chain->SetBranchStatus("mctrue.interactionType", 1);
  chain->SetBranchStatus("mctrue.nu*", 1);
  chain->SetBranchStatus("mctrue.fNue*", 1);
  chain->SetBranchStatus("mctrue.leptonMom*", 1);
  chain->SetBranchStatus("mctrue.emShower*", 1);
  chain->SetBranchStatus("mctrue.trueVisible*", 1);
  chain->SetBranchStatus("mctrue.hadronic*", 1);
  chain->SetBranchStatus("mctrue.atomic*", 1);
  chain->SetBranchStatus("mctrue.resonance*", 1);
  chain->SetBranchStatus("mctrue.q2", 1);
  chain->SetBranchStatus("mctrue.w2", 1);
  chain->SetBranchStatus("mctrue.bj*", 1);
  chain->SetBranchStatus("mctrue.target*", 1);
  chain->SetBranchStatus("mctrue.nonOscNuFlavor*", 1);
  chain->SetBranchStatus("mctrue.initialState", 1);
                                                                                
  chain->SetBranchStatus("srevent*", 0);
  chain->SetBranchStatus("srevent.vtx*", 1);
  chain->SetBranchStatus("srevent.vtxM*", 0);
  chain->SetBranchStatus("srevent.ph*", 1);
  chain->SetBranchStatus("srevent.triggerTime", 1);
  chain->SetBranchStatus("srevent.spillType", 1);
  chain->SetBranchStatus("srevent.coil*", 1);
  chain->SetBranchStatus("srevent.liTime", 1);
  chain->SetBranchStatus("srevent.eventTimeMin", 1);
  chain->SetBranchStatus("srevent.rcBoundary", 1);
  chain->SetBranchStatus("srevent.daveFD*", 1);
  chain->SetBranchStatus("fHeader*", 0);
  chain->SetBranchStatus("fHeader.fVld*", 1);
  chain->SetBranchStatus("fHeader.fRelease", 1);
  chain->SetBranchStatus("srevent.tracks",1);
  chain->SetBranchStatus("anainfo.abCCPID",1);

}

Double_t NueExtrapHelper::GetNueEnergy(Selection::Selection_t sel)
{
  return this->GetNueEnergy(fRecord, sel);
}

Double_t NueExtrapHelper::GetNueEnergy(NueRecord * nRecord, Selection::Selection_t sel)
{
  double ntupleEnergy = nRecord->srevent.phNueGeV;
  if(sel != Selection::kCC){
    return ntupleEnergy;
  }
  else{
   Double_t recoE = 0.0;
   if(nRecord->srtrack.phCCGeV > 0) recoE += nRecord->srtrack.phCCGeV;
   if(nRecord->srshower.phCCGeV > 0) recoE += nRecord->srshower.phCCGeV;
   return recoE;
  }

  // Killed by above at present
  Double_t fracDiff = 0;
  Double_t corEnergy = ntupleEnergy;
  if(fracDiff+1!=0) corEnergy /= (fracDiff+1);
  return corEnergy;
}
