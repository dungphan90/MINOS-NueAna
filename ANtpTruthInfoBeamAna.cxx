#include "CandNtupleSR/NtpSREvent.h"
#include "StandardNtuple/NtpStRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include "MCNtuple/NtpMCStdHep.h"
#include "TruthHelperNtuple/NtpTHEvent.h"
#include "TruthHelperNtuple/NtpTHTrack.h"
#include "TruthHelperNtuple/NtpTHShower.h"
#include "TruthHelperNtuple/NtpTHStrip.h"
#include "MessageService/MsgService.h"
#include "NueAna/ANtpTruthInfoBeamAna.h"
#include "AnalysisNtuples/Module/ANtpRecoNtpManipulator.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "NueAna/NueAnaTools/NueConvention.h"
#include "TVector3.h"
#include "OscProb/OscCalc.h"

CVSID("$Id: ANtpTruthInfoBeamAna.cxx,v 1.32 2014/02/18 03:44:20 rhatcher Exp $");

ANtpTruthInfoBeamAna::ANtpTruthInfoBeamAna(ANtpTruthInfoBeamNue &antib):
   fANtpTruthInfoBeam(antib)
{}

ANtpTruthInfoBeamAna::~ANtpTruthInfoBeamAna()
{}

void ANtpTruthInfoBeamAna::Analyze(int event, RecRecordImp<RecCandHeader> *srobj)
{
  if(srobj==0){
    return;
  }
  NtpStRecord *st = dynamic_cast<NtpStRecord *>(srobj);
  Analyze(event,st);
}

void ANtpTruthInfoBeamAna::Analyze(int evtn, NtpStRecord *srobj)
{
  if(srobj==0){
    return;
  }

    // This is mostly copied from AnalysisNtuples/Module/CondensedNtpModule.cxx
    NtpMCTruth *mctruth = 0;
    NtpMCStdHep *mcstdhep = 0;
    NtpTHEvent *thevent = 0;
    //NtpTHStrip *thstrip = 0;
    //NtpTHTrack *thtrack = 0;
    //NtpTHShower *thshower = 0;

    //instansiate a NtpHelper object to help you get the info you want
    ANtpRecoNtpManipulator ntpManipulator(srobj);
    fInfoFiller = new ANtpInfoObjectFillerBeam();
    fInfoFiller->SetStripArray(ntpManipulator.GetStripArray());
  
    thevent = ntpManipulator.GetMCManipulator()->GetNtpTHEvent(evtn);
    if(thevent){
      mctruth = ntpManipulator.GetMCManipulator()->GetNtpMCTruth(thevent->neumc);
      mcstdhep = ntpManipulator.GetMCManipulator()->GetNtpMCStdHep(thevent->neustdhep);
    }
				  
    // no strip info yet
    //if(ntpStrip) ntpTHStrip = ntpManipulator.GetNtpTHStrip(ntpStrip->index);
    //if(ntpTHStrip){
    //    ntpMCTruth = ntpManipulator.GetNtpMCTruth(ntpTHStrip->neumc);   
    //}
    
    if(mctruth){
      // fill the ANtpTruthInfo part
      fInfoFiller->FillMCTruthInformation(mctruth, &fANtpTruthInfoBeam);
      // fill the ANtpTruthInfoBeam part
      fInfoFiller->FillBeamMCTruthInformation(mctruth, ntpManipulator.GetStdHepArray(),&fANtpTruthInfoBeam);

      Int_t inu = fANtpTruthInfoBeam.nuFlavor;
      Int_t inunoosc = fANtpTruthInfoBeam.nonOscNuFlavor;
      Int_t iaction = fANtpTruthInfoBeam.interactionType;          

      fANtpTruthInfoBeam.fNueClass = GetNueClass(inu, inunoosc, iaction);
      fANtpTruthInfoBeam.fNueWeight = GetNueWeight(inu, inunoosc);
      if(srobj->GetHeader().GetVldContext().GetDetector() == Detector::kFar)
        GetOscProb();
      else
        fANtpTruthInfoBeam.fOscProb = 1.0;

      fANtpTruthInfoBeam.DirCosNeu = TrueLepDCosNeu(mctruth);
      fANtpTruthInfoBeam.DirCosZ_pan = TrueLepDCosZ(mctruth);

      fANtpTruthInfoBeam.istruckq = mctruth->istruckq;
      fANtpTruthInfoBeam.iflags = mctruth->iflags;
      fANtpTruthInfoBeam.sigmadiff = mctruth->sigmadiff;
      fANtpTruthInfoBeam.itg = mctruth->itg;

      fANtpTruthInfoBeam.eventCompleteness = thevent->completeall;
    }

    if(fInfoFiller){
      delete fInfoFiller;
      fInfoFiller=0;
    }                               	
}

void ANtpTruthInfoBeamAna::Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord * mcobj, NtpTHRecord * thobj)
{

  if(srobj==0){
    return;
  }
  if(mcobj==0){
    return;
  }
  if(thobj==0){
    return;
  }

    fInfoFiller = new ANtpInfoObjectFillerBeam();
    // This is mostly copied from AnalysisNtuples/Module/CondensedNtpModule.cxx
    NtpMCTruth *mctruth = 0;
    NtpMCStdHep *mcstdhep = 0;
    NtpTHEvent *thevent = 0;
    //NtpTHStrip *thstrip = 0;
    //NtpTHTrack *thtrack = 0;
    //NtpTHShower *thshower = 0;

    // Only fill if there is MC or Truth information
    if (mcobj && thobj){

        //instansiate a NtpHelper object to help you get the info you want
        ANtpRecoNtpManipulator ntpManipulator(srobj,mcobj,thobj);
        
        thevent = ntpManipulator.GetMCManipulator()->GetNtpTHEvent(evtn);
        if(thevent){
         mctruth = ntpManipulator.GetMCManipulator()->GetNtpMCTruth(thevent->neumc);
         mcstdhep = ntpManipulator.GetMCManipulator()->GetNtpMCStdHep(thevent->neustdhep);
	}
	
        // no strip info yet
        //if(ntpStrip) ntpTHStrip = ntpManipulator.GetNtpTHStrip(ntpStrip->index);
        //if(ntpTHStrip){
        //    ntpMCTruth = ntpManipulator.GetNtpMCTruth(ntpTHStrip->neumc);   
        //}

        if(mctruth){
	  // fill the ANtpTruthInfo part
	  fInfoFiller->FillMCTruthInformation(mctruth, &fANtpTruthInfoBeam);
	  // fill the ANtpTruthInfoBeam part
	  fInfoFiller->FillBeamMCTruthInformation(mctruth, ntpManipulator.GetStdHepArray(),&fANtpTruthInfoBeam);

          Int_t inu = fANtpTruthInfoBeam.nuFlavor;
          Int_t inunoosc = fANtpTruthInfoBeam.nonOscNuFlavor;
          Int_t iaction = fANtpTruthInfoBeam.interactionType;          

          fANtpTruthInfoBeam.fNueClass = GetNueClass(inu, inunoosc, iaction);
          fANtpTruthInfoBeam.fNueWeight = GetNueWeight(inu, inunoosc);
          if(srobj->GetHeader().GetVldContext().GetDetector() == Detector::kFar)
            GetOscProb();
          else
	    fANtpTruthInfoBeam.fOscProb = 1.0;
          
          fANtpTruthInfoBeam.DirCosNeu = TrueLepDCosNeu(mctruth);
          fANtpTruthInfoBeam.DirCosZ_pan = TrueLepDCosZ(mctruth);
        }                               
    }	

    if(fInfoFiller){
      delete fInfoFiller;
      fInfoFiller=0;
    }
}

Int_t ANtpTruthInfoBeamAna::GetNueClass(Int_t inu, Int_t inunoosc, Int_t iaction)
{
  Int_t cType=ANtpDefVal::kInt;
  cType = ClassType::DetermineClassType(inu, inunoosc, iaction);

  if(ANtpDefVal::IsDefault(cType)){
       MSG("ANtpTruthInfoBeamAna",Msg::kWarning)<< " GetNueClass called " 
                       << "eith invalid paramters - no class assigned" 
                       << endl;            
  }

  return cType;
}


Float_t ANtpTruthInfoBeamAna::GetNueWeight(Int_t inu, Int_t inunoosc)
{
                                                                          
  Float_t weight = ANtpDefVal::kFloat;

  if(ReleaseType::IsDaikon(release)) weight = 1.0;

  if(ReleaseType::IsCarrot(release)){
    if(abs(inu) == 12 || abs(inu) == 14 || abs(inu) == 16){
       if(abs(inunoosc) == 12 || abs(inunoosc) == 14 || abs(inunoosc) == 16){
          weight = 1.0;
          if(abs(inu) ==12 && abs(inunoosc) == 12)    weight = 0.5;
       }
    }
  }

  if(ANtpDefVal::IsDefault(weight)){
    // MSG("ANtpTruthInfoBeamAna",Msg::kWarning)<< " GetNueWeight called " 
    //						 << "with invalid paramters "
    //					 <<inu<<" "<<inunoosc
    //					 <<" - 0 weight assigned" 
    //					 << endl;            
        return 0.0;
  }
                          
  return weight;  
}

Float_t ANtpTruthInfoBeamAna::GetOscProb()
{
    fANtpTruthInfoBeam.Baseline = 735;
    fANtpTruthInfoBeam.DeltamSquared23 = 0.0024;
    double th23 = fANtpTruthInfoBeam.Theta23 = TMath::Pi()/4.0;
    fANtpTruthInfoBeam.Ue3Squared = 0.025;
    
// CC Results for MDC 
//    dm2 = 0.002175;  // +- 0.15e-3
    double th13 = TMath::ASin(TMath::Sqrt(fANtpTruthInfoBeam.Ue3Squared)); 

   fANtpTruthInfoBeam.fOscProb = NueConvention::Oscillate(&fANtpTruthInfoBeam);
  
   Double_t dm2_12 = 8.7e-5;
   double th12 = 0.816;
   Double_t dm2_23 = fANtpTruthInfoBeam.DeltamSquared23;

   static OscCalc fOscCalc;   
                                                                             
   Double_t par[9] = {0};
   par[OscPar::kL] = fANtpTruthInfoBeam.Baseline;
   par[OscPar::kTh23] = fANtpTruthInfoBeam.Theta23 = th23;
   par[OscPar::kTh12] = fANtpTruthInfoBeam.Theta12 = th12;
   par[OscPar::kTh13] = fANtpTruthInfoBeam.Theta13 = th13;
   par[OscPar::kDeltaM23] = dm2_23;
   par[OscPar::kDeltaM12] = fANtpTruthInfoBeam.DeltamSquared12 = dm2_12;
   par[OscPar::kDensity] = fANtpTruthInfoBeam.Density = 2.65; 
   par[OscPar::kDelta] = fANtpTruthInfoBeam.Delta = 0.0;
   par[OscPar::kNuAntiNu] = 1;
                                                                                
//  std::cout<<"About to call "<<dm2<<"  "<<ss13<<"  "<<delta<<std::endl;
   fOscCalc.SetOscParam(par);
 
   int nuFlavor = fANtpTruthInfoBeam.nuFlavor;
   int nonOscNuFlavor = fANtpTruthInfoBeam.nonOscNuFlavor;

   if(nonOscNuFlavor > 0) fOscCalc.SetOscParam(OscPar::kNuAntiNu, 1);
   if(nonOscNuFlavor < 0) fOscCalc.SetOscParam(OscPar::kNuAntiNu, -1);

   double Energy =  fANtpTruthInfoBeam.nuEnergy;

                                                                                
   fANtpTruthInfoBeam.fOscProbMatterNormal = fOscCalc.Oscillate(nuFlavor, nonOscNuFlavor, Energy);
    
   fOscCalc.SetOscParam(OscPar::kDeltaM23, -dm2_23);
   fANtpTruthInfoBeam.fOscProbMatterInverted = fOscCalc.Oscillate(nuFlavor, nonOscNuFlavor, Energy);


  return 1; 
}


Float_t ANtpTruthInfoBeamAna::TrueLepDCosNeu(NtpMCTruth *ntpTruth)
{ //ds_mu/ds_nu
  TVector3 *nuvec = new TVector3(ntpTruth->p4neu[0],
                                 ntpTruth->p4neu[1],
                                 ntpTruth->p4neu[2]);

  Float_t p4_0 = 0;
  Float_t p4_1 = 0;
  Float_t p4_2 = 0;

  Get3Momenta(ntpTruth, p4_0, p4_1, p4_2);

  TVector3 *muvec = new TVector3(p4_0, p4_1, p4_2);
  Float_t MuAng = nuvec->Angle(*muvec); //angle in rads
  delete nuvec;
  delete muvec;
  return TMath::Cos(MuAng);
}

Float_t ANtpTruthInfoBeamAna::TrueLepDCosZ(NtpMCTruth *ntpTruth){ //dz/ds

  Float_t p4_0 = 0;
  Float_t p4_1 = 0;
  Float_t p4_2 = 0;

  Float_t mom = Get3Momenta(ntpTruth, p4_0, p4_1, p4_2);

  if(mom==0) return 0;
  return p4_2/mom;
}

Float_t ANtpTruthInfoBeamAna::Get3Momenta(NtpMCTruth *ntpTruth,
                           Float_t &p4_0, Float_t &p4_1, Float_t &p4_2)
{ 
  p4_0 = 0;
  p4_1 = 0;
  p4_2 = 0;
                                                                                
  if(TMath::Abs(ntpTruth->p4mu1[3]) > 0){
    p4_0 = ntpTruth->p4mu1[0];
    p4_1 = ntpTruth->p4mu1[1];
    p4_2 = ntpTruth->p4mu1[2];
  }
  else if(TMath::Abs(ntpTruth->p4el1[3])>0){
    p4_0 = ntpTruth->p4el1[0];
    p4_1 = ntpTruth->p4el1[1];
    p4_2 = ntpTruth->p4el1[2];
  }
  else if(TMath::Abs(ntpTruth->p4tau[3])>0){
    p4_0 = ntpTruth->p4tau[0];
    p4_1 = ntpTruth->p4tau[1];
    p4_2 = ntpTruth->p4tau[2];
  }
  
  return (p4_0*p4_0 + p4_1*p4_1 + p4_2*p4_2);

}
