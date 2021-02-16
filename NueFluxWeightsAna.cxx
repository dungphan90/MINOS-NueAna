#include "StandardNtuple/NtpStRecord.h"
#include "NueAna/NueAnaTools/SntpHelpers.h"
#include "MCNtuple/NtpMCTruth.h"
#include "MessageService/MsgService.h"
#include "MCReweight/Zbeam.h"
#include "MCReweight/Zfluk.h"
#include "MCReweight/Kfluk.h"
#include "NueAna/NueFluxWeightsAna.h"

CVSID("$Id: NueFluxWeightsAna.cxx,v 1.14 2009/06/24 22:43:52 vahle Exp $");

NueFluxWeightsAna::NueFluxWeightsAna(NueFluxWeights &nuefw):
  fNueFluxWeight(nuefw),
//  beam(2),
  det(1)
{
  //   cfg = "PiMinus_CedarDaikon";
   cfg = "DetXs";
}

NueFluxWeightsAna::~NueFluxWeightsAna()
{
  MSG("NueFluxWeightsAns",Msg::kDebug)<<"in NueFluxWeightsAna destructor"<<endl;
}

void NueFluxWeightsAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
  NtpStRecord *st=dynamic_cast<NtpStRecord *>(srobj);
  if(st==0){
    MSG("NueFluxWeightsAna",Msg::kError)<<"Trying to do flux reweighting on an event"
					<<" that comes from a non NtpStRecord"<<endl
					<<"That's not a good idea"<<endl;
    return;
  }

  int thn = SntpHelpers::GetEvent2MCIndex(evtn,st);
  NtpMCTruth *mcrec = SntpHelpers::GetMCTruth(thn,st);
  if(mcrec==0){
    return;
  }
  MSG("NueFluxWeightsAns",Msg::kDebug)<<"in NueFluxWeightsAna::Analyze"<<endl;

  if(fi==0){
    MSG("NueFluxWeightsAna",Msg::kWarning)<<"No FluxInfo object set, "
					  <<"using the one from NtpMCTruth"
					  <<" Once the flux files are fixed, "
					  <<"and we no longer need the MuPi trees, "
					  <<"we can comment out this comment"<<endl;
    fi=&mcrec->flux;
  }

  MSG("NueFluxWeightsAna",Msg::kDebug)<<"starting flux reweight"<<endl;

  double pt = sqrt(fi->tpx*fi->tpx+fi->tpy*fi->tpy);
  double pz = 1.*fi->tpz;
  int tptype = fi->tptype;
  int inu = mcrec->inu; 
  int cc_nc = mcrec->iaction;

  float true_enu = mcrec->p4neu[3];

  zbeam = BeamType::ToZarko(beam);
  MSG("NueFluxWeightsAna",Msg::kDebug)<<"true_enu "<<true_enu<<" beam (Z) "<<zbeam<<" det "<<det<<endl;
  double kflukw = 1.;//kfluk->GetWeight(fi);
  fNueFluxWeight.kflukweight = kflukw;

  double newTrackE, newShwE;
  double bweight =  skzpCalc->GetBeamWeight(det,zbeam,tptype,pt,pz,true_enu,inu);
  double dweight =  skzpCalc->GetDetWeight(cc_nc,true_enu,inu,0,0,newTrackE,newShwE);
  
  fNueFluxWeight.totbeamweight = bweight;
  fNueFluxWeight.totskzpweight = dweight*bweight;
  fNueFluxWeight.detectorWeight = dweight;

  fNueFluxWeight.skzpTrkEnergy = 0.0;  //newTrackE;
  fNueFluxWeight.skzpShwEnergy = 0.0; //newShwE;
  fNueFluxWeight.skzpConfig = cfg;

  for(int rpit=SKZPWeightCalculator::kRunI;rpit<SKZPWeightCalculator::kEndOfList;rpit++){
    SKZPWeightCalculator::RunPeriod_t rp = skzpCalc->RunPeriodFromInt(rpit);
    double rpw = skzpCalc->GetBeamWeight(det,zbeam,tptype,pt,pz,true_enu,inu,rp);
    fNueFluxWeight.RPtotbeamweight.push_back(rpw);
  }


  MSG("NueFluxWeightsAna",Msg::kDebug)<<"alldone with flux reweight "<<fNueFluxWeight.totbeamweight<<endl;
}

