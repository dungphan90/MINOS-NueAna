#include "MessageService/MsgService.h"
#include "TVector3.h"
#include "StandardNtuple/NtpStRecord.h"
#include "MCNtuple/NtpMCRecord.h"
#include "CandNtupleSR/NtpSRRecord.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NueRecordAna.h"

CVSID("$Id: NueRecordAna.cxx,v 1.42 2009/06/26 19:15:00 xbhuang Exp $");

NueRecordAna::NueRecordAna(NueRecord &nr):
  sfa(nr.shwfit),
  hca(nr.hitcalc),
  aca(nr.angcluster),
  acf(nr.angclusterfit),
  msta(nr.mstvars),
  fva(nr.fracvars),
  sva(nr.subshowervars),
  hhv(nr.highhitvars),
  shrva(nr.shieldrejvars),
  anna(nr),
  anaia(nr.anainfo),
  aneia(nr.srevent),
  ansia(nr.srshower),
  antia(nr.srtrack),
  antiba(nr.mctrue),
//  bmona(nr.bmon),
  mda(nr,nr.mdadiscrim),
  tpa(nr,nr.treepid),
  fiana(nr.fluxinfo),
  nuefwa(nr.fluxweights),
  nuexsa(nr.xsecweights),
  shia(nr.shi),
  mria(nr.mri),
  cdia(nr.cdi),
  tva(nr.timingvars),
  bta(nr, nr.dtree),
  equala(nr)
{
  active3DHits = false;
}

void NueRecordAna::SetEventEnergyArray(float* array, float *ar2)
{
  evtstp0mip = array; evtstp1mip = ar2;
  sfa.SetEventEnergyArray(array, ar2);
  hca.SetEventEnergyArray(array, ar2);
  aca.SetEventEnergyArray(array, ar2);
  acf.SetEventEnergyArray(array, ar2);
  msta.SetEventEnergyArray(array, ar2);
  fva.SetEventEnergyArray(array, ar2);
  sva.SetEventEnergyArray(array, ar2);
  hhv.SetEventEnergyArray(array, ar2);
  shrva.SetEventEnergyArray(array, ar2);
//  anna.SetEventEnergyArray(array, ar2);
  anaia.SetEventEnergyArray(array, ar2);
  aneia.SetEventEnergyArray(array, ar2);
  ansia.SetEventEnergyArray(array, ar2);
  antia.SetEventEnergyArray(array, ar2);
  antiba.SetEventEnergyArray(array, ar2);
//  mda.SetEventEnergyArray(array, ar2);
//  tpa.SetEventEnergyArray(array, ar2);
  fiana.SetEventEnergyArray(array, ar2);
  nuefwa.SetEventEnergyArray(array, ar2);
  nuexsa.SetEventEnergyArray(array, ar2);
  shia.SetEventEnergyArray(array, ar2);
  mria.SetEventEnergyArray(array, ar2);
  cdia.SetEventEnergyArray(array, ar2);
  tva.SetEventEnergyArray(array, ar2);
  
}

void NueRecordAna::SetRelease(int rel)
{
  release = rel;
  sfa.SetRelease(rel);
  hca.SetRelease(rel);
  aca.SetRelease(rel);
  acf.SetRelease(rel);
  msta.SetRelease(rel);
  fva.SetRelease(rel);
  sva.SetRelease(rel);
  hhv.SetRelease(rel);
  shrva.SetRelease(rel);
//  anna.SetRecoRelease(rel);
  anaia.SetRelease(rel);
  aneia.SetRelease(rel);
  ansia.SetRelease(rel);
  antia.SetRelease(rel);
  antiba.SetRelease(rel);
//  mda.SetRecoRelease(rel);
//  tpa.SetRecoRelease(rel);
  fiana.SetRelease(rel);
  nuefwa.SetRelease(rel);
  nuexsa.SetRelease(rel);
  shia.SetRelease(rel);
  mria.SetRelease(rel);
  cdia.SetRelease(rel);
  tva.SetRelease(rel);
//  bta.SetRelease(rel);
  equala.SetRelease(rel);
}



NueRecordAna::~NueRecordAna()
{}

//void NueRecordAna::Analyze(int evtn, NtpSRRecord *srobj, NtpMCRecord *mcobj, NtpTHRecord *thobj)
void NueRecordAna::Analyze(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
  
  MSG("NueRecordAna",Msg::kDebug)<<"On Snarl "<<srobj->GetHeader().GetSnarl()
				 <<" event "<<evtn<<endl;
  
   DeqFloat_t x,y,z,e;
   Int_t primShow;
   DeqDeqInt_t clusterMap;
   TVector3 primDir;
   
   sfa.Analyze(evtn,srobj);
   MSG("NueRecordAna",Msg::kDebug)<<"done sfa"<<endl;
//   hca.Analyze(evtn,srobj);

   if(active3DHits){
     hca.Analyze(evtn, srobj);
     MSG("NueRecordAna",Msg::kDebug)<<"done hca"<<endl;
     hca.Get3DHit(x,y,z,e);    
     aca.Set3DHit(x,y,z,e);
     aca.Analyze(evtn,srobj);    // Must be preceded by hca.Get3DHit
                               //                     aca.Set3DHit
     MSG("NueRecordAna",Msg::kDebug)<<"done aca"<<endl;

     aca.GetAngCluster(primShow,clusterMap,primDir);
     acf.Set3DHit(x,y,z,e);
     acf.SetAngCluster(primShow,clusterMap,primDir);
     acf.Analyze(evtn,srobj);    // Must be preceded by:
                               // aca.GetAngCluster
                               // acf.Set3DHit
                               // acf.SetAngCluster
     MSG("NueRecordAna",Msg::kDebug)<<"done acf"<<endl;
   } 
   msta.Analyze(evtn,srobj);
   MSG("NueRecordAna",Msg::kDebug)<<"done msta"<<endl;

   fva.Analyze(evtn,srobj);
   MSG("NueRecordAna",Msg::kDebug)<<"done fva"<<endl;

   sva.Analyze(evtn,srobj);
   MSG("NueRecordAna",Msg::kDebug)<<"done sva"<<endl;


   hhv.Analyze(evtn,srobj);
   MSG("NueRecordAna",Msg::kDebug)<<"done hhv"<<endl;


   shrva.Analyze(evtn,srobj);
   MSG("NueRecordAna",Msg::kDebug)<<"done shrva"<<endl;

//   Is now being handled by a separate job module   
//   bmona.Analyze(evtn,srobj);
//   MSG("NueRecordAna",Msg::kDebug)<<"done bmona"<<endl;

//   anaia.Set3DHit(x,y,z,e);   //Must be preceded by hca.Get3Dhit   
//   anaia.Analyze(evtn,srobj); //Must be preceded by Set3DHit()   
   MSG("NueRecordAna",Msg::kDebug)<<"done anaia"<<endl;

   anna.Analyze();
   MSG("NueRecordAna",Msg::kDebug)<<"done anna"<<endl;

   mda.Analyze();
   MSG("NueRecordAna",Msg::kDebug)<<"done mda"<<endl;

   tpa.Analyze();
   MSG("NueRecordAna",Msg::kDebug)<<"done tpa"<<endl;
   
   tva.Analyze(evtn,srobj);
   MSG("NueRecordAna",Msg::kDebug)<<"done tva"<<endl;

   bta.Analyze();
   equala.Analyze(evtn,srobj);
   return;
}


void NueRecordAna::Analyze(int evtn, NtpMRRecord *mr, NtpStRecord *oldst)
{
  MSG("NueRecordAna",Msg::kDebug)<<"On Snarl "<<mr->GetHeader().GetSnarl()
				 <<" event "<<evtn<<endl;
  mria.Analyze(evtn,mr,oldst);
  MSG("NueRecordAna",Msg::kDebug)<<"done mri"<<endl;
  return;
}

void NueRecordAna::Analyze(UberRecord *ur)
{
  MSG("NueRecordAna",Msg::kDebug)<<"Time of current event "
				 <<ur->GetHeader().GetStartTime()
				 <<endl;
  cdia.Analyze(ur);
  MSG("NueRecordAna",Msg::kDebug)<<"done cdi"<<endl;
  return;
}

void NueRecordAna::FillTrue(int evtn, NtpSRRecord *srobj, NtpMCRecord *mcobj, NtpTHRecord *thobj)
{
  if(mcobj==0){
    return;
  }

  VldContext vc;
  if(srobj!=0){
    vc=srobj->GetHeader().GetVldContext();
  }
  else{
    vc=mcobj->GetHeader().GetVldContext();
  }

  if(vc.GetSimFlag()==SimFlag::kData){
    return;
  }
 
    antiba.Analyze(evtn,srobj,mcobj,thobj);
    return;
}

void NueRecordAna::FillTrue(int evtn, NtpStRecord *srobj)
{
  if(srobj==0){
    return;
  }

  VldContext vc=srobj->GetHeader().GetVldContext();
  if(vc.GetSimFlag()==SimFlag::kData){
    return;
  }

  antiba.Analyze(evtn,srobj);
  fiana.Analyze(evtn,srobj);
  shia.Analyze(evtn,srobj);
  Detector::Detector_t d=vc.GetDetector();
  if(d==Detector::kNear){
    nuefwa.SetDetector(1);
  }
  else if(d==Detector::kFar){
    nuefwa.SetDetector(2);
  }
  else{
    nuefwa.SetDetector(0);
  }
  nuefwa.SetFluxInfo(fiana.GetMCFluxInfo());
  nuefwa.Analyze(evtn,srobj);
  nuexsa.Analyze(evtn,srobj);
  return;
}

void NueRecordAna::FillReco(int evtn, RecRecordImp<RecCandHeader> *srobj)
{
  if(srobj==0){
    return;
  }
 
   aneia.Analyze(evtn,srobj);
   ansia.Analyze(evtn,srobj);
   antia.Analyze(evtn,srobj);
   anaia.Analyze(evtn,srobj);

   return;
}

void NueRecordAna::SetBeamType(int bType)
{
  return SetBeamType((BeamType::BeamType_t) bType);
}
                                                                                
void NueRecordAna::SetBeamType(BeamType::BeamType_t bType)
{
  sfa.SetBeamType(bType);
  hca.SetBeamType(bType);
  aca.SetBeamType(bType);
  acf.SetBeamType(bType);
  msta.SetBeamType(bType);
  fva.SetBeamType(bType);
  hhv.SetBeamType(bType);
  sva.SetBeamType(bType);
  shrva.SetBeamType(bType);
//  anna.SetRecoRelease(rel);
  anaia.SetBeamType(bType);
  aneia.SetBeamType(bType);
  ansia.SetBeamType(bType);
  antia.SetBeamType(bType);
  antiba.SetBeamType(bType);
//  mda.SetRecoRelease(rel);
//  tpa.SetRecoRelease(rel);
  fiana.SetBeamType(bType);
  nuefwa.SetBeamType(bType);
  nuexsa.SetBeamType(bType);
  shia.SetBeamType(bType);
  mria.SetBeamType(bType);
  cdia.SetBeamType(bType);
  tva.SetBeamType(bType);

}

