#include <iostream>
#include <cmath>
#include "NueAna/Reweight/NueRW.h"


ClassImp(NueRW)

const char * NueRW::AsString(FileType_t t){
  switch(t){
  case kUnknown: return "Unknown"; break;
  case kBEAM: return "Beam"; break;
  case kNUE: return "Nue"; break;
  case kTAU: return "Tau"; break;
  case kAGG: return "Aggregate"; break;
  default: return "Unknown"; break;
  }
  return "Unknown";
}

NueRW::NueRW():
   fRun(0),
   fSubRun(0),
   fDet(Detector::kUnknown),
   fFileType(NueRW::kUnknown),
   ftgt(0),
   qel_ma(0.),
   res_ma(0.),
   coh_ma(0.),
   qel_fa0(0.),
   qel_eta(0.),
   res_omega(0.),
   res_z(0.),
   coh_r0(0.),
   coh_rei(0.),
   kno_a1(0.),
   kno_a2(0.),
   kno_a3(0.),
   kno_a4(0.),
   kno_b(0.),
   kno_r112(0.),
   kno_r122(0.),
   kno_r132(0.),
   kno_r142(0.),
   kno_r113(0.),
   kno_r123(0.),
   kno_r133(0.),
   kno_r143(0.),

   dm2(0.),
   ss2th(0.),
   UE32(0.0),

   randrow(-1),
   nfiles(0),
   nsnarls(0),
   nevents(0),
   neventswpid(0),
   nacc(0),
   nsig(0.),
   nbg(0.),
   nnueb(0.),
   nnumu(0.),
   nnutau(0.),
   nnc(0.),
   EBINS(20),
   EBINW(1.)

{

    nsigE=new float[EBINS];
    nbgE=new float[EBINS];
    nnuebE=new float[EBINS];
    nnumuE=new float[EBINS];
    nnutauE=new float[EBINS];
    nncE=new float[EBINS];
  

  for(int i=0;i<EBINS;i++){
    nsigE[i]=0.;
    nbgE[i]=0.;
    nnuebE[i]=0.;
    nnumuE[i]=0.;
    nnutauE[i]=0.;
    nncE[i]=0.;
  }
}

NueRW::~NueRW()
{}

void NueRW::Print(Option_t * /*option*/) const
{
   NueRW::Print(std::cout);
   return;
}

std::ostream &NueRW::Print(std::ostream &os) const
{
   os<<"Run: "<<fRun<<" SubRun "<<fSubRun<<std::endl;
   os<<"Det "<<Detector::AsString(fDet)
     <<" File "<<AsString(fFileType)
     <<" itgt "<<ftgt<<std::endl
     <<"Rand Row: "<<randrow<<std::endl
     <<"NSnarls: "<<nsnarls<<" NEvents "<<nevents<<" NEvents w/PID "<<neventswpid<<std::endl
     <<"NAccepted: "<<nacc<<std::endl
     <<"Signal: "<<nsig<<" Total BG: "<<nbg<<std::endl
     <<"N nue beam "<<nnueb<<" N numu "<<nnumu<<" N nutau "<<nnutau<<" N NC "<<nnc<<std::endl;
   return os;
}

void NueRW::Reset() 
{
   fRun=0;
   fSubRun=0;
   fDet=Detector::kUnknown;
   fFileType=kUnknown;
   ftgt=0;
   randrow=0;
   nsnarls=0,
   nevents=0,
   neventswpid=0,
   nacc=0;
   nsig=0.;
   nbg=0.;
   nnueb=0.;
   nnumu=0.;
   nnutau=0.;
   nnc=0.;

   for(int i=0;i<EBINS;i++){
    nsigE[i]=0.;
    nbgE[i]=0.;
    nnuebE[i]=0.;
    nnumuE[i]=0.;
    nnutauE[i]=0.;
    nncE[i]=0.;
  }
     

}

int NueRW::FindEBin(float E)
{
  if(E>1.e10){
    return EBINS-1;
  }
  int bin = (int)(fabs(E)/EBINW);
  if(bin>=EBINS){
    bin=EBINS-1;
  }
  return bin;

}
