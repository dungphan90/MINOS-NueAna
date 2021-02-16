#include <iostream>
#include "TMath.h"
#include "NueAna/NueRWHelpers.h"
#include "AnalysisNtuples/ANtpTruthInfoBeam.h"
#include "MCNtuple/NtpMCTruth.h"
#include "NueAna/OscProb.h"
#include "NueAna/ANtpTruthInfoBeamNue.h"

float NueRWHelpers::Oscillate(NtpMCTruth *mcth,
			      float L, float dm2, float theta23, float UE32)
{
    return NueRWHelpers::Oscillate(mcth->inu, mcth->inunoosc, mcth->p4neu[3],
                                   L, dm2, theta23, UE32);
}  

float NueRWHelpers::Oscillate(ANtpTruthInfoBeam *ib,
			      float L, float dm2, float theta23, float UE32)
{
  return NueRWHelpers::Oscillate(ib->nuFlavor, ib->nonOscNuFlavor,ib->nuEnergy,
                                   L, dm2, theta23, UE32);
}

float NueRWHelpers::Oscillate(ANtpTruthInfoBeamNue *ib)
{
   return NueRWHelpers::Oscillate(ib->nuFlavor, ib->nonOscNuFlavor,ib->nuEnergy,
                                  ib->Baseline, ib->DeltamSquared23, 
				  ib->Theta23, ib->Ue3Squared);
}



float NueRWHelpers::Oscillate(int nuFlavor, int nonOscNuFlavor, float Energy, 
			      float L, float dm2, float theta23, float UE32)
{
   float oscterm = TMath::Sin(1.27*dm2*L/Energy);      

   float pmt=pow((1-UE32)*oscterm*TMath::Sin(2*theta23),2);
   float pme=pow(TMath::Sin(theta23),2)*4.*UE32*(1-UE32)*pow(oscterm,2);
   float pmm=1.-pmt-pme;
   
   float pet=4*(1-UE32)*UE32*pow(TMath::Cos(theta23)*oscterm,2);
   float pem=pow(TMath::Sin(theta23),2)*4.*UE32*(1-UE32)*pow(oscterm,2);
   float pee=1.-pet-pem;


   if(abs(nonOscNuFlavor)==14){
      if(abs(nuFlavor)==12){
	 return pme;
      }
      else if(abs(nuFlavor)==14){
	 return pmm;
      }
      else if(abs(nuFlavor)==16){
	 return pmt;
      }
   }
   else if(abs(nonOscNuFlavor)==12){
      if(abs(nuFlavor)==12){
	 return pee;
      }
      else if(abs(nuFlavor)==14){
	 return pem;
      }
      else if(abs(nuFlavor)==16){
	 return pet;
      }
   }
   else{
     std::cout<<"I don't know what to do with "<<nonOscNuFlavor
	 <<" "<<nuFlavor<<" "<<pee<<std::endl;
   }
   return 0.;
}

float NueRWHelpers::OscillateMatter(int nuFlavor, int nonOscNuFlavor, 
				    float Energy,
				    float L, float dm2, float th23, float UE32,
				    float delta,int hierarchy)
{

  Double_t x[1] = {};
  x[0] = Energy;
  Double_t th12 = 0.554;              //sinsq2theta_12=0.8
  Double_t ss2th13 = 4*UE32*(1-UE32); //sinsq2theta_13
  Double_t dm2_12 = 8.2e-5; //best fit SNO
  Double_t dm2_23 = dm2;

  Double_t par[9] = {0};
  par[0] = L;
  par[1] = th23;
  par[2] = th12;
  par[3] = TMath::ASin(TMath::Sqrt(ss2th13))/2.;
  par[4] = hierarchy*dm2_23;
  par[5] = hierarchy*dm2_12;
  par[6] = 2.65; //standard rock density
  par[7] = delta;
  par[8] = 1;
  if(nonOscNuFlavor < 0) par[8] = -1;
  
  if(nonOscNuFlavor==14) {
    if(nuFlavor==12) return ElecAppear(x,par);
    else if(nuFlavor==14) return MuSurvive(x,par) - ElecAppear(x,par);
    else if(nuFlavor==16) return MuDisappear(x,par);
  }
  else if(nonOscNuFlavor==12) {
    if(nuFlavor==12) return 1 - ElecAppear(x,par) - 
      ( ss2th13 * TMath::Power(TMath::Cos(th23) * 
			       TMath::Sin(1.27*dm2*L/x[0]),2) );
    else if(nuFlavor==14) return ElecAppear(x,par);
    else if(nuFlavor==16) return ( ss2th13 *
				   TMath::Power(TMath::Cos(th23) * 
						TMath::Sin(1.27*dm2*L/x[0]),2) );
  }
  return 0;
}

float NueRWHelpers::OscillateMatter(NtpMCTruth *mcth,
				    float L, float dm2, float theta23, 
				    float UE32,float delta,int hierarchy)
{
  return NueRWHelpers::OscillateMatter(mcth->inu, mcth->inunoosc, 
				       mcth->p4neu[3],
				       L, dm2, theta23, UE32,
				       delta,hierarchy);
}

float NueRWHelpers::OscillateMatter(ANtpTruthInfoBeam *ib,
				    float L, float dm2, float theta23, 
				    float UE32,float delta, int hierarchy)
{
  return NueRWHelpers::OscillateMatter(ib->nuFlavor, ib->nonOscNuFlavor,
				       ib->nuEnergy,
				       L, dm2, theta23, UE32,
				       delta,hierarchy);
}

