////////////////////////////////////////////////////////////////////////
// $Id: NueConvention.cxx,v 1.28 2015/03/02 19:08:16 wingmc Exp $
//
// The NueConvention File will hopefully organize
//  frequently used nue standards such as
//    event class (was ClassType)
//    fiducial Volume
//    file indexing
//
// Author: Josh Boehm
// Created: Sept. 14, 2006
////////////////////////////////////////////////////////////////////////

#include "NueConvention.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

#include "TMath.h"
#include "AnalysisNtuples/ANtpTruthInfoBeam.h"
#include "NueAna/ANtpTruthInfoBeamNue.h"
#include "MCNtuple/NtpMCTruth.h"
#include "NueAna/NueRecord.h"
#include "NueAna/NueStandard.h"

#include "MessageService/MsgService.h"
// ******************************************************************
// Class Type ID Function
// ******************************************************************
CVSID("");

Int_t NueConvention::DetermineClassType(Int_t inu, Int_t inunoosc, Int_t iaction)
{
     int cType=ANtpDefVal::kInt;

     if(iaction ==0){ cType= ClassType::NC; // NC class
      }
     else if(iaction >=1){
       if(inu ==14 || inu==-14){ 
          cType=ClassType::numu;  // CC numu class
       }
       else
       if(inu==12 || inu==-12){
         if(inunoosc==14 || inunoosc==-14){
            cType=ClassType::nue; // CC osc nue class
         }
         else if(inunoosc==12 || inunoosc==-12){
            cType=ClassType::bnue; // CC beam nue class
         }
       }
       else if(inu==16 || inu==-16){
          cType= ClassType::nutau;  // CC nutau class
       }
    }

    return cType;
}


// ********************************************************************
// Reweighting code
// ********************************************************************


float NueConvention::Oscillate(NtpMCTruth *mcth,
                              float L, float dm2, float theta23, float UE32)
{
    return NueConvention::Oscillate(mcth->inu, mcth->inunoosc, mcth->p4neu[3],
                                   L, dm2, theta23, UE32);
}

float NueConvention::Oscillate(ANtpTruthInfoBeam *ib,
                              float L, float dm2, float theta23, float UE32)
{
  return NueConvention::Oscillate(ib->nuFlavor, ib->nonOscNuFlavor,ib->nuEnergy,
                                   L, dm2, theta23, UE32);
}

float NueConvention::Oscillate(ANtpTruthInfoBeamNue *ib)
{
   return NueConvention::Oscillate(ib->nuFlavor, ib->nonOscNuFlavor,ib->nuEnergy,                                  ib->Baseline, ib->DeltamSquared23,
                                  ib->Theta23, ib->Ue3Squared);
}


float NueConvention::Oscillate(int nuFlavor, int nonOscNuFlavor, float Energy,
                              float L, float dm2, float theta23, float UE32)
{
   float oscterm = TMath::Sin(1.269*dm2*L/Energy);

//   std::cout<<oscterm<<"  "<<pow(TMath::Sin(2*theta23),2)<<"  "
//       <<pow((1-UE32),2)<<std::endl;

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

// ********************************************************************
// Fiducial Volume Code
// ********************************************************************

int NueConvention::IsInsideNearFiducial_Nue_Extended(float x, float y, float z)
{
  Float_t SuperModule1Beg = 0.50;
  Float_t SuperModule1End = 6.50;

  Float_t radialInner = 0;
  Float_t radialOuter = 1;
  Float_t xCenter = 1.4885;
  Float_t yCenter = 0.1397;
  Bool_t zContained = false;
  Bool_t xyContained = false;
  Float_t r = TMath::Sqrt((x-xCenter)*(x-xCenter) + (y-yCenter)*(y-yCenter));
  if( z >= SuperModule1Beg && z <=SuperModule1End)
     zContained = true;
  if( r >= radialInner && r <= radialOuter)
     xyContained = true;

  Int_t retVal = 0;
  if(zContained && xyContained) retVal = 1;
  if(!zContained) retVal = -1;
  if(!xyContained) retVal -= 2;
  return retVal;
}


int NueConvention::IsInsideFarFiducial_Nue_Extended(float x, float y, float z)
{
  Float_t SuperModule1Beg =  0.35;
  Float_t SuperModule2Beg = 16.20;
  Float_t SuperModule1End = 14.57;
  Float_t SuperModule2End = 29.62;

  Float_t radialInner = 0.40;
  Float_t radialOuter = 3.87;
  Bool_t zContained = false;
  Bool_t xyContained = false;
  Float_t r = TMath::Sqrt(x*x + y*y);

  if( (z >= SuperModule1Beg && z <=SuperModule1End) ||
      (z >= SuperModule2Beg && z <=SuperModule2End) )
     zContained = true;

  if( r >= radialInner && r <= radialOuter)
     xyContained = true;

  Int_t retVal = 0;
  if(zContained && xyContained) retVal = 1;
  if(!zContained) retVal = -1;
  if(!xyContained) retVal -= 2;

  return retVal;  //  1 contained, -1 out of bounds z
                  //  -2 oob xy, -3 oob both
}

int NueConvention::IsInsideNearFiducial_Nue_Standard(float x, float y, float z, bool /*isMC*/)
{
  Float_t SuperModule1Beg = 1.01080;  //Data and MC values (according to DataUtil/infid.h on 10/02/07
  Float_t SuperModule1End = 4.99059;

  Float_t radialInner = 0;
  Float_t radialOuter = 0.8;
  Float_t xCenter = 1.4885;
  Float_t yCenter = 0.1397;

  Bool_t zContained = false;
  Bool_t xyContained = false;

  Float_t r = TMath::Sqrt((x-xCenter)*(x-xCenter) + (y-yCenter)*(y-yCenter));

  if( z >= SuperModule1Beg && z <=SuperModule1End)
     zContained = true;
  if( r >= radialInner && r <= radialOuter)
     xyContained = true;
                                                                                                                                     
  Int_t retVal = 0;
  if(zContained && xyContained) retVal = 1;
  if(!zContained) retVal = -1;
  if(!xyContained) retVal -= 2;

  return retVal;
}

int NueConvention::IsInsideFarFiducial_Nue_Standard(float x, float y, float z, bool isMC){

  Float_t SuperModule1Beg =  0.49080;   // These are data values
  Float_t SuperModule2Beg = 16.27110;
  Float_t SuperModule1End = 14.29300;
  Float_t SuperModule2End = 27.98270;

  if(isMC){
    SuperModule1Beg =  0.47692;   // These are mc values
    SuperModule2Beg = 16.26470;
    SuperModule1End = 14.27860;
    SuperModule2End = 27.97240;
  }
                                                                              
  Float_t radialInner = 0.50;
  Float_t radialOuter = TMath::Sqrt(14.0);
  Bool_t zContained = false;
  Bool_t xyContained = false;

  Float_t r = TMath::Sqrt(x*x + y*y);
                                                                                
  if( (z >= SuperModule1Beg && z <=SuperModule1End) ||
      (z >= SuperModule2Beg && z <=SuperModule2End) )
     zContained = true;
                                                                                
  if( r >= radialInner && r <= radialOuter)
     xyContained = true;
                                                                                
  Int_t retVal = 0;
  if(zContained && xyContained) retVal = 1;
  if(!zContained) retVal = -1;
  if(!xyContained) retVal -= 2;
                                                                                
  return retVal;  //  1 contained, -1 out of bounds z
                  //  -2 oob xy, -3 oob both
}

int NueConvention::IsInsideNearFiducial_MRE_Standard(float x, float y, float z,
						     bool /*isMC*/)
{
  Float_t SuperModule1Beg = 0.5;  //Data and MC values 
                                  //  (according to DataUtil/infid.h on 10/02/07
  Float_t SuperModule1End = 5.5;
                                                                                
  Float_t radialInner = 0;
  Float_t radialOuter = 1.2;
  Float_t xCenter = 1.4885;
  Float_t yCenter = 0.1397;
                                                                                
  Bool_t zContained = false;
  Bool_t xyContained = false;
                                                                                
  Float_t r = TMath::Sqrt((x-xCenter)*(x-xCenter) + (y-yCenter)*(y-yCenter));
                                                                                
  if( z >= SuperModule1Beg && z <=SuperModule1End)
     zContained = true;
  if( r >= radialInner && r <= radialOuter)
     xyContained = true;
                                                                                
                                                                                
  Int_t retVal = 0;
  if(zContained && xyContained) retVal = 1;
  if(!zContained) retVal = -1;
  if(!xyContained) retVal -= 2;
                                                                                
  return retVal;
}

// ********************************************************************
// Specialty Functions
// ********************************************************************

Int_t NueConvention::InPartialRegion(UShort_t plane, UShort_t strip){
  // figure out if this (plane,strip) corresponds to something in
  // the partially instrumented region
  //
  // this region is defined as:
  // v planes: (strip<=4 || strip>=67)
  // partial u: (strip==0 || strip=63)
  // full u: (strip<=26 || strip>=88)
  //
  // if so, return 1
  // if not, return -1
  // if error, return 0


  // make a lookup ptype to hold the type of each plane
  // 1 = v partial   2 = u partial
  // 3 = v full   4 = u full
  // 0 = uninstrumented
  static bool first=true;
  static UShort_t ptype[282];
  if(first){
    ptype[0]=0;
    for(int i=1; i<=281; i++){
      if(i%2==0) ptype[i]=1; // a v plane
      else ptype[i]=2; // a u plane
      if((i-1)%5 == 0) ptype[i]+=2; // fully instrumented
      else if(i>120) ptype[i]=0; // not instrumented
    }
    first=false;
  }
  if(plane>281){
    //    std::cerr<<"InPartialRegion passed plane = "<<plane<<std::endl;
    return 0;
  }
  UShort_t pt = ptype[plane];

  Int_t result;
  switch(pt){
  case 1:
  case 3:
    if(strip<=4 || strip>=67) result=1;
    else result = -1;
    break;
  case 2:
    if(strip==0 || strip == 63) result=1;
    else result = -1;
    break;
  case 4:
    if(strip<=26 || strip>=88) result=1;
    else result = -1;
    break;
  case 0:
  default:
    result=0;
    break;
  }
  return result;

}



// ***************************************************************
// NueEnergyCorrection Code
// - used to modify energy calibration
// ***************************************************************
void NueConvention::NueEnergyCorrection(NueRecord* nr, bool isLinearityFixMC)
{
//scale all GeV values using MEU variables and constants contained in this method

//if these are RHC files, should be using that correction
if(nr->anainfo.isRHC)
{
  RHCNueEnergyCorrection(nr);
  return;
}

int release = nr->GetHeader().GetRelease();
int detector = nr->GetHeader().GetVldContext().GetDetector();
bool isMC = (nr->GetHeader().GetVldContext().GetSimFlag()==SimFlag::kMC);

try
{
if(nr->srevent.phMeu>-9999)nr->srevent.phNueGeV=NueConvention::NueEnergyCorrection(nr->srevent.phMeu,release,isMC,detector, isLinearityFixMC);
if(nr->srshower.phMeu>-9999)nr->srshower.phNueGeV=NueConvention::NueEnergyCorrection(nr->srshower.phMeu,release,isMC,detector, isLinearityFixMC);
if(nr->srtrack.phMeu>-9999)nr->srtrack.phNueGeV=NueConvention::NueEnergyCorrection(nr->srtrack.phMeu,release,isMC,detector, isLinearityFixMC);

}
catch (int e)
{
  MSG("NueConvention",Msg::kError)<<"NueEnergyCorrection - attempted to correct energy for unknown specification "<<release<<endl;
}

}


float NueConvention::NueEnergyCorrection(float meu, int type, bool isMC, int detector=2, bool isLinearityFixMC = false)
{

 //
 //Cedar Carrot
 //
 if(ReleaseType::IsCedar(type)&&ReleaseType::IsCarrot(type))//cedar carrot only
 {
  float offset =   -1.52;
  float slope  =   25.06;
  return (meu-offset)/slope;
 }

 //
 //Cedar Daikon/Data
 //
 else if(ReleaseType::IsCedar(type)&&(ReleaseType::IsDaikon(type)||ReleaseType::IsData(type))&&(ReleaseType::GetRecoSubVersion(type)==0 || ReleaseType::GetRecoSubVersion(type)==1 ) ) //cedar daikon or data
 {
  float offset =   -1.515;
  float slope  =   24.86;
  return (meu-offset)/slope;
 }
   
 //
 //Cedar_phy(_bhcurv) Daikon/Data
 //
 else if(ReleaseType::IsCedar(type)&&(ReleaseType::IsDaikon(type)||ReleaseType::IsData(type))&&(ReleaseType::GetRecoSubVersion(type)==2 || ReleaseType::GetRecoSubVersion(type)==3)) //cedar_phy or cedar_phy_bhcurv, daikon or data
 {
  if(isMC)
  {
   if(isLinearityFixMC)
   {
    //The below numbers are for quasielastic nue events from MC with the FIX for linearity bug
    if(detector==2) //far
    {
     float offset =   -1.17129;
     float slope  =   24.2273;
     return (meu-offset)/slope;
    }
    else if(detector==1) //near
    {
     float offset =   -0.820177;
     float slope  =   23.7963;
     return (meu-offset)/slope;
    }
   }
   else
   {
    //The below numbers are for quasielastic nue events from MC with the linearity bug
    //Hence the number are only correct for quasielastic nue events from MC with the linearity bug
    //The NC, numu CC, and pid selected nue events will be effected by the linearity bug differently
    //Hence the below numbers are just plain wrong for everything but QE nue events
    if(detector==2) //far
    {
     //Numbers from Greg based on PETrimmerTest_V5 files with linearity bug (August 26, 2008: docDB 5001-v1)
     float offset =   -1.29686;
     float slope  =   24.3332;
     return (meu-offset)/slope;
    }
    else if(detector==1) //near
    {
     //New Numbers from Greg based on PETrimmerTest_V5 files with linearity bug  (August 26, 2008: docDB 5001-v1)
     float offset =   -2.8176;
     float slope  =   25.7362;
     return (meu-offset)/slope;
    }
   }
  }//if it is MC do the above
  else
  {
   //The below numbers are for quasielastic nue events from MC with the linearity bug FIX
   //To the best of our knowledge these are the correct energy scales for data
   //There are also addition fudge factors to compensate for imperfect GeVPerMip tunning and incorrect DB MEU numbers
   if(detector==2) //far
   {
    //Numbers from Greg based on PETrimmerTest_V5 files with linearity bug (August 26, 2008: docDB 5001-v1)
    float offset =   -1.17129*(0.9978/0.9988);
    float slope  =   24.2273*(0.9978/0.9988);
    return (meu-offset)/slope;
   }
   else if(detector==1) //near
   {
    //New Numbers from Greg based on PETrimmerTest_V5 files with linearity bug  (August 26, 2008: docDB 5001-v1)
    float offset =   -1.17129*(0.9978/1.0039);
    float slope  =   24.2273*(0.9978/1.0039);
    return (meu-offset)/slope;
   }
  }//else is data
 }//if cedar_phy or cedar_phy_bhcurv, daikon or data
  
 //
 //Dogwood Daikon/Data
 //
 else if(ReleaseType::IsDogwood(type)&& (ReleaseType::IsDaikon(type)||ReleaseType::IsData(type)))
 { 
  //The below numbers are from Xiaobo using Dogwood1 Daikon04 MC
  //These are preliminary numbers.  We will need to run over the full MC samples
  //and will need to take into account Data/MC fudge factors determined with 
  //stopping muon study with dogwood1
  if(detector==2) //far
  { 
   float offset_mip =   -2.483;
   float slope_MipPerGev  =  23.56;
   return (meu-offset_mip)/slope_MipPerGev;
  }
  else if(detector==1) //near
  {
   float offset_mip =   -1.071;
   float slope_MipPerGev  =  22.87;
   return (meu-offset_mip)/slope_MipPerGev;
  }
 }

 //
 //Elm Daikon/Data calling MINOSPLUS CORRECTOR
 //
 else if(ReleaseType::IsElm(type)&& (ReleaseType::IsDaikon(type)||ReleaseType::IsData(type)))
 {
  return MINOSPLUSNueEnergyCorrection(meu, type, isMC, detector, isLinearityFixMC);
 }//CALL TO MINOSPLUS correction if Elm.
 //This function's dependence on release type seems really shortsighted.

 //msg asked for ecal change but can't determine which to use

 throw -1;
 return -1;
}
//
//


//-------------------------------------------------------------------------------------------------------


//RHC RHC RHC RHC RHC RHC RHC RHC RHC
//RHC RHC RHC RHC RHC RHC RHC RHC HRC
//NueEnergyCorrection for RHC nuebar analysis
void NueConvention::RHCNueEnergyCorrection(NueRecord* nr, bool isLinearityFixMC)
{
//scale all GeV values using MEU variables and constants contained in this method


int release = nr->GetHeader().GetRelease();
int detector = nr->GetHeader().GetVldContext().GetDetector();
bool isMC = (nr->GetHeader().GetVldContext().GetSimFlag()==SimFlag::kMC);

 try
  {
   if(nr->srevent.phMeu>-9999)nr->srevent.phNueGeV=NueConvention::RHCNueEnergyCorrection(nr->srevent.phMeu,release,isMC,detector, isLinearityFixMC);
   if(nr->srshower.phMeu>-9999)nr->srshower.phNueGeV=NueConvention::RHCNueEnergyCorrection(nr->srshower.phMeu,release,isMC,detector, isLinearityFixMC);
   if(nr->srtrack.phMeu>-9999)nr->srtrack.phNueGeV=NueConvention::RHCNueEnergyCorrection(nr->srtrack.phMeu,release,isMC,detector, isLinearityFixMC);

  }
 
  catch (int e)
  {
  MSG("NueConvention",Msg::kError)<<"NueEnergyCorrection - attempted to correct energy for unknown specification "<<release<<endl;
  }

}

float NueConvention::RHCNueEnergyCorrection(float meu, int type, bool /*isMC*/, int detector=2, bool /*isLinearityFixMC*/)
{
 if(ReleaseType::IsDogwood(type)&& (ReleaseType::IsDaikon(type)||ReleaseType::IsData(type)))
 { 
  //The below numbers are from Adam using Dogwood3 Daikon07 MC
  //We consider these numbers valid for Dogwood5 Data, as the two reconstructions are similar
  if(detector==2) //far
  { 
   float offset_mip =   -3.493;
   float slope_MipPerGev  =  23.46;
   return (meu-offset_mip)/slope_MipPerGev;
  }
  else if(detector==1) //near
  {
   float offset_mip =   -3.933;
   float slope_MipPerGev  =  23.28;
   return (meu-offset_mip)/slope_MipPerGev;
  }
 }

 throw -1;
 return -1;
}
//
//

// NOvA ERA ENERGY CORRECTION BEGINS HERE -----------------------------------------------------------

void NueConvention::MINOSPLUSNueEnergyCorrection(NueRecord* nr, bool isLinearityFixMC)
{
//scale all GeV values using MEU variables and constants contained in this method


int release = nr->GetHeader().GetRelease();
int detector = nr->GetHeader().GetVldContext().GetDetector();
bool isMC = (nr->GetHeader().GetVldContext().GetSimFlag()==SimFlag::kMC);

 try
  {
   if(nr->srevent.phMeu>-9999)nr->srevent.phNueGeV=NueConvention::MINOSPLUSNueEnergyCorrection(nr->srevent.phMeu,release,isMC,detector, isLinearityFixMC);
   if(nr->srshower.phMeu>-9999)nr->srshower.phNueGeV=NueConvention::MINOSPLUSNueEnergyCorrection(nr->srshower.phMeu,release,isMC,detector, isLinearityFixMC);
   if(nr->srtrack.phMeu>-9999)nr->srtrack.phNueGeV=NueConvention::MINOSPLUSNueEnergyCorrection(nr->srtrack.phMeu,release,isMC,detector, isLinearityFixMC);

  }
 
  catch (int e)
  {
  MSG("NueConvention",Msg::kError)<<"MINOSPLUSNueEnergyCorrection - attempted to correct energy for unknown specification "<<release<<endl;
  }

}

float NueConvention::MINOSPLUSNueEnergyCorrection(float meu, int type, bool /*isMC*/, int detector=2, bool /*isLinearityFixMC*/)
{
 if(ReleaseType::IsElm(type)&& (ReleaseType::IsDaikon(type)||ReleaseType::IsData(type)))
 { 
  //The below numbers are from Adam using Elm5 Daikon10 MC
  //Note, at the behest of Andy Blake, we have inverted to flipping the axes
  //This is a 2-12 GeV fit
  if(detector==2) //far
  { 
    float offset_gev =  0.539882;
    float slope_GevPerMip  = 0.0391404;
    return (meu*slope_GevPerMip)+offset_gev;
  }
  else if(detector==1) //near
  {
    float offset_gev = 0.734848;
    float slope_GevPerMip  = 0.0383259;
    return (meu*slope_GevPerMip)+offset_gev;
  }
 }

 throw -1;
 return -1;
}


//---------------------------------------------------------------------------------------------------


//
//THINGS WE SHOULD APPARENTLY NEVER USE
void NueConvention::NueEnergyCorrectionNeverUseThisFunction(NueRecord* nr)
{
//scale all GeV values using MEU variables and constants contained in this method

int release = nr->GetHeader().GetRelease();
int detector = nr->GetHeader().GetVldContext().GetDetector();
bool isMC = (nr->GetHeader().GetVldContext().GetSimFlag()==SimFlag::kMC);

try
{
if(nr->srevent.phMeu>-9999)nr->srevent.phNueGeV=NueConvention::NueEnergyCorrectionNeverUseThisFunction(nr->srevent.phMeu,release,isMC,detector);
if(nr->srshower.phMeu>-9999)nr->srshower.phNueGeV=NueConvention::NueEnergyCorrectionNeverUseThisFunction(nr->srshower.phMeu,release,isMC,detector);
if(nr->srtrack.phMeu>-9999)nr->srtrack.phNueGeV=NueConvention::NueEnergyCorrectionNeverUseThisFunction(nr->srtrack.phMeu,release,isMC,detector);

}
catch (int e)
{
	MSG("NueConvention",Msg::kError)<<"NueEnergyCorrectionNeverUseThisFunction - attempted to correct energy for unknown specification\n";
}

}



float NueConvention::NueEnergyCorrectionNeverUseThisFunction(float meu, int type, bool isMC, int detector=2)
{
	if(ReleaseType::IsCedar(type)&&(ReleaseType::IsDaikon(type)||ReleaseType::IsData(type))&&(ReleaseType::GetRecoSubVersion(type)==2 || ReleaseType::GetRecoSubVersion(type)==3)) //cedar_phy or cedar_phy_bhcurv, daikon or data
	{
	  if(isMC)
	  {
	    //The below numbers are the average E-scale from MC with and without the linearity bug
	    //This is just plain wrong for everything
		  if(detector==2) //far
		  {
       	float offset =   ( (-1.29686)+(-1.17129) )/2.0;
       	float slope  =   (  (24.3332)+(24.2273)  )/2.0;
       	return (meu-offset)/slope;
		  }
		  else if(detector==1) //near
		  {
        float offset =   ( (-2.8176)+(-0.820177) )/2.0;
        float slope  =   ( (25.7362)+(23.7963)   )/2.0;
        return (meu-offset)/slope;
		  }
	  }//if it is MC do the above
	  else
	  {
	    //The below numbers are for quasielastic nue events from MC with the linearity bug FIX
	    //To the best of our knowledge these are the correct energy scales for data
	    //There are also addition fudge factors to compensate for imperfect GeVPerMip tunning and incorrect DB MEU numbers
		  if(detector==2) //far
		  {
        //Numbers from Greg based on PETrimmerTest_V5 files with linearity bug (August 26, 2008: docDB 5001-v1)
       	float offset =   -1.17129*(0.9978/0.9988);
       	float slope  =   24.2273*(0.9978/0.9988);
       	return (meu-offset)/slope;
		  }
		  else if(detector==1) //near
		  {
     	  //New Numbers from Greg based on PETrimmerTest_V5 files with linearity bug  (August 26, 2008: docDB 5001-v1)
       	float offset =   -1.17129*(0.9978/1.0039);
       	float slope  =   24.2273*(0.9978/1.0039);
        return (meu-offset)/slope;
		  }
	  }//else is data
	}//if cedar_phy or cedar_phy_bhcurv, daikon or data

	//msg asked for ecal change but can't determine which to use

	throw -1;
	return -1;
}

//***************************************************************
// File Position Code
// - used for scrolling through the Files
//***************************************************************

bool operator>(FilePosition one, FilePosition two)
{
   bool result = false;
   if(one.Run > two.Run)  result = true;
   if(one.Run == two.Run && one.SubRun > two.SubRun) result = true;
   if(one.Run == two.Run && one.SubRun == two.SubRun && one.Snarl > two.Snarl) result = true;
   if(one.Run == two.Run && one.SubRun == two.SubRun && one.Snarl == two.Snarl &&
        one.Event > two.Event) result = true;


   return result;
}


bool operator<(FilePosition one, FilePosition two)
{
   bool result = false;
   if(one.Run < two.Run)  result = true;
   if(one.Run == two.Run && one.SubRun < two.SubRun) result = true;
   if(one.Run == two.Run && one.SubRun == two.SubRun && one.Snarl < two.Snarl) result = true;
   if(one.Run == two.Run && one.SubRun == two.SubRun && one.Snarl == two.Snarl &&
        one.Event < two.Event) result = true;


   return result;
}

bool operator==(FilePosition one, FilePosition two)
{
   bool result = false;
   if(one.Run == two.Run && one.SubRun == two.SubRun && one.Snarl == two.Snarl &&
              one.Event == two.Event) result = true;

   return result;
}

//***************************************************************
// End of File Position Code
//***************************************************************
