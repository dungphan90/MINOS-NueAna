////////////////////////////////////////////////////////////////////////
// $Id: NueStandard.cxx,v 1.91 2017/02/27 18:13:05 wingmc Exp $
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
#include <vector>
#include "NueAna/NueAnaTools/NueConvention.h"
#include "NueStandard.h"
#include "DcsUser/CoilTools.h"
#include "MCReweight/SKZPWeightCalculator.h"
#include "TMultiLayerPerceptron.h"
#include "MCNNAnalysis/LEMNNbarpid.h"
#include "MCNNAnalysis/LEMbarpid.h"
#include "MCNNAnalysis/LEM4pid.h"
#include "MCNNAnalysis/LEMAmby.h"
#include "MCNNAnalysis/LEMLSND.h"
#include "MCNNAnalysis/LEMAmbyE50N491.h"
#include "MCNNAnalysis/LEMAmbyE50N591.h"
#include "MCNNAnalysis/LEMAmbyE50N691.h"
#include "MCNNAnalysis/LEMAmbyE50N5111.h"
#include "MCNNAnalysis/LEMAmbyE50N6111.h"
#include "MCNNAnalysis/LEMAmbyE50S491.h"
#include "MCNNAnalysis/LEMAmbyE50S591.h"
#include "MCNNAnalysis/LEMAmbyE50S691.h"
#include "MCNNAnalysis/LEMAmbyE50S5111.h"
#include "MCNNAnalysis/LEMAmbyE50S6111.h"

#include "TFile.h"
#include "TMath.h"
#include <cmath>
#include <string>
#include <vector>
#include <iostream>

#include <fstream>
#include "MessageService/MsgService.h"

#include "NueAna/Extrapolation/Background.h"

CVSID("");

namespace NueStandard {
  static bool isNSI = false;
}

void NueStandard::SetNSI(bool nsiflag) { 
  isNSI = nsiflag; 
}

namespace NueStandard {
  static bool isLSND = false;
}

void NueStandard::SetLSND(bool lsndflag) { 
  isLSND = lsndflag; 
}

bool NueStandard::IsInFid(NueRecord *nr)
{    
  float x = nr->srevent.vtxX;
  float y = nr->srevent.vtxY;
  float z = nr->srevent.vtxZ;
  
  Detector::Detector_t fDet;
  fDet = nr->GetHeader().GetVldContext().GetDetector();
  if (fDet != Detector::kFar && fDet != Detector::kNear) return false;
  
  SimFlag::SimFlag_t sim = nr->GetHeader().GetVldContext().GetSimFlag();
  bool isMC = (sim == SimFlag::kMC);
  int retval = 0;
  
  if (fDet == Detector::kFar)
    retval = NueConvention::IsInsideFarFiducial_Nue_Standard(x,y,z, isMC);
    
  if (fDet == Detector::kNear)
     retval = NueConvention::IsInsideNearFiducial_Nue_Standard(x,y,z, isMC);
        
   return (retval == 1);

}

bool NueStandard::PassesPOTStandards(NueRecord *nr)
{
   SimFlag::SimFlag_t sim = nr->GetHeader().GetVldContext().GetSimFlag();
   bool isMC = (sim == SimFlag::kMC);
   VldContext vld = nr->GetHeader().GetVldContext();
   int sec = nr->GetHeader().GetVldContext().GetTimeStamp().GetSec(); 
   if (isMC) return true;
   
   //common to both
   if (nr->bmon.goodBeamMon != 1) return false;

   
  //FHC ND
  //if(ND && (Run1-3 || Run5-6 || Run8 || Run10+))
  if (nr->GetHeader().GetVldContext().GetDetector() == Detector::kNear && (sec < 1254248760 || (sec > 1269298800 && sec < 1278924900) || (sec > 1299011200 && sec < 1311993300) || (sec > 1317803440) ) )
   {
     bool goodCoil =  nr->eventq.coilQuality && (nr->eventq.coilDirection ==1);
                      //CoilTools::IsOK(vld) && !CoilTools::IsReverse(vld);
     int dpDQ = nr->eventq.passNearDetQuality;
     return goodCoil && (dpDQ == 1);
   }
   
  //RHC ND
  //if(ND && (Run4 || Run7 || Run9))
  if (nr->GetHeader().GetVldContext().GetDetector() == Detector::kNear && ( (sec > 1254248760 && sec < 1269298800)|| (sec > 1278924900 && sec < 1299011200) || (sec > 1311993300 && sec < 1317803440) ) )
   {
     bool goodCoil =  nr->eventq.coilQuality && (nr->eventq.coilDirection ==-1);
     int dpDQ = nr->eventq.passNearDetQuality;
     return goodCoil && (dpDQ == 1);
   }

   if (nr->GetHeader().GetVldContext().GetDetector() == Detector::kFar)
   {

    if (nr->eventq.spillType != 1) return false;
    int dpDQ = nr->eventq.passFarDetQuality;
    if (nr->GetHeader().GetEventNo() < 0) return (dpDQ == 1);

     return (dpDQ == 1) && (nr->eventq.passLI == 1);
   }

   return false;
}                                                                                                        
bool NueStandard::PassesDataQuality(NueRecord *nr)
{                                                 
   SimFlag::SimFlag_t sim = nr->GetHeader().GetVldContext().GetSimFlag();
   bool isMC = (sim == SimFlag::kMC);
   VldContext vld = nr->GetHeader().GetVldContext();

   if (isMC) return true;
   
   if(!PassesPOTStandards(nr)) return false;

   if(nr->GetHeader().GetVldContext().GetDetector() == Detector::kNear) return true;

   if (nr->GetHeader().GetVldContext().GetDetector() == Detector::kFar){

     int rcB = nr->eventq.rcBoundary;
     float tmin = nr->srevent.eventTimeMin;
     double spillT = nr->bmon.dt_stnd;   
     return PassesFarDataQuality(-1, rcB, 1, tmin, spillT);       
   }

   return false;
}

bool NueStandard::PassesNearDataQuality(int gbm, float cc, int st)
{
  //Assumes it is taking:
  // nr->bmon.goodBeamMon, nr->srevent.coilCurrent, nr->srevent.spillType
  if (gbm == 1 && cc < -1000.0 && st != 3)
     return true;
                                                                                
  return false;
}

bool NueStandard::PassesFarDataTiming(NueRecord *nr) 
{
 if(
       Detector::kFar == nr->GetHeader().GetVldContext().GetDetector()
    && SimFlag::kData == nr->GetHeader().GetVldContext().GetSimFlag()
  )
 {
  float tmin = nr->srevent.eventTimeMin;
  double spillT = nr->bmon.dt_stnd;
  return( (tmin - spillT) > -2e-6 && (tmin - spillT) < 12e-6 );
 }
  return true;
}

bool NueStandard::PassesFarDataQuality(float li, int rc, int dpfddq, 
                                          float tmin, double spillt)
{
  if (li == -1 && rc == 0 && dpfddq == 1
    && (tmin - spillt) > -20e-6 && (tmin - spillt) < 30e-6)
     return true;
                                                                                
  return false;
}

bool NueStandard::IsGoodRun(NueRecord *nr)
{
  if(nr->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC)
  {
    return true;
  }
  
  if(nr->GetHeader().GetVldContext().GetDetector() == Detector::kFar)
  {
    return NueStandard::IsGoodFarRun(nr);
  }
  
  return false;
}

bool NueStandard::IsGoodNearRun(NueRecord *nr)
{
  if(nr->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC)
  {
    return true;
  }
  
  if(nr->GetHeader().GetVldContext().GetDetector() == Detector::kFar)
  {
    return true;
  }
  
  //for near detector data
  if(nr->GetHeader().GetRun()==8165 || nr->GetHeader().GetRun()==11318)
  {
    return false;
  }
  
  if(nr->GetHeader().GetRun()==7942 && nr->GetHeader().GetSubRun()==0)
  {
    return false;
  }
  
  if(nr->GetHeader().GetRun()==7982 && nr->GetHeader().GetSubRun()==0)
  {
    return false;
  }
  
  if(nr->GetHeader().GetRun()==8214 && nr->GetHeader().GetSubRun()==0)
  {
    return false;
  }
  
  if(nr->GetHeader().GetRun()==9809 && nr->GetHeader().GetSubRun()==0)
  {
    return false;
  }
  
  return true;
}

bool NueStandard::IsGoodFarRun(NueRecord *nr)
{
  if(nr->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC)
  {
    return true;
  }
  
  if(nr->GetHeader().GetVldContext().GetDetector() == Detector::kNear)
  {
    return true;
  }
  
  //for far detector data
  static vector<FilePosition> badFDrunlist;
  static unsigned int listpos = 0;
  int Run, SubRun;
  
  // If this is the first call fill the file list
  if(badFDrunlist.size() == 0)
  {
    std::ifstream ins;
    char *srt_dir = getenv("SRT_PRIVATE_CONTEXT");
    ins.open(Form("%s/NueAna/badrunsFD_list.txt",srt_dir));
    if(!ins.is_open())
    {
      MSG("NueStandard", Msg::kError)<<"IsGoodFarRun(): Unable to open badrunsFD_list.txt"<<endl;
    }

    while(!ins.eof())
    {
      ins>>Run>>SubRun;
      FilePosition temp;
      if(!ins.eof()){
          temp.Run = Run;
          temp.SubRun = SubRun;
          temp.Snarl = -1;
          temp.Event = -1;
          badFDrunlist.push_back(temp);
      }
    }
    
    if(badFDrunlist.size() == 0)  return true;
   }
   
   FilePosition current;
   
   current.Snarl = -1;
   current.Event = -1;
   current.Run = nr->GetHeader().GetRun();
   current.SubRun = nr->GetHeader().GetSubRun();
   
   if(current < badFDrunlist[0]) return true;

   if(current > badFDrunlist[badFDrunlist.size() - 1]) return true;

   if(current > badFDrunlist[listpos] )
   {
      while(current > badFDrunlist[listpos] && listpos < badFDrunlist.size())
          listpos++;
   }
   else if(current < badFDrunlist[listpos] )
   {
      while(current < badFDrunlist[listpos] && listpos > 0)
          listpos--;
   }

   if(current == badFDrunlist[listpos])
   {
       return false;
   }

   
   return true;
}

//  Now to the PreSelection Cuts

bool NueStandard::PassesTrackPlaneCut(NueRecord *nr)
{
   int tp = TMath::Abs(nr->srtrack.endPlane - nr->srtrack.begPlane);

   return PassesTrackPlaneCut(tp);
}

bool NueStandard::PassesTrackPlaneCut(int trkplane)
{
   return (trkplane < 25);
}

bool NueStandard::PassesMinPlaneCut(NueRecord *nr)
{
   int planes = nr->shwfit.contPlaneCount050;
   return PassesMinPlaneCut(planes);
}

bool NueStandard::PassesMinPlaneCut(int planes)
{
   return (planes > 4);
}

bool NueStandard::PassesShowerCut(NueRecord *nr)
{
  int nshw = nr->srevent.showers;
  return PassesShowerCut(nshw);
}

bool NueStandard::PassesShowerCut(int nshw)
{
  return (nshw > 0);
}

bool NueStandard::PassesCosmicCutFunction(NueRecord *nr)
{
   bool result=true;	
   //cosmicy shower cut
   result = result && !(nr->shwfit.UVSlope > 10);

   //cosmicy track cut
   if (nr->srevent.tracks<1) return result;  //require a track to apply this cut

   float ex =  nr->srtrack.endX;
   float ey =  nr->srtrack.endY;
   float ez =  nr->srtrack.endZ;

   float vx =  nr->srtrack.vtxX;
   float vy =  nr->srtrack.vtxY;
   float vz =  nr->srtrack.vtxZ;

   float cosy = acos(fabs(vy-ey)/sqrt((vx-ex)*(vx-ex)+(vy-ey)*(vy-ey)+(vz-ez)*(vz-ez))); 
   float dy=fabs(vy-ey);

   result = result && !(dy>2 && cosy < 0.6);

   return result;
}

bool NueStandard::PassesCosmicCut(NueRecord *nr)
{
  if(nr->GetHeader().GetVldContext().GetDetector() == Detector::kNear)
       return true;

   bool ret=true;
   if(nr->eventq.passCosmicCut > -10)
     return (nr->eventq.passCosmicCut == 1);

   //or we are still in Ent
 
   if(nr->dtree.bt_var1==0)ret = false;
   else if(nr->dtree.bt_var1==1)ret = true;
   else{ NueStandard::FillCosmicCut(nr); ret = nr->dtree.bt_var1; }

   Detector::Detector_t fDet;
   fDet = nr->GetHeader().GetVldContext().GetDetector();
   if (fDet != Detector::kFar) ret = true;

   return ret;
}

void NueStandard::FillCosmicCut(NueRecord *nr)
{
   int bv = 0;
   if (NueStandard::PassesCosmicCutFunction(nr)) bv=1;
   nr->dtree.bt_var1=bv;
}

bool NueStandard::PassesTrackLikePlaneCut(NueRecord *nr)
{
   int tlp = nr->srtrack.trklikePlanes;
   return PassesTrackLikePlaneCut(tlp);
}

bool NueStandard::PassesTrackLikePlaneCut(int tlp)
{       
   return (tlp < 16);
}

bool NueStandard::PassesLowEnergyCut(NueRecord *nr)
{
   float energy = nr->srevent.phNueGeV;
   return PassesLowEnergyCut(energy);
}
 
bool NueStandard::PassesLowEnergyCut(float energy)
{
   return (energy > 1.0);
}

bool NueStandard::PassesHighEnergyCut(NueRecord *nr)
{
   float energy = nr->srevent.phNueGeV;
   return PassesHighEnergyCut(energy);
}
          
bool NueStandard::PassesHighEnergyCut(float energy)
{
   return (energy < 12.0);
}

// Now to combinations of Preselection Cuts                                                                                                          
bool NueStandard::PassesPreSelection(NueRecord *nr)
{
   bool temp = PassesNonHEPreSelection(nr);
   temp = temp && PassesHighEnergyCut(nr);
  
   return temp;
}

bool NueStandard::PassesPreSelection(int trkplane, int tlp, float energy)
{
   bool temp = PassesNonHEPreSelection(trkplane, tlp, energy);
   temp = temp && PassesHighEnergyCut(energy);
   return temp;
}

bool NueStandard::PassesNonHEPreSelection(NueRecord *nr)
{
   bool temp = PassesPreSelectionTrackCuts(nr);
   temp = temp && PassesLowEnergyCut(nr);
   temp = temp && PassesPreSelectionBasicCuts(nr);
   
   return temp;
}                                                                                                       
bool NueStandard::PassesNonHEPreSelection(int trkplane, int tlp, float energy)
{
   bool temp = PassesPreSelectionTrackCuts(trkplane, tlp);
   temp = temp && PassesLowEnergyCut(energy);
   return temp;
}

bool NueStandard::PassesPreSelectionTrackCuts(NueRecord *nr)
{
   bool temp = PassesTrackPlaneCut(nr);
   temp = temp && PassesTrackLikePlaneCut(nr);
   return temp;  
}

bool NueStandard::PassesPreSelectionTrackCuts(int trkplane, int tlp)
{
   bool temp = PassesTrackPlaneCut(trkplane);
   temp = temp && PassesTrackLikePlaneCut(tlp);
   return temp;
}

bool NueStandard::PassesPreSelectionBasicCuts(NueRecord *nr)
{
   bool one = PassesMinPlaneCut(nr);
   one = one && PassesCosmicCut(nr); 
   one = one && IsLargestEvent(nr);
   one = one && PassesShowerCut(nr);

   return one;
}

bool NueStandard::PassesMRCCFiducial(NueRecord *nr)
{
   float x = nr->mri.vtxx;
   float y = nr->mri.vtxy;
   float z = nr->mri.vtxz;
                                                                                
   Detector::Detector_t fDet;
   fDet = nr->GetHeader().GetVldContext().GetDetector();
   
   SimFlag::SimFlag_t sim = nr->GetHeader().GetVldContext().GetSimFlag();
   bool isMC = (sim == SimFlag::kMC);
   int retval = 0;
                                                                                
   if (fDet == Detector::kFar)
   {
     if(sqrt((x*x) + (y*y))>0.3 && sqrt((x*x) + (y*y))<sqrt(15.5) && ((z>0.3 && z<14.4)||(z>16.1 && z<28.2)))
     {
       retval=1;
     }
   }
       
   if (fDet == Detector::kNear)
     retval = NueConvention::IsInsideNearFiducial_MRE_Standard(x,y,z, isMC);

   return (retval == 1);
}

bool NueStandard::PassesMRCCPreSelection(NueRecord *nr)
{
    double pid = nr->mri.orig_roCCPID;
    int fitp = nr->mri.fitp;
    float bestC = nr->mri.best_complete;

    return NueStandard::PassesMRCCPreSelection(bestC, fitp, pid);
}

bool NueStandard::PassesMRCCPreSelection(float bestComp, int fitp, double pid)
{
   //currently assumes RO-PID
   bool good = true;
   good = good && (bestComp > -10);
   good = good && (fitp == 1);
   good = good && (pid > 0.3);
                                                                                
   return good;
}

bool NueStandard::PassesMREFiducial(NueRecord *nr)
{
    return NueStandard::PassesMRCCFiducial(nr);
}

bool NueStandard::PassesMREPreSelection(NueRecord *nr)
{
    return NueStandard::PassesMRCCPreSelection(nr);
}
                                                                                                              
bool NueStandard::PassesMREPreSelection(float bestC, int fitp, double pid)
{
    return NueStandard::PassesMRCCPreSelection(bestC, fitp, pid);
}
                                                                                             
bool NueStandard::PassesSysPreSelectionNoHE(NueRecord *nr)
{
   int tp = TMath::Abs(nr->srtrack.endPlane - nr->srtrack.begPlane);
   int tlp = nr->srtrack.trklikePlanes;
   float energy = nr->srevent.phNueGeV;

   return PassesSysPreSelectionNoHE(tp, tlp, energy);
}
                                                                                             
bool NueStandard::PassesSysPreSelectionNoHE(int trkplane, int tlp, float energy)
{
   bool good = true;
   good = good && (trkplane < 28);
   good = good && (tlp < 18);
   good = good && (energy > 0.5);
   return good;
}


bool NueStandard::PassesSysPreSelection(NueRecord *nr)
{
   int tp = TMath::Abs(nr->srtrack.endPlane - nr->srtrack.begPlane);
   int tlp = nr->srtrack.trklikePlanes;
   float energy = nr->srevent.phNueGeV;
                                                                                
   return PassesSysPreSelection(tp, tlp, energy);
}
                                                                                
bool NueStandard::PassesSysPreSelection(int trkplane, int tlp, float energy)
{
   bool good = true;
   good = good && (trkplane < 28);
   good = good && (tlp < 18);
   good = good && (energy > 0.5 && energy < 10);
   return good;
}

double NueStandard::GetPIDValue(NueRecord *nr, Selection::Selection_t sel)
{
  double retval = -9999;

  switch(sel) {
  case Selection::kCuts:  retval = nr->treepid.fCutPID;           break;
  case Selection::kANN6:  retval = nr->ann.pid_6inp;              break;
  case Selection::kANN30:  retval = nr->ann.pid_30inp;            break;
  case Selection::kANN2PE:  retval = nr->ann.pid_11inp;           break;
  case Selection::kANN2PE_DAIKON04:  retval = nr->ann.pid_11inp_daikon04;         break;
  case Selection::kANN14_DAIKON04:   retval = nr->ann.pid;        break;
//  case Selection::kANN14_DAIKON04:   retval = nr->ann.pid_14inp_daikon04;         break;
  case Selection::kANN4FHC: 
    retval = NueStandard::Calc4thAnaANN(nr, Selection::kANN4FHC); break;
  case Selection::kANN4RHC: 
    retval = NueStandard::Calc4thAnaANN(nr, Selection::kANN4RHC); break;
  case Selection::kSSPID: retval = nr->subshowervars.pid;  break;
  case Selection::kMCNN:
     if (nr->mcnnv.bestmatches < 50)  retval = -0.5;
     else        retval = nr->mcnnv.mcnn_var1;                         break;   
  
  case Selection::kLEMNNBAR:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
	//Merged Method
	retval = nr->mcnnv.mcnn_var2;
        
	//Stand Alone Selector Method
	//LEMNNbarpid lemnnbarpid;
	//retval = lemnnbarpid.GetLEMNNbarpid(nr);
           }                                                           break;

  case Selection::kLEMBAR:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
	LEMbarpid lembarpid;
        retval = lembarpid.GetLEMbarpid(nr);
      }                                                                break;

  case Selection::kLEMAmbyE50N491:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        retval = nr->mcnnv.mcnn_var3;       
        //Stand Alone Selector Method Obsolete at Trim
        //LEMAmbyE50N491 lemambye50n491;
	//retval = lemambye50n491.GetLEMAmbyE50N491(nr);
      }                                                                break;


  case Selection::kLEMAmbyE50S491:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        retval = nr->mcnnv.mcnn_var4;
        //Stand Alone Selector Method Obsolete at Trim
        //LEMAmbyE50S491 lemambye50s491;
        //retval = lemambye50s491.GetLEMAmbyE50S491(nr);
      }                                                                break;

  case Selection::kLEMAmbySAE50S491:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        LEMAmbyE50S491 lemambye50s491;
        retval = lemambye50s491.GetLEMAmbyE50S491(nr);
      }                                                                break;

      //ALL kLEMAmby/LSND's below this point are obsolete trials and should not be used in analysis

  case Selection::kLEM4:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
	//Stand Alone Selector Method
	LEM4pid lem4pid;
        retval = lem4pid.GetLEM4pid(nr);
      }                                                                break;

  case Selection::kLEMAmby:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
        LEMAmby lemamby;
        retval = lemamby.GetLEMAmby(nr);
      }                                                                break;

  case Selection::kLEMLSND:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
        LEMLSND lemlsnd;
        retval = lemlsnd.GetLEMLSND(nr);
      }                                                                break;

  case Selection::kLEMAmbyE50N591:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
        LEMAmbyE50N591 lemambye50n591;
        retval = lemambye50n591.GetLEMAmbyE50N591(nr);
      }                                                                break;

  case Selection::kLEMAmbyE50N691:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
        LEMAmbyE50N691 lemambye50n691;
        retval = lemambye50n691.GetLEMAmbyE50N691(nr);
      }                                                                break;

  case Selection::kLEMAmbyE50N5111:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
        LEMAmbyE50N5111 lemambye50n5111;
        retval = lemambye50n5111.GetLEMAmbyE50N5111(nr);
      }                                                                break;

   case Selection::kLEMAmbyE50N6111:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
        LEMAmbyE50N6111 lemambye50n6111;
        retval = lemambye50n6111.GetLEMAmbyE50N6111(nr);
      }                                                                break;


  case Selection::kLEMAmbyE50S591:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
        LEMAmbyE50S591 lemambye50s591;
        retval = lemambye50s591.GetLEMAmbyE50S591(nr);
      }                                                                break;

  case Selection::kLEMAmbyE50S691:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
        LEMAmbyE50S691 lemambye50s691;
        retval = lemambye50s691.GetLEMAmbyE50S691(nr);
      }                                                                break;

  case Selection::kLEMAmbyE50S5111:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
        LEMAmbyE50S5111 lemambye50s5111;
        retval = lemambye50s5111.GetLEMAmbyE50S5111(nr);
      }                                                                break;

   case Selection::kLEMAmbyE50S6111:
      if (nr->mcnnv.bestmatches < 50) retval = -0.5;
      else{
        //Stand Alone Selector Method
        LEMAmbyE50S6111 lemambye50s6111;
        retval = lemambye50s6111.GetLEMAmbyE50S6111(nr);
      }                                                                break;


  case Selection::kParticlePID:
  	if(PassesParticlePIDPreselectionCut(nr))retval = nr->precord.event.pidF;
 	break;	
  
  default:                                                        break;
  }

  return retval;
}

bool NueStandard::PassesPIDSelection(NueRecord *nr, Selection::Selection_t sel)
{
  bool pass = false;
  double val = GetPIDValue(nr, sel);
                                                                                                  
  switch(sel) {
  case Selection::kCuts:  pass = (int(val) == 1);       break;
  case Selection::kANN6:   pass = val > 0.5;            break;
  case Selection::kANN30:  pass = val > 0.5;            break;
  case Selection::kANN2PE:  pass = val > 0.7;           break;
  case Selection::kANN2PE_DAIKON04:  pass = val > 0.7;  break;
  case Selection::kANN14_DAIKON04:   pass = val > 0.75; break;
  case Selection::kANN4FHC: pass = val > 0.7;           break;
  case Selection::kANN4RHC: pass = val > 0.7;           break;
  case Selection::kSSPID: pass = val > 0.67;            break;
  case Selection::kMCNN:  pass = val > 0.80;            break;
  case Selection::kParticlePID:  pass = PassesParticlePIDCut(nr); break;
  case Selection::kLEMAmbyE50S491:  pass = val > 0.70;  break;
  default:                                              break;
  }
                                                                                                  
  return pass;
}

bool NueStandard::PassesParticlePIDCut(NueRecord *nr)
{
	bool passPID=true;
	
	passPID = passPID && nr->precord.event.pidF>0.7;

	if(!passPID)return false;


	return passPID && PassesParticlePIDPreselectionCut(nr);
	

}


bool NueStandard::PassesParticlePIDPreselectionCut(NueRecord *nr)
{
		
	//fiducial cut
	bool passFid=true;
	passFid = passFid && (nr->precord.event.inFiducial==1); 

	if(!passFid)return false;


	double visenergy = nr->precord.event.visenergy;
		
	//preselection cut
	int det = nr->GetHeader().GetVldContext().GetDetector();
	if(det == Detector::kFar)
	{
		float offset =   0.489987;
   	 	float slope  =  0.0387301;
        visenergy = visenergy*slope+offset;
	}else if(det == Detector::kNear)
	{
		float offset = 0.4803261;
		float slope = 0.03799819;
		visenergy = visenergy*slope+offset;
	}else{
		printf("please don't run in a detector that is not near or far (NueStandard::PassesParticlePIDCut)\n");
		exit(1);
	}
		
	bool passPre=true;
		
	double event_length=nr->precord.event.max_z-nr->precord.event.min_z;	
	passPre = passPre && event_length>0.1 && event_length<1.2;
	passPre = passPre && nr->precord.particles.longest_s_particle_s>0.1 && 
		nr->precord.particles.longest_s_particle_s<1.2;
	passPre = passPre && nr->precord.particles.ntot>0;
	passPre = passPre && visenergy>0.5 && visenergy<8;	

	return passPre && passFid;

}


bool NueStandard::IsLargestEvent(NueRecord *nr)
{
 if(nr->GetHeader().GetVldContext().GetDetector() == Detector::kFar){
    if(nr->srevent.largestEvent == 1) return true;
    return false;
 }
 return true;
}

bool NueStandard::PassesSelection(NueRecord *nr, Selection::Selection_t sel)
{
  bool pass = false;

  bool dq = NueStandard::PassesDataQuality(nr);
  bool infid = NueStandard::IsInFid(nr) && dq;
  
  bool isMR = nr->GetHeader().FoundMR();

  if (isMR){
    bool mrfid = PassesMRCCFiducial(nr);
    bool mrps = PassesMRCCPreSelection(nr);

    infid = infid && mrfid && mrps;
  } 
  
  bool basic = NueStandard::PassesPreSelectionBasicCuts(nr) && infid;
  bool presel = infid && NueStandard::PassesPreSelection(nr);
  bool pid = NueStandard::PassesPIDSelection(nr,sel); 

  switch(sel) {
  case Selection::kNone:     pass = true;                 break;
  case Selection::kDataQual: pass = dq;                   break;
  case Selection::kFid:      pass = infid;                break;
  case Selection::kBasic:    pass = basic;                break;
  case Selection::kPre:   pass = presel; 	          break;
  case Selection::kCuts:  pass = presel && pid;           break;
  case Selection::kANN6:  pass = presel && pid;           break;
  case Selection::kANN30: pass = presel && pid;           break;
  case Selection::kANN2PE: pass = presel && pid;          break;
  case Selection::kANN2PE_DAIKON04: pass = presel && pid; break;
  case Selection::kANN14_DAIKON04: pass = presel && pid;  break;
  case Selection::kANN4FHC: pass = presel && pid;         break;
  case Selection::kANN4RHC: pass = presel && pid;         break;
  case Selection::kSSPID: pass = presel && pid;           break;
  case Selection::kMCNN:  pass = presel && pid;           break;
  case Selection::kMDA:   pass = presel && false;         break;
  case Selection::kBDT:   pass = presel && false;         break;
  case Selection::kKNue:  pass = presel && false;         break;
  case Selection::kCC:    pass = infid && 
                     NueStandard::PassesCCSelection(nr);  break;
  default:                                                break;
  }

  return pass;
}

bool NueStandard::PassesSelection(NueRecord *nr, Selection::Selection_t sel, Selection::Selection_t sel2)
{
  bool pass = false;

  bool dq = NueStandard::PassesDataQuality(nr);
  bool infid = NueStandard::IsInFid(nr) && dq;
  bool basic = NueStandard::PassesPreSelectionBasicCuts(nr);
  bool presel = NueStandard::PassesPreSelection(nr);
  bool pid = NueStandard::PassesPIDSelection(nr,sel);

  switch(sel2) {
  case Selection::kDataQual:  dq = true;                  break;
  case Selection::kFid:       infid = true;               break;
  case Selection::kBasic:     basic = true;               break;
  case Selection::kPre:       presel = true;              break;
  default:                                                break;
  }
  
  bool isMR = nr->GetHeader().FoundMR();

  if(isMR){
    bool mrfid = true; //PassesMRCCFiducial(nr);
    bool mrps = PassesMRCCPreSelection(nr);

    infid = infid && mrfid && mrps;
  }
  
  basic = basic && infid;
  presel = dq && infid && presel;

  switch(sel) {
  case Selection::kNone:  pass = true;                    break;
  case Selection::kDataQual:  pass = dq;                  break;
  case Selection::kFid:   pass = infid;                   break;
  case Selection::kBasic: pass = basic;                   break;
  case Selection::kPre:   pass = presel;                  break;
  case Selection::kCuts:  pass = presel && pid;           break;
  case Selection::kANN6:   pass = presel && pid;          break;
  case Selection::kANN30:  pass = presel && pid;          break;
  case Selection::kANN2PE: pass = presel && pid;          break;
  case Selection::kANN2PE_DAIKON04: pass = presel && pid; break;
  case Selection::kANN14_DAIKON04: pass = presel && pid;  break;
  case Selection::kANN4FHC: pass = presel && pid;         break;
  case Selection::kANN4RHC: pass = presel && pid;         break;
  case Selection::kSSPID: pass = presel && pid;           break;
  case Selection::kMCNN:  pass = presel && pid;           break;
  case Selection::kMDA:   pass = presel && false;         break;
  case Selection::kBDT:   pass = presel && false;         break;
  case Selection::kKNue:  pass = presel && false;         break;
  case Selection::kCC:    pass = infid &&
                     NueStandard::PassesCCSelection(nr);  break;
  default:                                                break;
  }

  return pass;
}


bool NueStandard::PassesCCSelection(NueRecord *nr)
{
   int ntrack  = nr->srevent.tracks;
   int pass = nr->srtrack.passedFit;
   int endPlaneU = nr->srtrack.endPlaneU;
   int endPlaneV = nr->srtrack.endPlaneV;
   int deltaUVVtx = nr->srtrack.deltaUVVtx;   

   if (ntrack < 1){ pass = 0; endPlaneU = endPlaneV = deltaUVVtx = 300; }

   bool trackPass =  (pass == 1) || 
                (deltaUVVtx <= 5 && TMath::Abs(endPlaneU - endPlaneV) <= 40
                   && endPlaneU < 270 && endPlaneV < 270);


   return  (ntrack > 0) && trackPass && nr->anainfo.roCCPID > 0.3;
}


double NueStandard::GetRPWBeamWeight(NueRecord *nr, bool ismrcc)
{


  std::vector<double> pots;
  
  Detector::Detector_t det = nr->GetHeader().GetVldContext().GetDetector();
  BeamType::BeamType_t beam = nr->GetHeader().GetBeamType();

  if(beam != BeamType::kL010z185i
     &&beam != BeamType::kL010z185i_i124 
     &&beam != BeamType::kL010z185i_i191   
     &&beam != BeamType::kL010z185i_i213
     &&beam != BeamType::kL010z185i_i224
     &&beam != BeamType::kL010z185i_i232
     &&beam != BeamType::kL010z185i_i243
     &&beam != BeamType::kL010z185i_i257
     &&beam != BeamType::kL010z185i_i282
     &&beam != BeamType::kL010z185i_i303
     &&beam != BeamType::kL010z185i_i324
   
     &&beam != BeamType::kL010z000i
     &&beam != BeamType::kL010z000i_i209
     &&beam != BeamType::kL010z000i_i225
     &&beam != BeamType::kL010z000i_i232
     &&beam != BeamType::kL010z000i_i259
     &&beam != BeamType::kL010z000i_i300
     &&beam != BeamType::kL010z000i_i317
     &&beam != BeamType::kL010z000i_i326
     &&beam != BeamType::kL010z000i_i380){

    MAXMSG("NueStandard",Msg::kWarning,5)<<"Can't use NueStandard::GetRPWBeamWeight for non L010185 or L010000 beams, stop it!"<<endl;
    return 0;
  }

  if (det==Detector::kNear){
    if (ismrcc){
      for(int i=0;i<NRUNPERIODS;i++){
	pots.push_back(pot_ndmrcc[i]);
      }
    }
    else if (beam == BeamType::kL010z185i
     ||beam == BeamType::kL010z185i_i124
     ||beam == BeamType::kL010z185i_i191
     ||beam == BeamType::kL010z185i_i213
     ||beam == BeamType::kL010z185i_i224
     ||beam == BeamType::kL010z185i_i232
     ||beam == BeamType::kL010z185i_i243
     ||beam == BeamType::kL010z185i_i257
     ||beam == BeamType::kL010z185i_i282
     ||beam == BeamType::kL010z185i_i303
     ||beam == BeamType::kL010z185i_i324){
      for(int i=0;i<NRUNPERIODS;i++){
	pots.push_back(pot_nd[i]);
      }
    }
    else if (beam==BeamType::kL010z000i
     ||beam == BeamType::kL010z000i_i209
     ||beam == BeamType::kL010z000i_i225
     ||beam == BeamType::kL010z000i_i232
     ||beam == BeamType::kL010z000i_i259
     ||beam == BeamType::kL010z000i_i300
     ||beam == BeamType::kL010z000i_i317  
     ||beam == BeamType::kL010z000i_i326
     ||beam == BeamType::kL010z000i_i380){
      for(int i=0;i<NRUNPERIODS;i++){
	pots.push_back(pot_ndhornoff[i]);
      }
    }
  }
  else{
    if (ismrcc){
      for(int i=0;i<NRUNPERIODS;i++){
	pots.push_back(pot_fdmrcc[i]);
      }
    }
    //FD has no intensity samples
    else if (beam==BeamType::kL010z185i){
      for(int i=0;i<NRUNPERIODS;i++){
	pots.push_back(pot_fd[i]);
      }
    }
    //FD has no intensity samples
    else if (beam==BeamType::kL010z000i){
      MAXMSG("NueStandard",Msg::kWarning,5)<<"Can't use NueStandard::GetRPWBeamWeight for FD L010000 beams, stop it!"<<endl;
      return 0;
    }
  }

  if(nr->mctrue.nuFlavor < -1000){
    MAXMSG("NueStandard",Msg::kWarning,5)<<"No truth info for this event returning 0 for beamweight"<<endl;
    return 0;
  }

  double weight = GetRPWBeamWeight(nr->fluxweights.RPtotbeamweight,pots);
  return weight;
}

double NueStandard::GetRPWBeamWeight(std::vector<double> weights, std::vector<double> pots)
{

  double w=0.;
  double p=0.;
  for(unsigned int i=0;i<weights.size();i++){
    w+=weights[i]*pots[i];
    p+=pots[i];
  }

  if (p==0){
    cout<<"Total pots is set to zero, can not compute a exposure weighted average beam weight.  Returning 1"<<endl;
    return 1;
  }

  return w/p;
}

bool NueStandard::PassesNCCleaningCuts(NueRecord* nr)
{                                                                                
  // 40ns timing cut
  if (TMath::Abs(nr->eventq.minTimeSeparation)<40e-9)
    return false;
                                                                                
  // tiny events are mostly junk
  if (nr->srevent.totalStrips<5)    return false;
                                                                                
  // this cuts very steep showers, leaking in
  if ( (nr->srevent.totalStrips/(nr->srevent.planes*nr->srevent.planes))>1.0)
    return false;
                                                                                
  // this cuts leakage which leaves activity in partially instrumented
  // region
  if ( (nr->eventq.edgeActivityStrips>3)
       && (nr->eventq.edgeActivityPH>1000)
       && (nr->srevent.energyGeV<5.0)
       && (nr->srshower.planes>nr->srtrack.planes) ) // don't cut out CC
    return false;
  // also cut out events with too much activity in the opposite edge region
  if ( (nr->eventq.oppEdgeStrips>3)
       && (nr->eventq.oppEdgePH>1000)
       && (nr->srevent.energyGeV<5.0 )
         && (nr->srshower.planes>nr->srtrack.planes) ) // don't cut out CC
    return false;
                                                                                
  // make additional deltaZ cuts if (40ns<|minDeltaT|<120ns)
  if (TMath::Abs(nr->eventq.closeTimeDeltaZ)<1.0
      && TMath::Abs(nr->eventq.minTimeSeparation)<120e-9)
    return false;
                                                                                
  // Survived all cuts
  return true;
}

void NueStandard::SetDefaultOscParam()
{  
  Double_t dm2_12 = 8.0e-5;  // best fit SNO
  Double_t dm2_23 = 2.32e-3;

  Double_t par[OscPar::kNumParameters] = {};
  par[OscPar::kL] = 735.0;
  par[OscPar::kTh23] = 3.1415926/4.0;
  par[OscPar::kTh12] = 0.59365;    // Sin2(2Th12) = 0.86
  par[OscPar::kTh13] = 0.19885;  // Sin2(2Th13) = 0.15;
  par[OscPar::kDeltaM23] = dm2_23; //normal heirarchy
  par[OscPar::kDeltaM12] = dm2_12;
  par[OscPar::kDensity] = 2.75; //standard rock density
  par[OscPar::kDelta] = 0;
  par[OscPar::kNuAntiNu] = 1;

  // Add NSI parameters:
  par[OscPar::kEps_ee] =  0;
  par[OscPar::kEps_emu] =  0;
  par[OscPar::kEps_etau] =  0;
  par[OscPar::kEps_mumu] =  0;
  par[OscPar::kEps_mutau] = 0; 
  par[OscPar::kEps_tautau] = 0; 
  par[OscPar::kDelta_emu] = 0;
  par[OscPar::kDelta_etau] = 0;
  par[OscPar::kDelta_mutau] = 0;

  par[OscPar::kDm41] =  0;
  par[OscPar::kTh14] = 0; 
  par[OscPar::kTh24] = 0; 
  par[OscPar::kTh34] = 0;
  par[OscPar::kDelta14] = 0;
  par[OscPar::kDelta24] = 0;


  fOscGen.SetOscParam(par);  

}

void NueStandard::SetDefaultOscParamNoNue()
{
  NueStandard::SetDefaultOscParam();
  NueStandard::SetOscParam(OscPar::kTh13, 0);
}

void NueStandard::SetOscParamBestFitANN()
{
  NueStandard::SetDefaultOscParam();
  NueStandard::SetOscParam(OscPar::kTh13, 0.172088);
}

void NueStandard::SetOscParamBestFitkNN()
{ 
  NueStandard::SetDefaultOscParam(); 
  NueStandard::SetOscParam(OscPar::kTh13, 0.145518);
} 

void NueStandard::SetOscNoMatter()
{
   fOscGen.SetOscParam(OscPar::kDensity, 1e-9);
}

void NueStandard::SetOscParam(double *par)
{
   fOscGen.SetOscParam(par);
}

void NueStandard::SetOscParam(OscPar::OscPar_t pos, double val)
{
  fOscGen.SetOscParam(pos, val);
}

double NueStandard::GetOscWeight(int nuFlavor, int nonOsc, double E)
{
  if (isLSND){ return fOscGen.OscillateLSND(nuFlavor,nonOsc,E);}
  else if (isNSI){ return fOscGen.OscillateNSI(nuFlavor, nonOsc, E);}
  else{ return fOscGen.Oscillate(nuFlavor, nonOsc, E);}
}
double NueStandard::GetNSIOscWeight(int nuFlavor,int nonOsc, double E)
{
  return fOscGen.OscillateNSI(nuFlavor,nonOsc,E);
}
double NueStandard::GetLSNDOscWeight(int nuFlavor,int nonOsc, double E)
{
  return fOscGen.OscillateLSND(nuFlavor,nonOsc,E);
}
void NueStandard::FillDefaultOscParam(double* par)
{
  Double_t dm2_12 = 8.0e-5;  // best fit SNO
  Double_t dm2_23 = 2.32e-3;
                                                                                      
  par[OscPar::kL] = 735.0;
  par[OscPar::kTh23] = 3.1415926/4.0;
  par[OscPar::kTh12] = 0.59365;    // Sin2(2Th12) = 0.86
  par[OscPar::kTh13] = 0.19885;  // Sin2(2Th13) = 0.15;
  par[OscPar::kDeltaM23] = dm2_23; //normal heirarchy
  par[OscPar::kDeltaM12] = dm2_12;
  par[OscPar::kDensity] = 2.75; //standard rock density
  par[OscPar::kDelta] = 0;
  par[OscPar::kNuAntiNu] = 1;

  // Add NSI parameters:
  par[OscPar::kEps_ee] =  0;
  par[OscPar::kEps_emu] =  0;
  par[OscPar::kEps_etau] =  0;
  par[OscPar::kEps_mumu] =  0;
  par[OscPar::kEps_mutau] = 0; 
  par[OscPar::kEps_tautau] = 0; 
  par[OscPar::kDelta_emu] = 0;
  par[OscPar::kDelta_etau] = 0;
  par[OscPar::kDelta_mutau] = 0;

  par[OscPar::kDm41] =  0;
  par[OscPar::kTh14] = 0; 
  par[OscPar::kTh24] = 0; 
  par[OscPar::kTh34] = 0;
  par[OscPar::kDelta14] = 0;
  par[OscPar::kDelta24] = 0;

}

void NueStandard::GetOscParam(double* par){
  fOscGen.GetOscParam(par);
}

Bool_t NueStandard::IsRun1(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1145036420);
}


Bool_t NueStandard::IsRun2(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1145036420 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1190028111);
}

Bool_t NueStandard::IsRun3(NueRecord* nr)
{
//True Run3
// return(IsRun3NotPrime(nr));

//Run3Prime
 return(IsRun3Prime(nr));
}

Bool_t NueStandard::IsRun3Prime(NueRecord* nr)
{
  return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1190028111 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1370000000);
}

Bool_t NueStandard::IsRun3NotPrime(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1190028111 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1249594050);
}

Bool_t NueStandard::IsRun4(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1253039760 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1254248760);
}

Bool_t NueStandard::IsRun4RHC(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1254248760 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1269298800);
}

Bool_t NueStandard::IsRun5(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1269298800 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1273028400);
}

Bool_t NueStandard::IsRun6(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1273683600 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1278924900);
}

Bool_t NueStandard::IsRun7RHC(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1278924900 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1299011200);
}

Bool_t NueStandard::IsRun8(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1299011200 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1311993300);
}

Bool_t NueStandard::IsRun9RHC(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1311993300 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1316101520);
}

Bool_t NueStandard::IsRun10(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1316101520 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1370000000);
}


Bool_t NueStandard::IsRun11(NueRecord* nr)
{
 return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1378742400 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1409925600);
}

Bool_t NueStandard::IsRun12(NueRecord* nr)
{
  return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1414210860 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1436227200);
}
Bool_t NueStandard::IsRun13(NueRecord* nr)
{
  return(nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1445616000);
}

//RP XII end time set for July 7, 2015 (probably a few days into beam off)
//RP XIII start time set for Oct 23, 2015 - 4PM UTC (11 AM FNAL local -- roughly one hour before Sam Childress's noon note)

Bool_t NueStandard::IsSpecialRun(NueRecord* nr)
{
  return(TMath::Abs(nr->bmon.hornI) < 1 || (TMath::Abs(nr->bmon.hornI) > 193 && TMath::Abs(nr->bmon.hornI) < 200 && nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() < 1370000000 ) || ( nr->GetHeader().GetVldContext().GetTimeStamp().GetSec() > 1370000000 && TMath::Abs(nr->bmon.hornI) < 193 && TMath::Abs(nr->bmon.hornI) > 130 ) ); 
}

 Double_t NueStandard::GetNDDataWeights(NueRecord *nr){
  //Function to calculate Data event weights to correct for run proportions
   Double_t eventWeight = 1.0;
  
  //Check if it is MC
   if(nr->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kData &&
      nr->GetHeader().GetVldContext().GetDetector() == Detector::kNear)
   {
    
     //OLD ND DATA WEIGHTS FOR 3RD ANALYSIS
     //if(IsRun4(nr)) eventWeight *= 0.65625;
     //else if(IsRun5(nr)) eventWeight *= 0.701299;
     //else if(IsRun6(nr)) eventWeight *= 0.661466;
     //END OF OLD WEIGHTS

  if(!IsSpecialRun(nr)){
     //Period I 
    if(IsRun1(nr)) eventWeight *= 0.639286;
    else if(IsRun2(nr)) eventWeight *= 0.792391;
    else if(IsRun4(nr)) eventWeight *= 0.668901;
    else if(IsRun5(nr)) eventWeight *= 0.720630;
    else if(IsRun6(nr)) eventWeight *= 0.673991;
  
    //Period II
    else if(IsRun7RHC(nr)) eventWeight *= 1.13261095;
    else if(IsRun9RHC(nr)) eventWeight *= 0.91593815;

    //Period III
    else if(IsRun8(nr)) eventWeight *= 0.71409596;
    }
   }
   return(eventWeight);
 }

//This version was for the 5 LEM bin x 1 E bin extrapolation


Double_t NueStandard::GetPredWeights_DO_NOT_USE(NueRecord *nr){
  //Function to return prediction-to-raw ratio for nominal FD MC events
  Double_t eventWeight = 1.0;

  if(nr->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC &&
     nr->GetHeader().GetVldContext().GetDetector() == Detector::kFar &&
     NueStandard::PassesSelection(nr,Selection::kPre)){

    int iaction = nr->mctrue.interactionType;
    int inu = nr->mctrue.nuFlavor;
    int inunoosc = nr->mctrue.nonOscNuFlavor;
    int bgtype = (int)Background::TranslateFromMC(iaction,inu,inunoosc);
    
    double lemval = NueStandard::GetPIDValue(nr,Selection::kMCNN);
    int pidint = 0;
    if (lemval < 0) cout << "Got negative LEM PID value!!" << endl;
    if (lemval >= 0 && lemval < 0.5) pidint = 0;
    if (lemval >= 0.5 && lemval < 0.6) pidint = 1;
    if (lemval >= 0.6 && lemval < 0.7) pidint = 2;
    if (lemval >= 0.7 && lemval < 0.8) pidint = 3;
    if (lemval >= 0.8) pidint = 4;

    if (bgtype < 0 || bgtype > 4) cout << "Unknown background type. No prediction-to-raw weights available!" << endl;
    eventWeight = predwts[bgtype][pidint];

  }

  return eventWeight;
  
}


Double_t NueStandard::GetPredWeights(NueRecord *nr) {
  //Function to return prediction-to-raw ratio for nominal FD MC events
  Double_t eventWeight = 1.0;
  NueConvention::NueEnergyCorrection(nr);

  if(nr->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC &&
     nr->GetHeader().GetVldContext().GetDetector() == Detector::kFar &&
     NueStandard::PassesSelection(nr,Selection::kPre)){

    int iaction = nr->mctrue.interactionType;
    int inu = nr->mctrue.nuFlavor;
    int inunoosc = nr->mctrue.nonOscNuFlavor;
    int bgtype = (int)Background::TranslateFromMC(iaction,inu,inunoosc);
    
    double lemval = NueStandard::GetPIDValue(nr,Selection::kMCNN);
    int pidint = 0;
    if (lemval < 0) cout << "Got negative LEM PID value!! " << lemval << endl;
    if (lemval >= 0 && lemval < 0.5) pidint = 0;
    if (lemval >= 0.5 && lemval < 0.6) pidint = 1;
    if (lemval >= 0.6 && lemval < 0.7) pidint = 2;
    if (lemval >= 0.7 && lemval < 0.8) pidint = 3;
    if (lemval >= 0.8) pidint = 4;

    double recoE = nr->srevent.phNueGeV;
    int recoint = 0;
    if (recoE >= 1 && recoE < 2) recoint = 0;
    if (recoE >= 2 && recoE < 3) recoint = 1;
    if (recoE >= 3 && recoE < 4) recoint = 2;
    if (recoE >= 4 && recoE < 5) recoint = 3;
    if (recoE >= 5 && recoE < 8) recoint = 4;

    if (bgtype < 0 || bgtype > 4) cout << "Unknown background type. No prediction-to-raw weights available!" << endl;

    if (bgtype == 0) eventWeight = predwts_nuecc[recoint][pidint];
    if (bgtype == 1) eventWeight = predwts_nc[recoint][pidint];
    if (bgtype == 2) eventWeight = predwts_numucc[recoint][pidint];
    if (bgtype == 3) eventWeight = predwts_bnuecc[recoint][pidint];
    if (bgtype == 4) eventWeight = predwts_nutaucc[recoint][pidint];
    
  }

  return eventWeight;
  }



Double_t NueStandard::GetMCWeights(NueRecord *nr)
{
  //Function to calculate MC event weights to correct for various effect
  Double_t eventWeight = 1.0;
  
  //Check if it is MC
  if(nr->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC)
   {
     //Correct MC flux (i.e. SKZP and Helium) 
     eventWeight *= NueStandard::GetIntensityBeamWeight(nr);
     eventWeight *= NueStandard::GetSKZPBeamWeight(nr);
   }
  
  return(eventWeight);   
}

Double_t NueStandard::GetIntensityBeamWeight(NueRecord *nr){

   Double_t eventWeight = 1.0;

   if(nr->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC &&
      nr->GetHeader().GetVldContext().GetDetector() == Detector::kNear &&
      IsRun3(nr)){

     switch(nr->GetHeader().GetBeamType()){
       case BeamType::kL010z185i_i124: eventWeight = 1.24932; break;
       case BeamType::kL010z185i_i191: eventWeight = 1.62047; break;
       case BeamType::kL010z185i_i213: eventWeight = 1.19021; break;
       case BeamType::kL010z185i_i224: eventWeight = 1.69376; break;
       case BeamType::kL010z185i_i232: eventWeight = 1.28621; break;
       case BeamType::kL010z185i_i243: eventWeight = 1.26755; break;
       case BeamType::kL010z185i_i257: eventWeight = 1.07360; break;
       case BeamType::kL010z185i_i282: eventWeight = 1.11812; break;
       case BeamType::kL010z185i_i303: eventWeight = 1.06092; break;
       case BeamType::kL010z185i_i324: eventWeight = 2.60657; break;
       default:                                               break;
     }

   }

   return(eventWeight);
 }


 Double_t NueStandard::GetSKZPBeamWeight(NueRecord *nr)
 {
  //One time intialization
   string skzpVersion="DetXs";
   const int release = nr->GetHeader().GetRelease();
   if(ReleaseType::IsDaikon(release) && ReleaseType::GetMCSubVersion(release) == 7) 
   {
    skzpVersion="Dogwood1_Daikon07_v2";
   }
   else if(ReleaseType::IsDaikon(release) && ReleaseType::GetMCSubVersion(release) == 10)
   {
    skzpVersion="MINOSplus_2014_v2";
   }

 
   static SKZPWeightCalculator skzpWC(skzpVersion,true);  //static because we only want this created once
                                                       //and not every time it's called within a loop
   static bool firstTimeFunctionCalled = true;
  
   if(firstTimeFunctionCalled)
   {
    //print out the configuration information once
     skzpWC.PrintReweightConfig(std::cout);
     firstTimeFunctionCalled = false;
   }

  //Now get the weight
   Double_t eventWeight = 1.0;
   //cout << "NS1" << endl;
  //Check if it is MC
   if(nr->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC)
   {

    //Correct MC flux (i.e. SKZP and Helium) 
    double pt = (nr->fluxinfo.tpx * nr->fluxinfo.tpx);
           pt += (nr->fluxinfo.tpy * nr->fluxinfo.tpy);
           pt = sqrt(pt);
	   //cout << "NS2" << endl;
    eventWeight = skzpWC.GetBeamWeight(
                                nr->GetHeader().GetVldContext().GetDetector() /*int det*/,
                                BeamType::ToZarko( nr->GetHeader().GetBeamType() ) /*int Ibeam*/,

                                nr->fluxinfo.tptype /*int tptype*/,
                                pt /*double pt*/,
                                nr->fluxinfo.tpz /*double pz*/,
                      
                                nr->mctrue.nuEnergy /*double true_enu*/,
                                nr->mctrue.nonOscNuFlavor /*int inu*/,  //want weight before neutrino oscillations
                                nr->GetHeader().GetVldContext() /*VldContext vc*/
                               );
   }
   //cout << "NS3" << endl;
   return(eventWeight);   
 }

void NueStandard::ModifyANNPID(NueRecord *nr) {
 
 // This is only a temporary solution for the ANN14 PID.
 
  static bool firstcall = true;
  static TMultiLayerPerceptron* fneuralNet_14inp_daikon04 = 0;

  double ann14pid = -9999.9;

  char *srt_dir = getenv("SRT_PRIVATE_CONTEXT");
  char annfile[1000];

  if (firstcall) {
    srt_dir = getenv("SRT_PRIVATE_CONTEXT");
    sprintf(annfile,"%s/NueAna/data/ann14_400_14_9.root",srt_dir);
    ifstream Test2(annfile);
    if (!Test2){
      srt_dir = getenv("SRT_PUBLIC_CONTEXT");
      sprintf(annfile,"%s/NueAna/data/ann14_400_14_9.root",srt_dir);
      ifstream Test_again2(annfile);
      if (!Test_again2){
        cout<<"Couldn't find ANN object, blame Jiajie"<<endl;
        exit(0);
      }
    }

    static TFile *f = TFile::Open(annfile);
    fneuralNet_14inp_daikon04 = (TMultiLayerPerceptron*) f->Get("mlp");
  }

  if (!fneuralNet_14inp_daikon04) {
    cout << "Couldn't find ANN object, blame Jiajie" << endl;
    exit(0);
  }

  if (firstcall) {
    cout << "Get ann14 value from : " << annfile << endl;
    firstcall = false;
  }

  NueConvention::NueEnergyCorrection(nr); 

  double paras[14];
  paras[0] = nr->shwfit.par_a;
  paras[1] = nr->shwfit.par_b;
  paras[2] = nr->shwfit.uv_rms_9s_2pe_dw;
  paras[3] = nr->fracvars.fract_road;
  paras[4] = nr->shwfit.LongE;
  paras[5] = nr->fracvars.shw_max;
  paras[6] = nr->fracvars.shw_slp;
  paras[7] = nr->srevent.phNueGeV;
  paras[8] = nr->fracvars.dis2stp;
  paras[9] = nr->fracvars.fract_1_plane;
  paras[10] = nr->fracvars.fract_5_planes;
  paras[11] = nr->fracvars.fract_6_counters;
  paras[12] = nr->fracvars.fract_20_counters;
  paras[13] = TMath::Abs(nr->srevent.endPlane - nr->srevent.begPlane);

  int pass = 1;

  if (nr->shwfit.par_a<-1000) pass = 0;
  if (nr->shwfit.par_b<-1000) pass = 0;
  if (nr->shwfit.uv_rms_9s_2pe_dw<-1000) pass = 0;
  if (nr->fracvars.fract_1_plane<-1000) pass = 0;
  if (nr->fracvars.fract_5_planes<-1000) pass = 0;
  if (nr->fracvars.fract_6_counters<-1000) pass = 0;
  if (nr->fracvars.fract_20_counters<-1000) pass = 0;
  if (nr->fracvars.fract_road<-1000) pass = 0;
  if (nr->shwfit.LongE<-1000) pass = 0;
  if (nr->shwfit.LongE>1000) pass = 0;
  if (nr->fracvars.shw_max<-1000) pass = 0;
  if (nr->fracvars.dis2stp<-1000) pass = 0;
  if (nr->fracvars.shw_slp<0) pass = 0;
  if (nr->srevent.endPlane<-1000) pass = 0;
  if (nr->srevent.begPlane<-1000) pass = 0;
  
  if (pass) {
    ann14pid = fneuralNet_14inp_daikon04->Evaluate(0, paras);
  } 

  // set the ann.pid value to be ANN14 
  nr->ann.pid = ann14pid;
}

Double_t NueStandard::Calc4thAnaANN(NueRecord *nr, Selection::Selection_t sel) {

  // This is only a temporary solution for the ANN11 PID for the 4th Analysis.

  static bool firstcall = true;
  static TMultiLayerPerceptron* fneuralNet_11ann4 = 0;

  double ann4pid = -9999.9;
  char annfile[1000];

  if (firstcall) {
    char *srt_dir = getenv("SRT_PRIVATE_CONTEXT");

    switch(sel) {
    case Selection::kANN4FHC:  
      sprintf(annfile,"%s/NueAna/data/ann11_ana4_fhc.root",srt_dir);  break;
    case Selection::kANN4RHC:
      sprintf(annfile,"%s/NueAna/data/ann11_ana4_rhc.root",srt_dir);  break;
    default:      break;
    }
      
    ifstream Test(annfile);
    if (!Test){
      srt_dir = getenv("SRT_PUBLIC_CONTEXT");

      switch(sel) {
      case Selection::kANN4FHC:
	sprintf(annfile,"%s/NueAna/data/ann11_ana4_fhc.root",srt_dir);	break;
      case Selection::kANN4RHC:
	sprintf(annfile,"%s/NueAna/data/ann11_ana4_rhc.root",srt_dir);	break;
      default:	break;
      }

      ifstream Test_again(annfile);
      if (!Test_again){
        cout<<"Couldn't find ANN root file, ask Xinjie"<<endl;
        exit(0);
      }
    }

    static TFile *f = TFile::Open(annfile);
    fneuralNet_11ann4 = (TMultiLayerPerceptron*) f->Get("mlp");
  }

  if (!fneuralNet_11ann4) {
    cout << "Couldn't find ANN object, ask Xinjie" << endl;
    exit(0);
  }

  if (firstcall) {
    cout << "Get ann11 value from : " << annfile << endl;
    firstcall = false;
  }

  double pars[11];

  pars[0] = nr->shwfit.par_a;
  pars[1] = nr->shwfit.par_b;
  pars[2] = nr->shwfit.uv_molrad_peak_9s_2pe_dw;
  pars[3] = nr->shwfit.uv_rms_9s_2pe_dw;
  pars[4] = nr->mstvars.e4w+nr->mstvars.o4w;
  pars[5] = nr->fracvars.fract_2_planes;
  pars[6] = nr->fracvars.fract_4_planes;
  pars[7] = nr->fracvars.fract_6_planes;
  pars[8] = nr->fracvars.fract_8_counters;
  pars[9] = nr->fracvars.fract_road;
  pars[10]= nr->shwfit.LongE;
  bool pass = true;
  if (nr->shwfit.par_a<-1000) pass = false;
  if (nr->shwfit.par_b<-1000) pass = false;
  if (nr->shwfit.uv_molrad_peak_9s_2pe_dw<-1000) pass = false;
  if (nr->shwfit.uv_rms_9s_2pe_dw<-1000) pass = false;
  if (nr->mstvars.e4w<-1000) pass = false;
  if (nr->mstvars.e4w>500) pass = false;
  if (nr->mstvars.o4w<-1000) pass = false;
  if (nr->mstvars.o4w>500) pass = false;
  if (nr->fracvars.fract_2_planes<-1000) pass = false;
  if (nr->fracvars.fract_4_planes<-1000) pass = false;
  if (nr->fracvars.fract_6_planes<-1000) pass = false;
  if (nr->fracvars.fract_8_counters<-1000) pass = false;
  if (nr->fracvars.fract_road<-1000) pass = false;
  if (nr->shwfit.LongE<-1000) pass = false;
  if (nr->shwfit.LongE>1000) pass = false;

  if (pass) {
      ann4pid = fneuralNet_11ann4->Evaluate(0,pars);
  }

  return ann4pid;
}





void NueStandard::SetE50PID(NueRecord *nr) {
 
  NueConvention::NueEnergyCorrection(nr); 
  
  double retval = -9999;
  double secval = -9999;

  //Standard nue 
  LEMAmbyE50N491 setlemambye50n491;
  retval = setlemambye50n491.GetLEMAmbyE50N491(nr);

  //Sterile search
  LEMAmbyE50S491 setlemambye50s491;
  secval = setlemambye50s491.GetLEMAmbyE50S491(nr);
  

  // set the E50 LEM shtuff
  nr->mcnnv.mcnn_var3 = retval;
  nr->mcnnv.mcnn_var4 = secval;
}
