////////////////////////////////////////////////////////////////////////
// $Id: NueConvention.h,v 1.16 2015/02/05 20:36:24 wingmc Exp $
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
#ifndef  NUECONVENTIONS
#define NUECONVENTIONS
                                                                                
// This file contains a lot of the standards used in the NueAnalysis
#include "TObject.h"
#include <iostream>
                                                                                
class NtpMCTruth;
class ANtpTruthInfoBeam;
class ANtpTruthInfoBeamNue;
                       
class NueRecord;

namespace NueConvention
{

    // Old ClassType code
    static const Int_t NC = 0;
    static const Int_t numu = 1;
    static const Int_t nue = 2;
    static const Int_t nutau = 3;
    static const Int_t bnue = 4;                                                                                
    Int_t DetermineClassType(Int_t inu, Int_t inunoosc, Int_t iaction);

    //***********************************************
    typedef enum ENueRelease {
      kUnknown  = 0x00,
      kEnt      = 0x01,
      kFirebird = 0x02,
      kGriffin  = 0x03,
      kHydra    = 0x04,    //Production began December 2006
      kImp      = 0x05,
      kJaberwocky = 0x06,
      kKraken   = 0x07,
      kLeviathan = 0x08,
      kManticore = 0x09
    } NueRel_t;

   
    //*********************************************************
    //Fiducial Volume Code
                                                                                
    int IsInsideNearFiducial_Nue_Extended(float x, float y, float z);
    int IsInsideFarFiducial_Nue_Extended(float x, float y, float z);
                                                                                
    int IsInsideNearFiducial_Nue_Standard(float x, float y, float z, bool isMC);
    int IsInsideFarFiducial_Nue_Standard(float x, float y, float z, bool isMC);
    int IsInsideNearFiducial_MRE_Standard(float x, float y, float z, bool isMC);


    //Special Case Functions 
    Int_t InPartialRegion(UShort_t plane, UShort_t strip);




    // **********************************************************
    // Oscillation Functions       
    float Oscillate(NtpMCTruth *mcth,float L, float dm2, float theta23, float UE32);
    float Oscillate(ANtpTruthInfoBeam *ib,float L, float dm2, float theta23, float UE32);
    float Oscillate(int nuFlavor, int nonOscNuFlavor, float Energy,
                  float L, float dm2, float theta23, float U);
    float Oscillate(ANtpTruthInfoBeamNue *ib);

    // **********************************************************
    // Energy Correction Functions        

    void NueEnergyCorrection(NueRecord* nr, bool isLinearityFixMC = false);
    float NueEnergyCorrection(float meu, int type, bool isMC, int detector, bool isLinearityFixMC);

    void RHCNueEnergyCorrection(NueRecord* nr, bool isLinearityFixMC = false);
    float RHCNueEnergyCorrection(float meu, int type, bool isMC, int detector, bool isLinearityFixMC = false);

    void MINOSPLUSNueEnergyCorrection(NueRecord* nr, bool isLinearityFixMC = false);
    float MINOSPLUSNueEnergyCorrection(float meu, int type, bool isMC, int detector, bool isLinearityFixMC = false);

    void NueEnergyCorrectionNeverUseThisFunction(NueRecord* nr);
    float NueEnergyCorrectionNeverUseThisFunction(float meu, int type, bool isMC, int detector);

}

//Keep things compliant to the old standard
namespace ClassType = NueConvention;

// Code For File Position Reading
struct FilePosition
{
   int Run;
   int SubRun;
   int Snarl;
   int Event;
};
                                                                                
bool operator>(FilePosition one, FilePosition two);
bool operator<(FilePosition one, FilePosition two);
bool operator==(FilePosition one, FilePosition two);

#endif
