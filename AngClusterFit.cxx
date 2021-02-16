/// $Id: AngClusterFit.cxx,v 1.5 2005/06/02 00:17:00 asousa Exp $
///
/// class AngClusterFit
///
/// NueAna package
///
/// Purpose: Return 3D Hit angular clustering shower fitting variables.
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Fri May 06 2005

#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "NueAna/AngClusterFit.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"


ClassImp(AngClusterFit)

AngClusterFit::AngClusterFit():
    fACluFitParA(ANtpDefVal::kFloat),
    fACluFitParB(ANtpDefVal::kFloat),
    fACluFitParLongE0(ANtpDefVal::kFloat),
    fACluFitShwMax(ANtpDefVal::kFloat),
    fACluFitE0EnergyRatio(ANtpDefVal::kFloat),    
    fACluFitParL1(ANtpDefVal::kFloat),
    fACluFitParL2(ANtpDefVal::kFloat),
    fACluFitParC12(ANtpDefVal::kFloat),
    fACluFitParTransE0(ANtpDefVal::kFloat),
    fACluFitLongChiSq(ANtpDefVal::kFloat),
    fACluFitLongConv(ANtpDefVal::kInt),
    fACluFitLongNDF(ANtpDefVal::kFloat),
    fACluFitTransChiSq(ANtpDefVal::kFloat),  
    fACluFitTransConv(ANtpDefVal::kInt), 
    fACluFitTransNDF(ANtpDefVal::kFloat),
    fACluFitAsymPeak(ANtpDefVal::kFloat),
    fACluFitAsymVert(ANtpDefVal::kFloat),
    fACluFitMolRadPeak(ANtpDefVal::kFloat),
    fACluFitMolRadVert(ANtpDefVal::kFloat),
    fACluFitMean(ANtpDefVal::kFloat),
    fACluFitRMS(ANtpDefVal::kFloat),
    fACluFitSkew(ANtpDefVal::kFloat),
    fACluFitKurt(ANtpDefVal::kFloat)
{}



AngClusterFit::~AngClusterFit() 
{}


void AngClusterFit::Draw(Option_t */*option*/) 
{

/// To be filled later.

}

void AngClusterFit::Print(Option_t */*option*/) const
{
    std::cout<<"There's a lot of stuff to print here!"<<std::endl;
}

void AngClusterFit::Reset()
{
    fACluFitParA=ANtpDefVal::kFloat;
    fACluFitParB=ANtpDefVal::kFloat;
    fACluFitParLongE0=ANtpDefVal::kFloat;
    fACluFitShwMax=ANtpDefVal::kFloat;
    fACluFitE0EnergyRatio=ANtpDefVal::kFloat;
    
    fACluFitParL1=ANtpDefVal::kFloat;
    fACluFitParL2=ANtpDefVal::kFloat;
    fACluFitParC12=ANtpDefVal::kFloat;
    fACluFitParTransE0=ANtpDefVal::kFloat;

    fACluFitLongChiSq=ANtpDefVal::kFloat;
    fACluFitLongConv=ANtpDefVal::kInt;
    fACluFitLongNDF=ANtpDefVal::kFloat;
  
    fACluFitTransChiSq=ANtpDefVal::kFloat;
    fACluFitTransConv=ANtpDefVal::kInt;    
    fACluFitTransNDF=ANtpDefVal::kFloat;

    fACluFitAsymPeak=ANtpDefVal::kFloat;
    fACluFitAsymVert=ANtpDefVal::kFloat;
    fACluFitMolRadPeak=ANtpDefVal::kFloat;
    fACluFitMolRadVert=ANtpDefVal::kFloat;
    fACluFitMean=ANtpDefVal::kFloat;
    fACluFitRMS=ANtpDefVal::kFloat;
    fACluFitSkew=ANtpDefVal::kFloat;
    fACluFitKurt=ANtpDefVal::kFloat;


}
