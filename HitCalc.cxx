/// $Id: HitCalc.cxx,v 1.3 2005/05/17 20:53:30 boehm Exp $
///
/// class HitCalc
///
/// NueAna package
///
/// Purpose: Calculate 3D Hits, Hit angular clustering and 
///          Tufts Analysis variables.
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Tue Mar 29 2005

#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "NueAna/HitCalc.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(HitCalc)

HitCalc::HitCalc():
    
    fHitTotalEnergy(ANtpDefVal::kFloat),
    fHitTransEnergy(ANtpDefVal::kFloat),
    fHitLongEnergy(ANtpDefVal::kFloat),
    fHitTransCMEnergy(ANtpDefVal::kFloat),
    fHitTransEnergyRatio(ANtpDefVal::kFloat),
    fHitLongEnergyRatio(ANtpDefVal::kFloat),
    fHitTransLongEnergyRatio(ANtpDefVal::kFloat),
    fHitTransCMEnergyRatio(ANtpDefVal::kFloat),
    fHitFarMomBalance(ANtpDefVal::kFloat),
    fHitPeakMomBalance(ANtpDefVal::kFloat),
    fHitFarAngle(ANtpDefVal::kFloat),
    fHitPeakAngle(ANtpDefVal::kFloat)

{}

HitCalc::~HitCalc()
{

}

// Note: the methods below should be kept in alphabetical order.

//......................................................................

void HitCalc::Draw(Option_t */*option*/) 
{

/// To be filled later.

}

void HitCalc::Print(Option_t */*option*/) const
{
    std::cout<<"There's a lot of stuff to print here!"<<std::endl;
}

void HitCalc::Reset()
{
    fHitTotalEnergy=ANtpDefVal::kFloat;
    fHitTransEnergy=ANtpDefVal::kFloat;
    fHitLongEnergy=ANtpDefVal::kFloat;
    fHitTransCMEnergy=ANtpDefVal::kFloat;
    fHitTransEnergyRatio=ANtpDefVal::kFloat;
    fHitLongEnergyRatio=ANtpDefVal::kFloat;
    fHitTransLongEnergyRatio=ANtpDefVal::kFloat;
    fHitTransCMEnergyRatio=ANtpDefVal::kFloat;
    fHitFarMomBalance=ANtpDefVal::kFloat;
    fHitPeakMomBalance=ANtpDefVal::kFloat;
    fHitFarAngle=ANtpDefVal::kFloat;
    fHitPeakAngle=ANtpDefVal::kFloat;

}
