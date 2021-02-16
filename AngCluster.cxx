/// $Id: AngCluster.cxx,v 1.4 2007/05/18 21:01:09 danche Exp $
///
/// class AngCluster
///
/// NueAna package
///
/// Purpose: Calculate 3D Hits, Hit angular clustering and
///          Tufts Analysis variables.
///
/// Author: Alex Sousa <asousa@minos.phy.tufts.edu>
///
/// Created on: Thu April 28 2005

#include <iostream>
#include "TH1F.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "NueAna/AngCluster.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(AngCluster)

AngCluster::AngCluster():
    fACluRmsShwAxis(ANtpDefVal::kFloat),
    fACluRmsZAxis(ANtpDefVal::kFloat),
    fACluPrimEnergy(ANtpDefVal::kFloat),
    fACluPrimEnergyRatio(ANtpDefVal::kFloat),
    fACluShwDirX(ANtpDefVal::kFloat),
    fACluShwDirY(ANtpDefVal::kFloat),
    fACluShwDirZ(ANtpDefVal::kFloat),
    weightedPH0(ANtpDefVal::kFloat),
    weightedPH1(ANtpDefVal::kFloat),
    weightedPH2(ANtpDefVal::kFloat),
    weightedPH3(ANtpDefVal::kFloat)

{}

AngCluster::~AngCluster()
{}

void AngCluster::Draw(Option_t */*option*/)
{

/// To be filled later.

}

void AngCluster::Print(Option_t */*option*/) const
{
    std::cout<<"There's a lot of stuff to print here!"<<std::endl;
}


void AngCluster::Reset()
{

    fACluRmsShwAxis=ANtpDefVal::kFloat;
    fACluRmsZAxis=ANtpDefVal::kFloat;
    fACluPrimEnergy=ANtpDefVal::kFloat;
    fACluPrimEnergyRatio=ANtpDefVal::kFloat;
    fACluShwDirX=ANtpDefVal::kFloat;
    fACluShwDirY=ANtpDefVal::kFloat;
    fACluShwDirZ=ANtpDefVal::kFloat;
    weightedPH0=ANtpDefVal::kFloat;
    weightedPH1=ANtpDefVal::kFloat;
    weightedPH2=ANtpDefVal::kFloat;
    weightedPH3=ANtpDefVal::kFloat;

}
