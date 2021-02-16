#include "NueAna/Extrapolation/MatrixHists.h"

MatrixHists::MatrixHists(std::string name) :
  fRecoVsTrueEnergy_ND(0),fRecoVsTrueEnergy_FD(0),
  fFDVsNDMatrix(0),fFDVsNDMatrixRW(0),
  fFDVsNDMatrixXSec(0),fFDVsNDMatrixXSecRW(0),
  fEfficiency_ND(0),fEfficiency_FD(0),fPurity_ND(0),
  fPurity_FD(0),fXSec_CC(0),fFracErrOnPred(0),
  fRecoEnergyAllEvents_ND(0),fRecoEnergyCCOnlyEvents_ND(0),
  fTrueEnergyCCOnlyEvents_ND(0),fTrueEnergyTrueCCFidEvents_ND(0),
  fTrueEnergyNuFlux_ND(0),fTrueEnergyNuFluxRW_ND(0),
  fTrueEnergyNuFlux_FD(0),fTrueEnergyNuFluxRW_FD(0),
  fTrueEnergyCCFlux_ND(0),fTrueEnergyCCFluxRW_ND(0),
  fTrueEnergyCCFlux_FD(0),fTrueEnergyCCFluxRW_FD(0),
  fTrueEnergyTrueCCFidEvents_FD(0),fTrueEnergyCCOnlyEvents_FD(0),
  fRecoEnergyCCOnlyEvents_FD(0),fRecoEnergyAllEvents_FD(0)
{
  fDirectory = new TDirectory(name.c_str(),name.c_str());
}

MatrixHists::~MatrixHists()
{
}
