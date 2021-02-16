#ifndef MATRIXHISTS_H
#define MATRIXHISTS_H
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"

class MatrixHists
{

 public:

  MatrixHists(std::string);
  ~MatrixHists();

  TDirectory *fDirectory;

  //Helper histos:
  TH2D *fRecoVsTrueEnergy_ND; // conversion from reco->true 
  TH2D *fRecoVsTrueEnergy_FD; // for ND/FD
  TH2D *fFDVsNDMatrix;        // ND->FD number of events conversion
  TH2D *fFDVsNDMatrixRW;      // ND->FD number of events conversion+new nearwt
  TH2D *fFDVsNDMatrixXSec;    // ND->FD number of events conversion with xsec folded in
  TH2D *fFDVsNDMatrixXSecRW;  // ND->FD number of events conversion with xsec + new nearwt
  TH1D *fEfficiency_ND;       // FD and ND efficiency (with true energy)
  TH1D *fEfficiency_FD;
  TH1D *fPurity_ND;           // FD and ND purity (with reco energy)
  TH1D *fPurity_FD;
  TH1D *fXSec_CC;             // numu CC cross-section with energy
  TH1D *fFracErrOnPred;       // fractional error on predicted FD spectrum with energy
      
  //Check list histos:
  TH1D *fRecoEnergyAllEvents_ND;
  TH1D *fRecoEnergyCCOnlyEvents_ND;
  TH1D *fTrueEnergyCCOnlyEvents_ND;
  TH1D *fTrueEnergyTrueCCFidEvents_ND;
  TH1D *fTrueEnergyNuFlux_ND;
  TH1D *fTrueEnergyNuFluxRW_ND;
  TH1D *fTrueEnergyNuFlux_FD;
  TH1D *fTrueEnergyNuFluxRW_FD;
  TH1D *fTrueEnergyCCFlux_ND;
  TH1D *fTrueEnergyCCFluxRW_ND;
  TH1D *fTrueEnergyCCFlux_FD;
  TH1D *fTrueEnergyCCFluxRW_FD;
  TH1D *fTrueEnergyTrueCCFidEvents_FD;
  TH1D *fTrueEnergyCCOnlyEvents_FD;
  TH1D *fRecoEnergyCCOnlyEvents_FD;
  TH1D *fRecoEnergyAllEvents_FD;

};
#endif //MATRIXHISTS_H
