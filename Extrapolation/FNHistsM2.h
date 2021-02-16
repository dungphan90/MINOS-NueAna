#ifndef FNHISTSM2_H
#define FNHISTSM2_H
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"

class FNHistsM2
{

 public:

  FNHistsM2(std::string);
  ~FNHistsM2();

  TDirectory *fDirectory;

  TH1D *fND_RecoEnergy;
  TH1D *fFD_RecoEnergy;
  TH1D *fND_TrueEnergy;
  TH1D *fFD_TrueEnergy;

  TH2D *fFD_RecoVsTrue;

  TH1D* fFD_numu_TrueEnergyBase;
  TH2D* fFD_RecoVsTrue_BNue;

};
#endif //FNHISTSM2_H
