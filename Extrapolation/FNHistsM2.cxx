#include "NueAna/Extrapolation/FNHistsM2.h"

FNHistsM2::FNHistsM2(std::string name) :
  fND_RecoEnergy(0),fFD_RecoEnergy(0),
  fND_TrueEnergy(0),fFD_TrueEnergy(0),
  fFD_RecoVsTrue(0),
  fFD_numu_TrueEnergyBase(0), fFD_RecoVsTrue_BNue(0)
{
  fDirectory = new TDirectory(name.c_str(),name.c_str());
}

FNHistsM2::~FNHistsM2()
{
}
