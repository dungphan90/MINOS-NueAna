#include "NueAna/Extrapolation/FNHists.h"

FNHists::FNHists(std::string name) :
  fND_RecoEnergy(0),fFD_RecoEnergy(0),
  fND_TrueEnergy(0),fFD_TrueEnergy(0)
{
  fDirectory = new TDirectory(name.c_str(),name.c_str());
}

FNHists::~FNHists()
{
}
