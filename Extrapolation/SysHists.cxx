#include "NueAna/Extrapolation/SysHists.h"

SysHists::SysHists(std::string name) :
  fND_RecoEnergy(0),fFD_RecoEnergy(0),
  fND_RecoVsTrue(0),
  fFD_RecoVsTrue(0)
{
  fDirectory = new TDirectory(name.c_str(),name.c_str());
}

SysHists::~SysHists()
{
}
