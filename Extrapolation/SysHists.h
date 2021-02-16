#ifndef SYSHISTS_H
#define SYSHISTS_H
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"

class SysHists
{

 public:

  SysHists(std::string);
  ~SysHists();

  TDirectory *fDirectory;

  TH1D *fND_RecoEnergy;
  TH1D *fFD_RecoEnergy;

  TH2D* fND_RecoVsTrue;
  TH2D *fFD_RecoVsTrue;

};
#endif //SYSHISTS_H
