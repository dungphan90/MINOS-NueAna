#ifndef FNHISTS_H
#define FNHISTS_H
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"

class FNHists
{

 public:

  FNHists(std::string);
  ~FNHists();

  TDirectory *fDirectory;

  TH1D *fND_RecoEnergy;
  TH1D *fFD_RecoEnergy;
  TH1D *fND_TrueEnergy;
  TH1D *fFD_TrueEnergy;
      
};
#endif //FNHISTS_H
