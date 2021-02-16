#ifndef NUEFLUXWEIGHTS_H
#define NUEFLUXWEIGHTS_H

#include "TObject.h"
#include <vector>

class NueFluxWeights: public TObject
{
 public:
  NueFluxWeights();
  NueFluxWeights(const NueFluxWeights *nuefw);
  virtual ~NueFluxWeights();

  virtual void Print(Option_t* option="") const;
  void Reset();
  void Clear(Option_t* option = "");

  double totbeamweight;
  double detectorWeight;
  double kflukweight;
  double totskzpweight;

  double skzpTrkEnergy;
  double skzpShwEnergy;
  std::string skzpConfig;

  std::vector<double> RPtotbeamweight;

  private:
   ClassDef(NueFluxWeights,4)
};

#endif //NUEFLUXWEIGHTS_H
