#ifndef NUEFLUXWEIGHTSANA_H
#define NUEFLUXWEIGHTSANA_H

#include "NueAna/NueAnaBase.h"
#include "NueAna/NueFluxWeights.h"
#include "MCReweight/SKZPWeightCalculator.h"

class NtpMCFluxInfo;
class Zbeam;
class Zfluk;
class Kfluk;

class NueFluxWeightsAna: public NueAnaBase
{

 public:
  NueFluxWeightsAna(NueFluxWeights &fw);
  virtual ~NueFluxWeightsAna();

//  void SetBeam(int b) {beam = b;}
  void SetDetector(int d) {det = d;}
  void SetFluxInfo(NtpMCFluxInfo *q) {fi=q;}
  void Analyze(int evt, RecRecordImp<RecCandHeader> *st);
  void SetKfluk(Kfluk *k) {kfluk=k;}
  void SetSKZPCalc(SKZPWeightCalculator *sc, std::string config) {skzpCalc = sc; cfg = config;};

 private:
  NueFluxWeights &fNueFluxWeight;

  int zbeam;
  int det;

  NtpMCFluxInfo *fi;
  Kfluk *kfluk;
  SKZPWeightCalculator *skzpCalc;
  std::string cfg;
};


#endif //NUEFLUXWEIGHTSANA_H
