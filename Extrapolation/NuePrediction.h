#ifndef NUEPREDICTION_H
#define NUEPREDICTION_H
#include <vector>
#include <map>
#include "NueAna/Extrapolation/Extrapolation.h"
#include "NueAna/Extrapolation/NueBackground.h"
#include "NueAna/Extrapolation/NueExtrapolation.h"

class NuePrediction
{
  
 public:
  
  NuePrediction();
  ~NuePrediction();
  Bool_t AddExtrapolation(NueExtrapolation *,Extrapolation::Extrapolation_t);
  Bool_t AddBackground(NueBackground *,Extrapolation::Extrapolation_t);
  TH1D *GetPrediction(Detector::Detector_t,Double_t);
  void Clear();
  void ClearPred();
  void ClearExtrap();
  void Draw();

 private:

  TH1D *fFinalPred;
  std::vector<TH1*> fCurrentPred;
  std::vector<NueBackground*> fBgVec;
  std::vector<Extrapolation::Extrapolation_t> fExtrapVec;
  std::map<Extrapolation::Extrapolation_t,NueExtrapolation*> fExtrapMap;

};
#endif //NUEPREDICTION_H
