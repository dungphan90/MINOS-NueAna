#ifndef MCNNBESTMATCH_H
#define MCNNBESTMATCH_H

#include "TObject.h"
#include <vector>

class MCNNBestMatch;
class TClonesArray;

std::ostream &operator << (std::ostream& os, const MCNNBestMatch& mcnnbmatch);

class MCNNBestMatch: public TObject
{
 public: 
  MCNNBestMatch();
  MCNNBestMatch(const MCNNBestMatch *nuefw);
  virtual ~MCNNBestMatch();

  virtual void Print(Option_t* option="") const;
  void Reset();
  void Clear(Option_t* option = "");

  Int_t run; //<-- MC run number
  Int_t snarl; //<-- MC snarl number
  Int_t interactionType; //<-- 0==NC; 1==CC
  Int_t nuFlavor; //<-- true flavor of neutrino best match
  Float_t true_enu; //<-- true energy of neutrino best match
  Float_t true_y; //<-- true y of neutrino best match
  Int_t resonanceCode; //<-- resonance code of neutrino best match
  Float_t dlnL; //<-- delta log(L) with input event
  Float_t fracQmatched; //<-- mean fractional charge matched
  Int_t stdhepsize;//<-- size of stdhep array

  TClonesArray *stdhep;//<-- stdhep array

 private:
  ClassDef(MCNNBestMatch,1);


};

#endif //MCNNBESTMATCH_H
