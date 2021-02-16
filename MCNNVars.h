#ifndef MCNNVARS_H
#define MCNNVARS_H

#include "TObject.h"
#include <vector>
#include "NueAna/MCNNBestMatch.h"

class TClonesArray;

class MCNNVars: public TObject
{
 public:
  MCNNVars();
  MCNNVars(const MCNNVars *nuefw);
  virtual ~MCNNVars();

  virtual void Print(Option_t* option="") const;
  void Reset();
  void Clear(Option_t* option = "");

  Int_t bestmatches; // <-- number of best matches in bmatch array of MCNNBestMatch
  Float_t meanU; // <--charge weighted mean U
  Float_t meanV; //<---charge weighted mean V
  Int_t meanPlane; //<--- charge weighted mean Plane
  Float_t meanfracQmatched; // <-- mean fractional charge thatwas
                            // matched; tells us if event was properly matched to libraries or not.
  Bool_t mcpresel; // <-- pass library pre-selection?
  Float_t fracCC;  // <-- fraction of best 20 matches that were nue
  Float_t fracCCy; // <-- fraction of best 20 matches that were nue with y<0.5
  Float_t ymean;   // <--mean y of nue matches (among 20 best matches).
  Float_t ncymean;
  Float_t ncmeanfracQmatched;

  Int_t qtot; //<-- total charge of event (PEs) as calculated in the mcnn.
  Float_t mcnn_var1; //<--- open variable.
  Float_t mcnn_var2; //<-- open variable.
  Float_t mcnn_var3; //<--- open variable.
  Float_t mcnn_var4; // <--open variable.

  TClonesArray* bmatch; //-> array of best matches

  private:
   static TClonesArray *fgBmatch;
   ClassDef(MCNNVars,2)
};

#endif //MCNNVARS_H
