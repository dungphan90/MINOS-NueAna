#ifndef TIMINGVARS_H
#define TIMINGVARS_H

#include "TObject.h"

class TH1F;
class TF1;
class TPad;
class TPaveText;

class TimingVars : public TObject
{

public:
   TimingVars();
   TimingVars(const TimingVars &tv);
   virtual ~TimingVars();


   void Reset();

   //TimingVars variables

   Int_t ShwMaxBin;
   Int_t TrkMaxBin;
   Int_t TotalMaxBin;
   Float_t ShwMaxBinCont;
   Float_t TrkMaxBinCont;
   Float_t TotalMaxBinCont;
   Float_t ShwSecondMaxBin;
   Float_t TrkSecondMaxBin;
   Float_t TotalSecondMaxBin;
   Float_t ShwEventPulseHeight;
   Float_t TrkEventPulseHeight;
   Float_t TotalEventPulseHeight;
   Float_t ShwPercentofTotal;
   Float_t TrkPercentofTotal;
   Float_t TotalPercentofTotal;
   Int_t ShwbinsPrior;
   Int_t TrkbinsPrior;
   Int_t TotalbinsPrior;
   Int_t ShwBinsOver40;
   Int_t TrkBinsOver40;
   Int_t TotalBinsOver40;

    TH1F *lenepl; //!  don't write this to the tree
    TH1F *tenestu; //!  don't write this to the tree
    TH1F *tenestv; //!  don't write this to the tree

    TF1 *efit; //!  don't write this to the tree
    TF1 *hfit; //!  don't write this to the tree

    TPaveText *info1; //! don't write this to the tree

private:

ClassDef(TimingVars,3)
};

#endif// TIMINGVARS_H
