#include <iostream>
#include "TClonesArray.h"
#include "NueAna/MCNNVars.h"
#include "NueAna/MCNNBestMatch.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(MCNNVars)

TClonesArray *MCNNVars::fgBmatch = 0;


MCNNVars::MCNNVars() : 
    bmatch(0)
{
  Reset(); 
}

MCNNVars::MCNNVars(const MCNNVars *mv):
  meanU(mv->meanU),
  meanV(mv->meanV),
  meanPlane(mv->meanPlane),
  meanfracQmatched(mv->meanfracQmatched),
  mcpresel(mv->mcpresel),
  fracCC(mv->fracCC),
  fracCCy(mv->fracCCy),
  ymean(mv->ymean),
  ncymean(mv->ncymean),
  ncmeanfracQmatched(mv->ncmeanfracQmatched),
  qtot(mv->qtot),
  mcnn_var1(mv->mcnn_var1),
  mcnn_var2(mv->mcnn_var2),
  mcnn_var3(mv->mcnn_var3),
  mcnn_var4(mv->mcnn_var4),
  bmatch(mv->bmatch)
{
}

MCNNVars::~MCNNVars()
{
  //Clear();
  //if(bmatch) delete bmatch;
}

void MCNNVars::Clear(Option_t* /* option */)
{
  if(bmatch) {bmatch->Clear("C");}
}

void MCNNVars::Reset()
{
  bestmatches = 0; //<-- default value should be 0 so nobody tries to read in an empty TClonesArray of best matches
  meanU = ANtpDefVal::kFloat; // <--charge weighted mean U
  meanV = ANtpDefVal::kFloat; //<---charge weighted mean V
  meanPlane= ANtpDefVal::kInt; //<--- charge weighted mean Plane
  meanfracQmatched= ANtpDefVal::kFloat; // <-- mean fractional charge thatwas
                            // matched; tells us if event was properly matched to libraries or not.
  mcpresel= ANtpDefVal::kBool; // <-- pass library pre-selection?
  fracCC= ANtpDefVal::kFloat;  // <-- fraction of best 20 matches that were nue
  fracCCy= ANtpDefVal::kFloat; // <-- fraction of best 20 matches that were nue with y<0.5
  ymean= ANtpDefVal::kFloat;   // <--mean y of nue matches (among 20 best matches).
  ncymean=  ANtpDefVal::kFloat; // <-- mean y of NC matches
  ncmeanfracQmatched= ANtpDefVal::kFloat; // <-- mean frac charge of NC matches

  qtot= ANtpDefVal::kInt; //<-- total charge of event (PEs) as calculated in the mcnn.
  mcnn_var1= ANtpDefVal::kFloat; //<--- open variable.
  mcnn_var2= ANtpDefVal::kFloat; //<-- open variable.
  mcnn_var3= ANtpDefVal::kFloat; //<--- open variable.
  mcnn_var4= ANtpDefVal::kFloat; //<--- open variable.

  if (!fgBmatch) fgBmatch = new TClonesArray("MCNNBestMatch");
  bmatch = fgBmatch;
  Clear();
}

void MCNNVars::Print(Option_t * /*option*/) const
{
}



