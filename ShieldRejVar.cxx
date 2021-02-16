#include "NueAna/ShieldRejVar.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(ShieldRejVar)

  ShieldRejVar::ShieldRejVar():ShieldHit(ANtpDefVal::kInt){

}

ShieldRejVar::~ShieldRejVar(){}

void ShieldRejVar::Reset()
{
  ShieldHit = ANtpDefVal::kInt;

}
