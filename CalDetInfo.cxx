#include "NueAna/CalDetInfo.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(CalDetInfo)

CalDetInfo::CalDetInfo():
  beamp(ANtpDefVal::kDouble),
  inct(ANtpDefVal::kInt),
  pid(ANtpDefVal::kInt),
  olchi2(ANtpDefVal::kDouble),
  p0stripmaxmip(ANtpDefVal::kInt)
{}

CalDetInfo::~CalDetInfo(){}

void CalDetInfo::Zero()
{
  beamp = 0.;
  inct = 0;
  pid = 0;
  olchi2 = 0.;
  p0stripmaxmip = 0;
}

void CalDetInfo::Reset()
{
  beamp = ANtpDefVal::kDouble;
  inct = ANtpDefVal::kInt;
  pid = ANtpDefVal::kInt;
  olchi2 = ANtpDefVal::kDouble;
  p0stripmaxmip = ANtpDefVal::kInt;
}
