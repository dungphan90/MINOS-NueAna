#include "NueAna/FracVar.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(FracVar)

FracVar::FracVar():
  fract_1_plane(ANtpDefVal::kFloat),
  fract_2_planes(ANtpDefVal::kFloat),
  fract_3_planes(ANtpDefVal::kFloat),
  fract_4_planes(ANtpDefVal::kFloat),
  fract_5_planes(ANtpDefVal::kFloat),
  fract_6_planes(ANtpDefVal::kFloat),
  fract_2_counters(ANtpDefVal::kFloat),
  fract_4_counters(ANtpDefVal::kFloat),
  fract_6_counters(ANtpDefVal::kFloat),
  fract_8_counters(ANtpDefVal::kFloat),
  fract_10_counters(ANtpDefVal::kFloat),
  fract_12_counters(ANtpDefVal::kFloat),
  fract_14_counters(ANtpDefVal::kFloat),
  fract_16_counters(ANtpDefVal::kFloat),
  fract_18_counters(ANtpDefVal::kFloat),
  fract_20_counters(ANtpDefVal::kFloat),
  fract_road(ANtpDefVal::kFloat),
  fract_asym(ANtpDefVal::kFloat),
  shw_max(ANtpDefVal::kInt),
  dis2stp(ANtpDefVal::kInt),
  shw_nstp(ANtpDefVal::kInt),
  shw_npl(ANtpDefVal::kInt),
  shw_slp(ANtpDefVal::kFloat),
  vtxph(ANtpDefVal::kFloat),
  passcuts(0),//0 is good
  pid(ANtpDefVal::kFloat),
  pid1(ANtpDefVal::kFloat)
{}

FracVar::~FracVar(){}

void FracVar::Reset()
{
  fract_1_plane = ANtpDefVal::kFloat;
  fract_2_planes = ANtpDefVal::kFloat;
  fract_3_planes = ANtpDefVal::kFloat;
  fract_4_planes = ANtpDefVal::kFloat;
  fract_5_planes = ANtpDefVal::kFloat;
  fract_6_planes = ANtpDefVal::kFloat;
  fract_2_counters = ANtpDefVal::kFloat;
  fract_4_counters = ANtpDefVal::kFloat;
  fract_6_counters = ANtpDefVal::kFloat;
  fract_8_counters = ANtpDefVal::kFloat;
  fract_10_counters = ANtpDefVal::kFloat;
  fract_12_counters = ANtpDefVal::kFloat;
  fract_14_counters = ANtpDefVal::kFloat;
  fract_16_counters = ANtpDefVal::kFloat;
  fract_18_counters = ANtpDefVal::kFloat;
  fract_20_counters = ANtpDefVal::kFloat;
  fract_road = ANtpDefVal::kFloat;
  fract_asym = ANtpDefVal::kFloat;
  shw_max = ANtpDefVal::kInt;
  dis2stp = ANtpDefVal::kInt;
  shw_nstp = ANtpDefVal::kInt;
  shw_npl = ANtpDefVal::kInt;
  shw_slp = ANtpDefVal::kFloat;
  vtxph   = ANtpDefVal::kFloat;
  passcuts = 0;
  pid = ANtpDefVal::kFloat;
  pid1 = ANtpDefVal::kFloat;
}
