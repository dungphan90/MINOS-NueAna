#include "NueAna/HighHitVars.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(HighHitVars)

HighHitVars::HighHitVars():

  // Pulse height of highest hit in shower in sigcor (and 2nd highest, and 3rd highest, 4rth, 5th, 6th)
  high_hit_1_ph(ANtpDefVal::kFloat),
  high_hit_2_ph(ANtpDefVal::kFloat),
  high_hit_3_ph(ANtpDefVal::kFloat),
  high_hit_4_ph(ANtpDefVal::kFloat),
  high_hit_5_ph(ANtpDefVal::kFloat),
  high_hit_6_ph(ANtpDefVal::kFloat),

  // Pulse height of highest hit in shower in MIP (and 2nd highest, and 3rd highest, 4rth, 5th, 6th)
  high_hit_1_mip(ANtpDefVal::kFloat),
  high_hit_2_mip(ANtpDefVal::kFloat),
  high_hit_3_mip(ANtpDefVal::kFloat),
  high_hit_4_mip(ANtpDefVal::kFloat),
  high_hit_5_mip(ANtpDefVal::kFloat),
  high_hit_6_mip(ANtpDefVal::kFloat),


  //Variables added at Jenny's request - find all the strips above a certain threshold and count them or total their ph:

  hitsabove1MIP_total_ph(ANtpDefVal::kFloat),
  hitsabove2MIP_total_ph(ANtpDefVal::kFloat),
  hitsabove3MIP_total_ph(ANtpDefVal::kFloat),
  hitsabove4MIP_total_ph(ANtpDefVal::kFloat),
  hitsabove5MIP_total_ph(ANtpDefVal::kFloat),
  hitsabove6MIP_total_ph(ANtpDefVal::kFloat),

  hitsabove1MIP_total_mip(ANtpDefVal::kFloat),
  hitsabove2MIP_total_mip(ANtpDefVal::kFloat),
  hitsabove3MIP_total_mip(ANtpDefVal::kFloat),
  hitsabove4MIP_total_mip(ANtpDefVal::kFloat),
  hitsabove5MIP_total_mip(ANtpDefVal::kFloat),
  hitsabove6MIP_total_mip(ANtpDefVal::kFloat),

  hitsabove1MIP_total_strips(ANtpDefVal::kInt),
  hitsabove2MIP_total_strips(ANtpDefVal::kInt),
  hitsabove3MIP_total_strips(ANtpDefVal::kInt),
  hitsabove4MIP_total_strips(ANtpDefVal::kInt),
  hitsabove5MIP_total_strips(ANtpDefVal::kInt),
  hitsabove6MIP_total_strips(ANtpDefVal::kInt),


  // position variables of those high hits
  high_hit_1_plane(ANtpDefVal::kInt),
  high_hit_2_plane(ANtpDefVal::kInt),
  high_hit_3_plane(ANtpDefVal::kInt),
  high_hit_4_plane(ANtpDefVal::kInt),
  high_hit_5_plane(ANtpDefVal::kInt),
  high_hit_6_plane(ANtpDefVal::kInt),
  high_hit_1_strip(ANtpDefVal::kInt),
  high_hit_2_strip(ANtpDefVal::kInt),
  high_hit_3_strip(ANtpDefVal::kInt),
  high_hit_4_strip(ANtpDefVal::kInt),
  high_hit_5_strip(ANtpDefVal::kInt),
  high_hit_6_strip(ANtpDefVal::kInt),
  high_hit_1_planeview(7),
  high_hit_2_planeview(7),
  high_hit_3_planeview(7),
  high_hit_4_planeview(7),
  high_hit_5_planeview(7),
  high_hit_6_planeview(7),
  high_hit_1_tpos(ANtpDefVal::kFloat),
  high_hit_2_tpos(ANtpDefVal::kFloat),
  high_hit_3_tpos(ANtpDefVal::kFloat),
  high_hit_4_tpos(ANtpDefVal::kFloat),
  high_hit_5_tpos(ANtpDefVal::kFloat),
  high_hit_6_tpos(ANtpDefVal::kFloat),
  high_hit_1_zpos(ANtpDefVal::kFloat),
  high_hit_2_zpos(ANtpDefVal::kFloat),
  high_hit_3_zpos(ANtpDefVal::kFloat),
  high_hit_4_zpos(ANtpDefVal::kFloat),
  high_hit_5_zpos(ANtpDefVal::kFloat),
  high_hit_6_zpos(ANtpDefVal::kFloat),
  showervtx_plane(ANtpDefVal::kInt),
  showervtx_u_pos(ANtpDefVal::kFloat),
  showervtx_v_pos(ANtpDefVal::kFloat),
  showerbeg_plane(ANtpDefVal::kInt)

{}

HighHitVars::~HighHitVars(){}

void HighHitVars::Reset()
{
  high_hit_1_ph = ANtpDefVal::kFloat;
  high_hit_2_ph = ANtpDefVal::kFloat;
  high_hit_3_ph = ANtpDefVal::kFloat;
  high_hit_4_ph = ANtpDefVal::kFloat;
  high_hit_5_ph = ANtpDefVal::kFloat;
  high_hit_6_ph = ANtpDefVal::kFloat;

  high_hit_1_mip = ANtpDefVal::kFloat;
  high_hit_2_mip = ANtpDefVal::kFloat;
  high_hit_3_mip = ANtpDefVal::kFloat;
  high_hit_4_mip = ANtpDefVal::kFloat;
  high_hit_5_mip = ANtpDefVal::kFloat;
  high_hit_6_mip = ANtpDefVal::kFloat;


  //Variables added at Jenny's request - find all the strips above a certain threshold and count them or total their ph:

  hitsabove1MIP_total_ph = ANtpDefVal::kFloat;
  hitsabove2MIP_total_ph = ANtpDefVal::kFloat;
  hitsabove3MIP_total_ph = ANtpDefVal::kFloat;
  hitsabove4MIP_total_ph = ANtpDefVal::kFloat;
  hitsabove5MIP_total_ph = ANtpDefVal::kFloat;
  hitsabove6MIP_total_ph = ANtpDefVal::kFloat;

  hitsabove1MIP_total_mip = ANtpDefVal::kFloat;
  hitsabove2MIP_total_mip = ANtpDefVal::kFloat;
  hitsabove3MIP_total_mip = ANtpDefVal::kFloat;
  hitsabove4MIP_total_mip = ANtpDefVal::kFloat;
  hitsabove5MIP_total_mip = ANtpDefVal::kFloat;
  hitsabove6MIP_total_mip = ANtpDefVal::kFloat;

  hitsabove1MIP_total_strips = ANtpDefVal::kInt;
  hitsabove2MIP_total_strips = ANtpDefVal::kInt;
  hitsabove3MIP_total_strips = ANtpDefVal::kInt;
  hitsabove4MIP_total_strips = ANtpDefVal::kInt;
  hitsabove5MIP_total_strips = ANtpDefVal::kInt;
  hitsabove6MIP_total_strips = ANtpDefVal::kInt;


  // position variables of those high hits
  high_hit_1_plane = ANtpDefVal::kInt;
  high_hit_2_plane = ANtpDefVal::kInt;
  high_hit_3_plane = ANtpDefVal::kInt;
  high_hit_4_plane = ANtpDefVal::kInt;
  high_hit_5_plane = ANtpDefVal::kInt;
  high_hit_6_plane = ANtpDefVal::kInt;
  high_hit_1_strip = ANtpDefVal::kInt;
  high_hit_2_strip = ANtpDefVal::kInt;
  high_hit_3_strip = ANtpDefVal::kInt;
  high_hit_4_strip = ANtpDefVal::kInt;
  high_hit_5_strip = ANtpDefVal::kInt;
  high_hit_6_strip = ANtpDefVal::kInt;
  high_hit_1_planeview = 7;
  high_hit_2_planeview = 7;
  high_hit_3_planeview = 7;
  high_hit_4_planeview = 7;
  high_hit_5_planeview = 7;
  high_hit_6_planeview = 7;
  high_hit_1_tpos = ANtpDefVal::kFloat;
  high_hit_2_tpos = ANtpDefVal::kFloat;
  high_hit_3_tpos = ANtpDefVal::kFloat;
  high_hit_4_tpos = ANtpDefVal::kFloat;
  high_hit_5_tpos = ANtpDefVal::kFloat;
  high_hit_6_tpos = ANtpDefVal::kFloat;
  high_hit_1_zpos = ANtpDefVal::kFloat;
  high_hit_2_zpos = ANtpDefVal::kFloat;
  high_hit_3_zpos = ANtpDefVal::kFloat;
  high_hit_4_zpos = ANtpDefVal::kFloat;
  high_hit_5_zpos = ANtpDefVal::kFloat;
  high_hit_6_zpos = ANtpDefVal::kFloat;
  showervtx_plane = ANtpDefVal::kInt;
  showervtx_u_pos = ANtpDefVal::kFloat;
  showervtx_v_pos = ANtpDefVal::kFloat;
  showerbeg_plane = ANtpDefVal::kInt;


}
