///////////////////////////////////////////////////////////////////////////////////////
//                                                                                   //
// HighHitVars: variables that look at the individual high hits in the shower and    //
//              where they are located so as to see where they are in relation to    //
//              each other and the shower vertex                                     //
//                                                                                   //
// Author: Anna Holin                                                                //
// University College London                                                         //
// October 2007                                                                      //
//                                                                                   //
// Added variables in January 2013
///////////////////////////////////////////////////////////////////////////////////////

#ifndef HIGHHITVARS_H
#define HIGHHITVARS_H

#include "TObject.h"

class HighHitVars : public TObject
{

 public:
  HighHitVars();
  virtual ~HighHitVars();

  void Reset();

  //HighHitVars variables

  // Pulse height of highest hit in shower in sigcor (and 2nd highest, and 3rd highest, 4rth, 5th, 6th)
  Float_t high_hit_1_ph;
  Float_t high_hit_2_ph;
  Float_t high_hit_3_ph;
  Float_t high_hit_4_ph;
  Float_t high_hit_5_ph;
  Float_t high_hit_6_ph;

  // Pulse height of highest hit in shower in MIP (and 2nd highest, and 3rd highest, 4rth, 5th, 6th)
  Float_t high_hit_1_mip;
  Float_t high_hit_2_mip;
  Float_t high_hit_3_mip;
  Float_t high_hit_4_mip;
  Float_t high_hit_5_mip;
  Float_t high_hit_6_mip;

  //Variables added at Jenny's request - find all the strips above a certain threshold and count them or total their ph:                                                   
  Float_t hitsabove1MIP_total_ph;
  Float_t hitsabove2MIP_total_ph;
  Float_t hitsabove3MIP_total_ph;
  Float_t hitsabove4MIP_total_ph;
  Float_t hitsabove5MIP_total_ph;
  Float_t hitsabove6MIP_total_ph;

  Float_t hitsabove1MIP_total_mip;
  Float_t hitsabove2MIP_total_mip;
  Float_t hitsabove3MIP_total_mip;
  Float_t hitsabove4MIP_total_mip;
  Float_t hitsabove5MIP_total_mip;
  Float_t hitsabove6MIP_total_mip;

  Int_t hitsabove1MIP_total_strips;
  Int_t hitsabove2MIP_total_strips;
  Int_t hitsabove3MIP_total_strips;
  Int_t hitsabove4MIP_total_strips;
  Int_t hitsabove5MIP_total_strips;
  Int_t hitsabove6MIP_total_strips;

  // position variables of those high hits
  Int_t high_hit_1_plane;
  Int_t high_hit_2_plane;
  Int_t high_hit_3_plane;
  Int_t high_hit_4_plane;
  Int_t high_hit_5_plane;
  Int_t high_hit_6_plane;
  Int_t high_hit_1_strip;
  Int_t high_hit_2_strip;
  Int_t high_hit_3_strip;
  Int_t high_hit_4_strip;
  Int_t high_hit_5_strip;
  Int_t high_hit_6_strip;
  Char_t high_hit_1_planeview;
  Char_t high_hit_2_planeview;
  Char_t high_hit_3_planeview;
  Char_t high_hit_4_planeview;
  Char_t high_hit_5_planeview;
  Char_t high_hit_6_planeview;
  Float_t high_hit_1_tpos;
  Float_t high_hit_2_tpos;
  Float_t high_hit_3_tpos;
  Float_t high_hit_4_tpos;
  Float_t high_hit_5_tpos;
  Float_t high_hit_6_tpos;
  Float_t high_hit_1_zpos;
  Float_t high_hit_2_zpos;
  Float_t high_hit_3_zpos;
  Float_t high_hit_4_zpos;
  Float_t high_hit_5_zpos;
  Float_t high_hit_6_zpos;
  Int_t showervtx_plane;
  Float_t showervtx_u_pos;
  Float_t showervtx_v_pos;
  Int_t showerbeg_plane;



 private:

  ClassDef(HighHitVars,6)
};

#endif// HIGHHITVARS_H
