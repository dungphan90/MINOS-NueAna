///////////////////////////////////////////////////////////////////////////
// 
// FracVar: variables for TJ's analysis
//
// References: Adam Para: NuMI-NOTE-SIM-0284
//             Stan Wojcicki: NuMI-NOTE-SIM-0290
//             http://www.slac.stanford.edu/~tjyang/talks/nue20050318/nue20050318.pdf
// tjyang@stanford.edu
//
///////////////////////////////////////////////////////////////////////////
#ifndef FRACVAR_H
#define FRACVAR_H

#include "TObject.h"

class FracVar : public TObject
{

 public:
  FracVar();
  virtual ~FracVar();

  //virtual void Draw(Option_t *option);
  //virtual void Print(Option_t *option) const;
  void Reset();

  //FracVar variables
  //maximum fractions of the observed total_energy deposited in 1,2...6 plane(s).
  Float_t fract_1_plane;
  Float_t fract_2_planes;
  Float_t fract_3_planes;
  Float_t fract_4_planes;
  Float_t fract_5_planes;
  Float_t fract_6_planes;

  //fractions of the observed total_energy deposited in 1,2...6 channel(s) (i.e. strips) with the highest pulse height
  Float_t fract_2_counters;
  Float_t fract_4_counters;
  Float_t fract_6_counters;
  Float_t fract_8_counters;
  Float_t fract_10_counters;
  Float_t fract_12_counters;
  Float_t fract_14_counters;
  Float_t fract_16_counters;
  Float_t fract_18_counters;
  Float_t fract_20_counters;

  Float_t fract_road; //fraction of energy deposited in a narrow road
  Float_t fract_asym; //asymmetry of u and v plane charges

  Int_t shw_max;//position (with respect to the first plane) of the plane with the highest observed energy in the shower core
  Int_t dis2stp;   //distance between the 2 strips that have the highest PH.

  Int_t shw_nstp;  //no. of strips
  Int_t shw_npl;   //no. of planes

  Float_t shw_slp;   //slope of shower dir
  Float_t vtxph;     //fraction of ph deposited around vertex

  Int_t passcuts;  //pass the crude cuts?
  Float_t pid;       //ann output, e-like close to 1
  Float_t pid1;      //updated ann output  3/15/06

 private:

  ClassDef(FracVar,6)
};

#endif// FRACVAR_H
