///////////////////////////////////////////////////////////////////////////
// 
// ShieldRejVar: variables from Pedro's Shield package
//
//
///////////////////////////////////////////////////////////////////////////
#ifndef SHIELDREJVAR_H
#define SHIELDREJVAR_H

#include "TObject.h"

class ShieldRejVar : public TObject
{

 public:
  ShieldRejVar();
  virtual ~ShieldRejVar();

  void Reset();

  //ShieldRejVar variables
  Int_t ShieldHit;

 private:
  ClassDef(ShieldRejVar,3)
};

#endif// SHIELDREJVAR_H
