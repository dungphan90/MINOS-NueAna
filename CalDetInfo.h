///////////////////////////////////////////////////////////////////////////
// 
// CalDetInfo: variables from Mike and Trish's CalDetDST package
//
///////////////////////////////////////////////////////////////////////////
#ifndef CALDETINFO_H
#define CALDETINFO_H

#include "TObject.h"

class CalDetInfo : public TObject
{

 public:
  CalDetInfo();
  virtual ~CalDetInfo();

  void Zero();
  void Reset();

  //CalDetInfo variables
  Float_t beamp;
  Int_t inct;
  Int_t pid;
  Float_t olchi2;
  Int_t p0stripmaxmip;

 private:

  ClassDef(CalDetInfo,1)
};

#endif// CALDETINFO_H
