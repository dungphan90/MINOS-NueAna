#ifndef NUEXSECWEIGHT_H
#define NUEXSECWEIGHT_H
#include "TObject.h"

class NueXsecWeight: public TObject
{
 public:
  NueXsecWeight();
  NueXsecWeight(const NueXsecWeight *nuexs);
  virtual ~NueXsecWeight();

  virtual void Print(Option_t* option="") const;
  void Reset();

  double xsecweight;

 private:

  ClassDef(NueXsecWeight,1)

};

#endif //NUEXSECWEIGHT_H
