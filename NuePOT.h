#ifndef NUEPOT_H
#define NUEPOT_H
#include "Record/RecRecordImp.h"
#include "Record/RecHeader.h"

class NuePOT : public TObject
{
 public:
  NuePOT();
  virtual ~NuePOT();

  double pot;
  double pot_nocut;
  int nruns;
  int nsnarls;
  int beamtype;

 private:
  ClassDef(NuePOT,2)

};

#endif //NUEPOT_H
