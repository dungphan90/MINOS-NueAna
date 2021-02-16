#ifndef POT_H
#define POT_H
#include "Record/RecRecordImp.h"
#include "Record/RecHeader.h"

class POT : public TObject
{
 public:
  POT();
  virtual ~POT();

  double pot;
  int nruns;
  int nsnarls;
  int beamtype;
  
  void Reset();

 private:
  ClassDef(POT,2)

};

#endif //POT_H
