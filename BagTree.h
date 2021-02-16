#ifndef BagTree_H
#define BagTree_H

#include "TObject.h"
#include <vector>

class TClonesArray;

class BagTree: public TObject
{
 public:
  BagTree();
  BagTree(const BagTree *bt);
  virtual ~BagTree();

  virtual void Print(Option_t* option="") const;
  void Reset();
  void Clear(Option_t* option = "");

  double hePID;
  double bt_var1;
  double bt_var2;
  double bt_var3;
 

  private:
   ClassDef(BagTree,1)
};

#endif //MCNNVARS_H
