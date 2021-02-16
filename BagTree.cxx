#include <iostream>
#include "TClonesArray.h"
#include "NueAna/BagTree.h"
#include "NueAna/MCNNBestMatch.h"
#include "AnalysisNtuples/ANtpDefaultValue.h"

ClassImp(BagTree)

BagTree::BagTree() 
{
  Reset(); 
}

BagTree::BagTree(const BagTree *mv):
  hePID(mv->hePID),
  bt_var1(mv->bt_var1),
  bt_var2(mv->bt_var2),
  bt_var3(mv->bt_var3)
{
}

BagTree::~BagTree()
{
  //  std::cout<<"in BagTree destructor"<<std::endl;
}

void BagTree::Clear(Option_t* /* option */)
{
}

void BagTree::Reset()
{
  hePID= ANtpDefVal::kDouble; //<-- open variable.
  bt_var1= ANtpDefVal::kDouble;
  bt_var2= ANtpDefVal::kDouble;
  bt_var3= ANtpDefVal::kDouble;

}

void BagTree::Print(Option_t * /*option*/) const
{
}



