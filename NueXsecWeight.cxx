#include <iostream>
#include "NueAna/NueXsecWeight.h"

ClassImp(NueXsecWeight)

NueXsecWeight::NueXsecWeight():
  xsecweight(1.)
{}

NueXsecWeight::NueXsecWeight(const NueXsecWeight *nuexs):
  xsecweight(nuexs->xsecweight)
{}

NueXsecWeight::~NueXsecWeight()
{}

void NueXsecWeight::Reset()
{
  xsecweight=1.0;
}

void NueXsecWeight::Print(Option_t * /*option*/) const
{

  std::cout<<"NueXsecWeight: totxsecweight "<<xsecweight<<std::endl;
}

