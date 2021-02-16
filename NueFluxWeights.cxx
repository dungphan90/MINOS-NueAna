#include <iostream>
#include "NueAna/NueFluxWeights.h"

ClassImp(NueFluxWeights)

NueFluxWeights::NueFluxWeights():
  totbeamweight(1.),
  detectorWeight(1.0),
  kflukweight(1.),
  totskzpweight(1.),
  skzpTrkEnergy(1.0),
  skzpShwEnergy(1.0),
  RPtotbeamweight()
{
  skzpConfig = "";
}

NueFluxWeights::NueFluxWeights(const NueFluxWeights *nuefw):
  totbeamweight(nuefw->totbeamweight),
  detectorWeight(nuefw->detectorWeight),
  kflukweight(nuefw->kflukweight),
  totskzpweight(nuefw->totskzpweight),
  skzpTrkEnergy(nuefw->skzpTrkEnergy),
  skzpShwEnergy(nuefw->skzpShwEnergy)
{
  skzpConfig = nuefw->skzpConfig;
  for(unsigned int i=0;i<nuefw->RPtotbeamweight.size();i++){
    RPtotbeamweight.push_back(nuefw->RPtotbeamweight[i]);
  }
}

NueFluxWeights::~NueFluxWeights()
{
  //  std::cout<<"in NueFluxWeights destructor"<<std::endl;
   RPtotbeamweight.clear();
}

void NueFluxWeights::Clear(Option_t* /* option */)
{
}

void NueFluxWeights::Reset()
{
  kflukweight=1.;
  totbeamweight=1.;
  totskzpweight=1.;
  detectorWeight = 1.0;
  RPtotbeamweight.clear();
}

void NueFluxWeights::Print(Option_t * /*option*/) const
{
  std::cout<<"Det weight "<<detectorWeight<<std::endl;
  std::cout<<"Beam weight "<<totbeamweight<<std::endl;
  std::cout<<"Kfluk weight "<<kflukweight<<std::endl;
  std::cout<<"Tot Beam weight "<<totskzpweight<<std::endl;

  for(unsigned int i=0;i<RPtotbeamweight.size();i++){
    std::cout<<"Run Period "<<i<<" weight "<<RPtotbeamweight[i]<<std::endl;
  }


}



