#include "NueAna/ParticlePID/ParticleFinder/POT.h"

ClassImp(POT)

POT::POT()
{
	Reset();
}

POT::~POT()
{}

  
void POT::Reset()
{
	pot=0.;
  	nruns=0;
  	nsnarls=0;
  	beamtype=0;
}


