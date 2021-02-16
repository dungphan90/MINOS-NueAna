#include "NueAna/ParticlePID/ParticleFinder/MCTrue.h"

// MINOS headers



ClassImp(MCTrue)

MCTrue::MCTrue() 
{
  this -> Init();  

}

MCTrue::~MCTrue()
{
 stdhep.clear();
}



void MCTrue::Clear(Option_t* /* option */) {
  // Purpose: Clear memory allocated to arrays so that record can
  // be reused.  

	Init();

}
 
 
 
void MCTrue::Init() {
  // 
  // Purpose: Initialize ntuple TClonesArrays
  //

	vtx_u=0;
	vtx_v=0;
	vtx_z=0;
	inu=0;
	iresonance=-1;
	iaction=-1;
	inunoosc=0;
	nuenergy=0;
	oscprob=1;
	totbeamweight=1;
	osc_L=0;
	osc_dm2=0;
	osc_sinth23=0;
	osc_sin2th13=0;
	type=-1;
	visible_energy=0;

  


  stdhep.clear();
  
}

