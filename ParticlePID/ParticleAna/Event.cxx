#include "NueAna/ParticlePID/ParticleAna/Event.h"

ClassImp(Event)

Event::Event()
{
	Init();

}

Event::~Event()
{}


void Event::Clear(Option_t* /* option */) {
  // Purpose: Clear memory allocated to arrays so that record can
  // be reused.  

	Init();

}
 
 
 
void Event::Init() {
  // 
  // Purpose: Initialize ntuple TClonesArrays
  //
	foundlongmuon=0;
	foundprimaryshower=0;
	pidA=-1;
	pidB=-1;
	pidC=-1;
	pidD=-1;
	pidE=-1;
	pidF=-1;


	contained=0;
	vtx_u=0;
	vtx_v=0;
	vtx_z=0;
	sr_vtx_u=0;
	sr_vtx_v=0;
	sr_vtx_z=0;
	min_u=0;
	min_v=0;
	min_z=0; 
	max_u=0;
	max_v=0;
	max_z=0;
	visenergy=0;
	nstrips=0;
	nclusters=0; 
	large_minz=0;
	large_minu=0;
	large_minv=0;
	large_maxz=0;
	large_maxu=0;
	large_maxv=0;	
	unused_e=0;
	unused_strips=0;
	unused_e_avg=0;
	unused_e_rms=0;
	inFiducial=0;
}



