#include "ParticleEvent.h"
#include "TMath.h"

ClassImp(ParticleEvent)

ParticleEvent::ParticleEvent()
{
	vtx_u=0.0;
	vtx_v=0.0;
	vtx_z=0.0;

        visenergy=0.0;
        sr_vtx_u=0.0;
	sr_vtx_v=0.0;
	sr_vtx_z=0.0;

        nstrips=0;
	nclusters=0;
        minz=0;
	minu=0;
	minv=0;
	maxz=0;
	maxu=0;
	maxv=0;
               
	large_minz=100000;
	large_minu=100000;
	large_minv=100000;
	large_maxz=-100000;
	large_maxu=-100000;
	large_maxv=-100000;
                   
	unused_e=0;
	unused_strips=0;
	unused_e_avg=0;
	unused_e_rms=0;
	inFiducial=0;


	contained=0;
}


ParticleEvent::~ParticleEvent()
{

}

