#include "NueAna/ParticlePID/ParticleAna/EventAna.h"


EventAna::EventAna()
{}

EventAna::~EventAna()
{}


void EventAna::ana(ParticleObjectHolder * poh, Event * e)
{
	e->vtx_u = poh->event.vtx_u;	
	e->vtx_v = poh->event.vtx_v;
	e->vtx_z = poh->event.vtx_z;
	
	e->sr_vtx_u = poh->event.sr_vtx_u;
	e->sr_vtx_v = poh->event.sr_vtx_v;
	e->sr_vtx_z = poh->event.sr_vtx_z;
	
	
	e->min_u=poh->event.minu;
	e->min_v=poh->event.minv;
	e->min_z=poh->event.minz;
	e->max_u=poh->event.maxu;
	e->max_v=poh->event.maxv;
	e->max_z=poh->event.maxz;
	
	e->visenergy = poh->event.visenergy;
	
	e->nstrips=poh->event.nstrips;
	e->nclusters=poh->event.nclusters;
	
	e->large_minz=poh->event.large_minz;
	e->large_minu=poh->event.large_minu;
	e->large_minv=poh->event.large_minv;
	e->large_maxz=poh->event.large_maxz;
	e->large_maxu=poh->event.large_maxu;
	e->large_maxv=poh->event.large_maxv;
		
	e->unused_e=poh->event.unused_e;
	e->unused_strips=poh->event.unused_strips;
	e->unused_e_avg=poh->event.unused_e_avg;
	e->unused_e_rms=poh->event.unused_e_rms;
	e->inFiducial = poh->event.inFiducial;
	e->contained = poh->event.contained;
	e->foundlongmuon=poh->eventquality.foundlongmuon;
	e->foundprimaryshower=poh->eventquality.foundprimaryshower;
	
}


