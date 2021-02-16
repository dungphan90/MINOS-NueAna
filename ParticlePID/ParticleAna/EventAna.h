#ifndef EVENTANA_H
#define EVENTANA_H

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleAna/Event.h"

class EventAna
{
	public:
	
		EventAna();
		~EventAna();


		void ana(ParticleObjectHolder * poh, Event * e);

};

#endif

