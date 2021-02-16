#ifndef PRECORDANA_H
#define PRECORDANA_H


#include "NueAna/ParticlePID/ParticleAna/PRecord.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleBeamMon.h"

#include <map>


class PRecordAna
{
	public:

		PRecordAna();
		~PRecordAna();
		
		void ana(ParticleObjectHolder * h, PRecord * r,ParticleBeamMon *bmon);
	private:

};

#endif

