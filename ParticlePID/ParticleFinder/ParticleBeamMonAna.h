#ifndef PARTICLEBEAMMONANA_H
#define PARTICLEBEAMMONANA_H

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleBeamMon.h"
#include <string>

class ParticleBeamMonAna
{
	public:
	
		ParticleBeamMonAna();
		~ParticleBeamMonAna();


		void ana(ParticleObjectHolder * poh, ParticleBeamMon * e);
		BeamType::BeamType_t DetermineBeamType(VldContext vc);
		BeamType::BeamType_t DetermineBeamType(std::string file);

		std::string fname;

};

#endif

