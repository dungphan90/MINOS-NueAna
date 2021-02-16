#ifndef PARTICLESANA_H
#define PARTICLESANA_H

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleAna/Particles.h"

#include "TH1F.h"

class ParticlesAna
{
	public:
	
		ParticlesAna();
		~ParticlesAna();


		void ana(ParticleObjectHolder * poh, Particles * e);

		static TH1F * rhist;
		static TH1F * viewTheta;
		static TH1F * viewPhi;


};

#endif

