#ifndef MRCCAna_H
#define MRCCAna_H

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleAna/MRCC.h"

class MRCCAna
{
	public:
	
		MRCCAna();
		~MRCCAna();


		void ana(ParticleObjectHolder * poh, MRCC * e);

};

#endif

