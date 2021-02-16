#ifndef DetailedParticle_H
#define DetailedParticle_H

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObject.h"
#include "NueAna/ParticlePID/ParticleAna/PRecord.h"
#include "TChain.h"


class DetailedParticle
{

	public:
		DetailedParticle(const char*fname);
		~DetailedParticle();
		void Run();
int PassesPreselec();

	TChain * input;
	ParticleObjectHolder *poh;

	TChain * inputPA;
	PRecord *prec;	

};


#endif

