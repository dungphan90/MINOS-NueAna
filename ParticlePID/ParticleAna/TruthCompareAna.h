#ifndef TruthCompareAna_H
#define TruthCompareAna_H

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleAna/TruthCompare.h"

#include <vector>

class TruthCompareAna
{
	public:
	
		TruthCompareAna();
		~TruthCompareAna();


		void ana(ParticleObjectHolder * poh, TruthCompare * t);
		void Reset();
		
	private:	
		std::vector<int> MatchRecoParticle(Particle3D *p, ParticleObjectHolder *poh);
		void ProcessTrueParticles(ParticleObjectHolder *poh);
		

		std::vector<int>true_idx;//index into stdhep array for the particles that we conside (final state, sufficient energy,etc)
		std::vector<int>matched;//allready used?
		
		std::vector<int> MatchRecoElectron(Particle3D *p,ParticleObjectHolder *poh);
		std::vector<int> MatchRecoMuon(Particle3D *p,ParticleObjectHolder *poh);
		std::vector<int> MatchRecoProton(Particle3D *p,ParticleObjectHolder *poh);

		void ComputeStatistics(ParticleObjectHolder *poh, TruthCompare *t);	


		double reco_matched_visible_e;
		double reco_visible_e;
		
		double matchangle;
		
		double meupergev;
		double sigcorpermip;
};

#endif

