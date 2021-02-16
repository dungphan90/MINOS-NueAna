#ifndef MCTRUEANA_H
#define MCTRUEANA_H

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleAna/MCTrue.h"



#include "MCReweight/SKZPWeightCalculator.h"
#include "MCNtuple/NtpMCTruth.h"
#include "MCReweight/Zbeam.h"
#include "MCReweight/Zfluk.h"
#include "MCReweight/Kfluk.h"

#include "NueAna/ParticlePID/ParticleFinder/ParticleBeamMon.h"


#include "TF1.h"

class MCTrueAna
{
	public:
	
		MCTrueAna();
		~MCTrueAna();


		void ana(ParticleObjectHolder * poh, MCTrue * e,ParticleBeamMon * bmon);

		double OscillationProb(TF1* f2, int ntype, double NuE, double sinth23, double sin2th13);

		static TF1 * osceq;
		static SKZPWeightCalculator *skzpCalc;
		
		
		
	private:	
		double osc_dm2;
 		double osc_L; 
		double osc_sinth23;
		double osc_sin2th13;

};

#endif

