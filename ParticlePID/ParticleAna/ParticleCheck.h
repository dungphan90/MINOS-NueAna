#ifndef ParticleCheck_H
#define ParticleCheck_H

#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObject.h"
#include "NueAna/ParticlePID/ParticleAna/PRecord.h"
#include "TChain.h"
#include "TH2F.h"

#define nCut 2
#define nType 5

class ParticleCheck
{

	public:
		ParticleCheck(const char*fname);
		~ParticleCheck();
		void Run();
int PassesPreselec();

	TChain * input;
	ParticleObjectHolder *poh;

	TChain * inputPA;
	PRecord *prec;	
		
	void MakeHistos();
	void SaveHistos();
	void Process(int cut, int type, double weight);

	TH2F *elec_single_match_eres[nCut][nType];
	TH2F *elec_single_match_angle[nCut][nType];

	TH2F *elec_single_pi0_eres[nCut][nType];
	TH2F *elec_single_pi0_angle[nCut][nType];

	TH2F *elec_matched_elec_eres[nCut][nType];
	TH2F *elec_matched_pi0_eres[nCut][nType];

	TH2F *elec_matched_type_vs_recoE[nCut][nType];
	TH2F *elec_matched_typefrac_vs_recoE[nCut][nType];

};


#endif

