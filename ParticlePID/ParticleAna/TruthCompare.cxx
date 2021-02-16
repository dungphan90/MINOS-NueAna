#include "NueAna/ParticlePID/ParticleAna/TruthCompare.h"


ClassImp(TruthCompare)

TruthCompare::TruthCompare() 
{
	Init();
}

TruthCompare::~TruthCompare()
{
	Init(); //to empty vectors...
}


void TruthCompare::Reset()
{

	Init();
}



void TruthCompare::Clear(Option_t* /* option */) {
  // Purpose: Clear memory allocated to arrays so that record can
  // be reused.  

	Init();

}
 
 
 
void TruthCompare::Init() {
  // 
  // Purpose: Initialize ntuple TClonesArrays
  //
	reco_to_true_match.clear();
	true_not_recoed.clear();
	reco_not_matched.clear();
	reco_matched.clear();

	frac_e_found=0;
	frac_particles_found=0;
		
	true_visible_e=0;
	reco_matched_visible_e=0;
	reco_visible_e=0;
		
	particles_matched_to_true=0; 
	possible_true_particles=0;
	matchangle=0;

	stage=0;
	
	truelepE=0;
	truelepType=0;
	emenergy=0;
  
}
