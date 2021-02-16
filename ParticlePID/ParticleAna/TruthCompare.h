#ifndef TruthCompare_H
#define TruthCompare_H

#include <vector>
#include "TObject.h"


class TruthCompare :public TObject
{

	public:
	
		TruthCompare();
		virtual ~TruthCompare();
		void Reset();

		virtual void Clear(Option_t* option = "");


		
		double frac_e_found;
		double frac_particles_found;
		
		double true_visible_e;
		double reco_matched_visible_e;
		double reco_visible_e;
		double matchangle;

		//find the lepton and its energy from CC events
		double truelepE;

		double emenergy;		

		

		int particles_matched_to_true; //number of reco particles matched to true particles
		int possible_true_particles; //number of true particles of large enough energies to see
		int truelepType;	
		int stage;


		//there is an entry in this vector for each recoed particle
		//each entry is a vector of stdhep indicies to which this reco particle is matched
		std::vector< std::vector<int> > reco_to_true_match;
		std::vector<int> reco_matched;
		
		//list of true particles that are large enough to be seen that are not matched
		std::vector<int> true_not_recoed;
		
		//list of recoed particles that we can't find a particle in true to match to
		std::vector<int> reco_not_matched;
  
	private:
		void Init();
		ClassDef(TruthCompare,4)

};

#endif

