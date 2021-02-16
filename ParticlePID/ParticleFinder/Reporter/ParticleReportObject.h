#ifndef PARTICLEREPORTOBJECT_H
#define PARTICLEREPORTOBJECT_H

#include "TObject.h"

#include <vector>

class ParticleReportObject : public TObject
{
	public:

		ParticleReportObject();
		virtual ~ParticleReportObject();



//		std::vector<int> found_index;
//		std::vector<int> best_strip_count_match_truth_index;
//		std::vector<int> best_energy_match_truth_index;

//		double frac_e_correct;
//		double frac_strips_correct;
		
//		double frac_e_wrong;
//		double frac_strips_wrong;

        	int totstrips;
        	int matched_strips;
        	double tote;
        	double matched_e;







	private:
	ClassDef(ParticleReportObject,1)

};

#endif
