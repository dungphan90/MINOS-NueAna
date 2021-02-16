#ifndef PARTICLEOBJECT_H
#define PARTICLEOBJECT_H

#include "TObject.h"

#include <vector>

class ParticleObject : public TObject
{
	public:

		ParticleObject();
		//ParticleObject(size_t a){};
		virtual ~ParticleObject();



		Int_t a;

		int particle_id;
		int begstrip;
		int endstrip;
		int begplane;
		int endplane;

		std::vector<int> plane;
		std::vector<int> strip;
		std::vector<double> energy;
		std::vector<double> t;
		std::vector<double> z;
	
		double sumenergy;
		double fitenergy;

		double begt;
		double endt;
		double begz;
		double endz;
		double begu;
		double begv;


	private:
	ClassDef(ParticleObject,1)

};

#endif

