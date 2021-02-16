#ifndef PARTICLETYPE_H
#define PARTICLETYPE_H

#include "TObject.h"

class ParticleType : public TObject
{

	public:
		
			ParticleType();
			~ParticleType();
	
	
typedef enum Eparticletype
{
	em,
	emshort,
	muon,
	prot,
	pi0,
	uniquemuon


}particletype_t;


	
			particletype_t type;
			int start;
			int stop;
			
			
			
						ClassDef(ParticleType,1)
};

#endif

