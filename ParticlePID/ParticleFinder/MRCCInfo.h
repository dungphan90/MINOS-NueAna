#ifndef MRCCInfo_H
#define MRCCInfo_H


#include "TObject.h"


#include <vector>
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"



class MRCCInfo : public TObject

{

	public:
	
		MRCCInfo();
		virtual ~MRCCInfo();

		void Reset();
		
		Particle3D removedmuon;
		int stage;
     private:
        ClassDef(MRCCInfo,2)


};

#endif

