#ifndef PARTICLETRUTHOBJECT_H
#define PARTICLETRUTHOBJECT_H

#include "TObject.h"

#include <vector>

class ParticleTruthObject : public TObject
{
	public:

		ParticleTruthObject();
		virtual ~ParticleTruthObject();

		int pid;
		int strip;
		int plane;
		
		double t;
		double z;

		double sumenergy;

		int trkid;

		int view;

	private:
	ClassDef(ParticleTruthObject,1)

};

#endif
