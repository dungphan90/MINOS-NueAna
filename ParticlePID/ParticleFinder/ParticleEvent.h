#ifndef PARTICLEEVENT_H
#define PARTICLEEVENT_H

#include "TObject.h"

#include "MCNtuple/NtpMCFluxInfo.h"

#include <vector>

class ParticleEvent : public TObject
{
	public:

		ParticleEvent();
		virtual ~ParticleEvent();



		double vtx_u;
		double vtx_v;
		double vtx_z;
		
		double visenergy;
		
		double sr_vtx_u;
		double sr_vtx_v;
		double sr_vtx_z;
		
		


		double minz;
		double minu;
		double minv;
		double maxz;
		double maxu;
		double maxv;
		
		double large_minz;
		double large_minu;
		double large_minv;
		double large_maxz;
		double large_maxu;
		double large_maxv;	
		
		double unused_e;
		double unused_e_avg;
		double unused_e_rms;	

		int nstrips;
		int nclusters;
		int unused_strips;

		int inFiducial;
		int contained;
		
	private:
	ClassDef(ParticleEvent,1)

};

#endif

