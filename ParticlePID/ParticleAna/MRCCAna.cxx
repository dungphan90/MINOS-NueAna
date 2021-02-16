#include "NueAna/ParticlePID/ParticleAna/MRCCAna.h"
#include "NueAna/ParticlePID/ParticleFinder/MRCCInfo.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"


MRCCAna::MRCCAna()
{}

MRCCAna::~MRCCAna()
{}


void MRCCAna::ana(ParticleObjectHolder * poh, MRCC * e)
{
	MRCCInfo * mi = poh->mrcc;
	if(!mi)
	{
		e->hasMRCC=0;
		return;
	}
	
	Particle3D *p=&mi->removedmuon;
	if(!p)
	{
		e->hasMRCC=0;
		return;	
	}
	
	e->hasMRCC=1;
	
	e->stage = mi->stage;
	e->start_u=p->start_u;
	e->end_u=p->end_u;
	e->start_v=p->start_v;
	e->end_v=p->end_v;
	e->start_z=p->start_z;
	e->end_z=p->end_z;
		
	e->sum_e=p->sum_e;

	e->muonfrac=p->muonfrac;
	e->entries=p->entries;
	
	e->theta=p->theta; 	
	e->phi=p->phi; 	
				
	e->calibrated_energy=p->calibrated_energy;
	e->avg_rms_t=p->avg_rms_t;



		double particle_s=0;		
		for(int j=0;j<p->entries;j++)
		{
			//double u=p->u[j]-poh->event.vtx_u;
			//double v=p->v[j]-poh->event.vtx_v;
	
			
			if(j>0)
			{
				double du = p->u[j]-p->u[j-1];
				double dv = p->v[j]-p->v[j-1];
				double dz = p->z[j]-p->z[j-1];
				
				double dd=sqrt(du*du+dv*dv+dz*dz);
				if(dd<100)particle_s+=dd; //in case we have a single bogus value, we will skip it
			}
		
		}
		
	e->particle_s=particle_s;	
	
}


