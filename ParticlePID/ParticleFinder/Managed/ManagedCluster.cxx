#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedCluster.h"
#include <math.h>

using namespace Managed;

int ManagedCluster::idcounter=1;

ClassImp(ManagedCluster)

ManagedCluster::ManagedCluster() :z(0),dz(0),t(0),dt(0),view(0),inuse(0),rms_t(0),e(0)
{
	id=idcounter++;
	Reset();
}


ManagedCluster::~ManagedCluster()
{
	hitplane.clear();
	hitstrip.clear();
	hite.clear();
	hitz.clear();
	hitt.clear();
	hit_id.clear();
	tsortmap.clear();
}


void ManagedCluster::Reset()
{
	z=0;
	dz=0;
	t=0;
	dt=0;
	view=0;
	inuse=0;
	rms_t=0;
	e=0;
	hitplane.clear();
	hitstrip.clear();
	hite.clear();
	hitz.clear();
	hitt.clear();
	hit_id.clear();
	tsortmap.clear();
	
	status=0;
	tmin=0;
	tmax=0;
	zmin=0;
	zmax=0;
}

void ManagedCluster::AdvanceID()
{
	id=idcounter++;
}

void ManagedCluster::Insert(double z,double t,double energy,int plane,int strip, int my_hit_id)
{
	hitplane.push_back(plane);
	hitstrip.push_back(strip);
	hite.push_back(energy);
	hitz.push_back(z);
	hitt.push_back(t);
	hit_id.push_back(my_hit_id);
	tsortmap.insert(std::make_pair(t,hitt.size()-1));
}

void ManagedCluster::ResetIDCounter()
{
	idcounter=1;
}
	
	
void ManagedCluster::Finalize()
{
	//printf("Finalizing cluster %d\n\thits (z,t,e)",id);

	double wz=0;
	double wt=0;
	e=0;
	for(unsigned int i=0;i<hitplane.size();i++)
	{
		//printf("(%f %f %f)",hitz[i],hitt[i],hite[i]);
		wz+=hitz[i]*hite[i];
		wt+=hitt[i]*hite[i];
		e+=hite[i];
	}
	if(e<0.001)return;
	wz/=e;
	wt/=e;
	
	
	t=wt;
	z=wz;
	
	double mdz=0;
	double mdt=0;
	for(unsigned int i=0;i<hitplane.size();i++)
	{
		if(fabs(z-hitz[i])>mdz)mdz=fabs(z-hitz[i]);
		if(fabs(t-hitt[i])>mdt)mdt=fabs(t-hitt[i]);
		
	}	
	
	dt=mdt;
	dz=mdz;
	
	
	double rmt=0;
	for(unsigned int i=0;i<hitplane.size();i++)
	{
		rmt+=(wt-hitt[i])*(wt-hitt[i])*hite[i]*hite[i];
	}
	rmt=sqrt(rmt);
	rms_t=rmt/e;
		
		
		
	if(hitt.size()>0)
	{
		tmin=hitt[0];
		tmax=hitt[0];
		zmin=hitz[0];
		zmax=hitz[0];
		for(unsigned int i=1;i<hitt.size();i++)
		{
			tmin=tmin<hitt[i]?tmin:hitt[i];
			tmax=tmax>hitt[i]?tmax:hitt[i];		
			zmin=zmin<hitz[i]?zmin:hitz[i];
			zmax=zmax>hitz[i]?zmax:hitz[i];			
		
		}
	
	
	}	
		
	//("\n\ttrms: %f\n",rms_t);	
	
}

