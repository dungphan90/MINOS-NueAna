#include "NueAna/ParticlePID/ParticleFinder/Truth/ParticleTruthHelper.h"

#include "Plex/PlexPlaneId.h"

#include "UgliGeometry/UgliGeomHandle.h"
#include "UgliGeometry/UgliScintPlnHandle.h"
#include "Plex/PlexStripEndId.h"
#include "UgliGeometry/UgliStripHandle.h"

#include "TVector3.h"

ParticleTruthHelper::ParticleTruthHelper()
{

//	sorter_map.clear();
}

ParticleTruthHelper::~ParticleTruthHelper()
{
}


void ParticleTruthHelper::AddStrip(int pid, int trkid, int plane, int strip, double myenergy, int view)
{



	if(test_map[view][plane][strip][trkid]!=1)	
	main_map[view][plane][strip][trkid]=0;   //reset if not already in use... there must be a better way to do this, if we can can guarente the default value for a double, otherwise, we may need to use an object instead of a double....


	main_map[view][plane][strip][trkid]+=myenergy;
	test_map[view][plane][strip][trkid]=1;
	pid_map[trkid]=pid;


}


void ParticleTruthHelper::Process(ParticleObjectHolder & p)
{


	VldContext vldc(Detector::kFar , SimFlag::kData , VldTimeStamp());
	UgliGeomHandle geo(vldc);            


	//iterate over map, and for each entry, make a new entry to add in p

	std::map<int, std::map<int, std::map<int, std::map<int, double> > > >::iterator it;
	std::map<int, std::map<int, std::map<int, double> > >::iterator it1;
	std::map<int, std::map<int, double> > ::iterator it2;
	std::map<int, double>::iterator it3;




	for(it=main_map.begin();it!=main_map.end();it++)
	for(it1=it->second.begin();it1!=it->second.end();it1++)
	{
		double thisz=0;

		PlexPlaneId pl(Detector::kFar, it1->first);
		UgliScintPlnHandle  h = geo.GetScintPlnHandle(pl); 
		thisz = h.GetCenter()[2];

		for(it2=it1->second.begin();it2!=it1->second.end();it2++)
		{
			double thist=0;
			PlexStripEndId s(pl,it2->first);
			UgliStripHandle a = h.GetStripHandle(s);   
			thist = a.GetTPos();

	
			for(it3=it2->second.begin();it3!=it2->second.end();it3++)
			{	
/*
			        ParticleTruthObject *d = (ParticleTruthObject*)p.particletruth->New(p.particletruth->GetEntries());
				d->plane=it1->first;
				d->strip=it2->first;
				d->trkid=it3->first;
				d->pid=pid_map[it3->first];
				d->sumenergy=it3->second;
				d->view=it->first;
				d->t=thist;
				d->z=thisz;
*/			}
		}
	}
}


