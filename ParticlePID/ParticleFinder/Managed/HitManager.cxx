#include "NueAna/ParticlePID/ParticleFinder/Managed/HitManager.h"

#include <math.h>
#include <map>

#include "MessageService/MsgService.h"

#include <sstream>

#include <cstdlib>

CVSID("$Id: HitManager.cxx,v 1.3 2009/09/11 05:00:40 gmieg Exp $");



using namespace Managed;

ClassImp(HitManager)

HitManager::HitManager()
{
	Reset();
}

HitManager::~HitManager(){}


Managed::ManagedHit  * HitManager::FindHit(int id)
{
	for(unsigned int i=0;i<hits.size();i++)
		if(hits[i].id==id)return &hits[i];

	return 0;
}

void HitManager::Reset()
{
	hits.clear();
	ManagedHit::ResetIDCounter();
}

int HitManager::InsertHit(int view, int plane, int strip, double z, double t, double e)
{
	Managed::ManagedHit  h(view,plane,strip,z,t,e);
	hits.push_back(h);
	return h.id;
}

Managed::ManagedHit  *HitManager::FindHit(int /*view*/, int plane, int strip)
{
	for(unsigned int i=0;i<hits.size();i++)
		if(hits[i].GetPlane()==plane && hits[i].GetStrip()==strip)return &hits[i];
	return 0;
}


Managed::ManagedHit  *HitManager::FindHit(int /*view*/, double z, double t)
{
	for(unsigned int i=0;i<hits.size();i++)
		if(fabs(hits[i].GetZ()-z)<0.01 && fabs(hits[i].GetT()-t)<0.01)return &hits[i];
	return 0;
}


std::vector<ManagedHit> HitManager::GetAvailableHits()
{
	std::vector<ManagedHit> ret;
	for(unsigned int i=0;i<hits.size();i++)
		if(hits[i].GetERemaining()>0)ret.push_back(hits[i]);
		
	return ret;

}


void HitManager::ClearXTalk()
{

	
		std::vector<int> view;
		std::vector<int> plane;
		std::vector<int> strip;
		std::vector<double> t;
		std::vector<double> z;
		std::vector<double> energy;
		
		
		std::multimap<double,int> sorter_map;

		std::map<int, std::map<int, int> > loc_map;  //plane, strip, index in vectors


	for(unsigned int i=0;i<hits.size();i++)
	{
		view.push_back(hits[i].GetView());
		plane.push_back(hits[i].GetPlane());
		strip.push_back(hits[i].GetStrip());
		t.push_back(hits[i].GetT());
		z.push_back(hits[i].GetZ());
		energy.push_back(hits[i].GetERemaining());
		
		
	
	
	}


	double thresh = 3; //number of mips below which we don't care about looking for xtalk....

	for(unsigned int i=0;i<plane.size();i++)
	{
		sorter_map.insert(std::pair <double,int>(energy[i],i));
		loc_map[plane[i]][strip[i]]=i;
	}
	
	ostringstream os;
	
	double reassignedxtalke=0.0;
	
	int foundxtalk=0;
	std::map<int, std::map<int, int> >::iterator p_iter;
	for(p_iter=loc_map.begin();p_iter!=loc_map.end(); p_iter++)
    {
    	std::map<int, int>::iterator s_iter;
        for(s_iter=p_iter->second.begin();s_iter!=p_iter->second.end(); s_iter++)
        {
			if (energy[s_iter->second]<thresh)continue;
			
			std::map<int, int>::iterator s_iter2;
			for(s_iter2=p_iter->second.begin();s_iter2!=p_iter->second.end(); s_iter2++)
        	{
				if(s_iter==s_iter2)continue;
				int sp=abs(s_iter->first-s_iter2->first);
				os<<"sep "<<sp<<" view "<<view[s_iter->second]<<" plane "<<plane[s_iter->second]<<" high "<<energy[s_iter->second]<<" low "<<energy[s_iter2->second]<<"\n";
				
				if( sp>9 && sp<14)//sp == 13)
				{
					if(energy[s_iter2->second] < 0.3 * energy[s_iter->second]) //suspected crosstalk hit is < 30% of the high hit
					{
						reassignedxtalke+=energy[s_iter2->second];
						
						energy[s_iter->second]+=energy[s_iter2->second];
						energy[s_iter2->second]=0.0;
						foundxtalk=1;
					}
				}
			
			}
		}
	}


//	if(reassignedxtalke>0)printf("XTALK FILTER --- %f reassigned!\n",reassignedxtalke);


	if(foundxtalk)
	{
		std::vector<int>tplane;
		std::vector<int>tstrip;
		std::vector<int>tview;
		std::vector<double>tt;
		std::vector<double>tz;
		std::vector<double>tenergy;
	
	
		for(unsigned int i=0;i<energy.size();i++)
		{
			if(energy[i]>0.001)
			{
				tplane.push_back(plane[i]); 
				tstrip.push_back(strip[i]); 
				tenergy.push_back(energy[i]); 
				tt.push_back(t[i]); 
				tz.push_back(z[i]);
				tview.push_back(view[i]);
			}
		}
		
		plane=tplane;
		strip=tstrip;
		view=tview;
		t=tt;
		z=tz;
		energy=tenergy;
		
		os<<"Removing XTalk ---- "<<hits.size() <<" hits merged to make ";
		Reset();
		for(unsigned int i=0;i<plane.size();i++)
		{
			Managed::ManagedHit  h(view[i],plane[i],strip[i],z[i],t[i],energy[i]);
			hits.push_back(h);
		}
		os << (int)hits.size() <<" hits\n";
		
	}


	sorter_map.clear();
	loc_map.clear();
	
	MSG("HitManager",Msg::kDebug)<<os;

}



