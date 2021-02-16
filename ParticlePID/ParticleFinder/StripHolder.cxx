#include "MessageService/MsgService.h"
#include "NueAna/ParticlePID/ParticleFinder/StripHolder.h"
#include <math.h>
#include <map>
#include <algorithm>

#include <sstream>

ClassImp(StripHolder)


//using namespace clusterer;

StripHolder::StripHolder():maxz(-10000),minz(10000),maxu(-10000),minu(10000),maxv(-10000),minv(10000)
{
		hits.clear();
}






StripHolder::~StripHolder()
{
}


void StripHolder::AddStrip(int myplane, int mystrip, double myenergy, double st, double sz, int view)
{
	StripHit h;
	h.plane=myplane;
	h.strip=mystrip;
	h.e=myenergy;
	h.t=st;
	h.z=sz;
	h.view=view;
	
	hits.push_back(h);
	
	minz=minz<sz?minz:sz;
	maxz=maxz>sz?maxz:sz;
	
	if(view==2)
	{
		minu=minu<st?minu:st;
		maxu=maxu>st?maxu:st;	
	}else if(view==3)
	{
		minv=minv<st?minv:st;
		maxv=maxv>st?maxv:st;		
	}
	
		
}



void StripHolder::Reset()
{
		hits.clear();
}
