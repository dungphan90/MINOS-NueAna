#ifndef HITMANAGER_H
#define HITMANAGER_H

#include <vector>
#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedHit.h"

#include "TObject.h"

namespace Managed
{

class HitManager : public TObject
{
	public:
		HitManager();
		~HitManager();

		void Reset();
		int InsertHit(int view, int plane, int strip, double z, double t, double e);
		ManagedHit *FindHit(int view, int plane, int strip);
		ManagedHit *FindHit(int view, double z, double t );
		ManagedHit *FindHit(int id);

		void ClearXTalk();

		std::vector<Managed::ManagedHit> GetAvailableHits();
		int GetHitCount(){return hits.size();};
		
	private:
		std::vector<Managed::ManagedHit> hits;

	ClassDef(HitManager,1);

};


}

#endif

