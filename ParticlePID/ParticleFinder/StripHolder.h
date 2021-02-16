#ifndef STRIPHOLDER_H
#define STRIPHOLDER_H

#include <vector>
//#include "NueAna/ParticlePID/ParticleFinder/Cluster.h"
#include <map>



#include "TObject.h"
#include "StripHit.h"

class StripHolder : public TObject
{
	public:
	
		StripHolder();
		~StripHolder();
		
		void Reset();
		void AddStrip(int myplane, int mystrip, double myenergy, double st, double sz, int view);


		double maxz;
		double minz;
		double maxu;
		double minu;
		double maxv;
		double minv;
		std::vector<StripHit> hits;
		
		

     private:
        ClassDef(StripHolder,1)		
		
};



#endif

