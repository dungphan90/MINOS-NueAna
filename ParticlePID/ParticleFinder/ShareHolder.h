#ifndef SHAREHOLDER_H
#define SHAREHOLDER_H

#include <vector>
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"

struct hit
{

	hit():chain(-1), chainhit(-1), e(0), nshared(0){};

	int chain;
	int chainhit;
	double e;
	int nshared;
};


class ShareHolder{
	public:
		ShareHolder();
		~ShareHolder();



		void Insert(int chain, int chainhit, double e);
		int GetTotRemaining();
		int GetNumShared(int chain, int chainhit);
		double GetEShared(int chain, int chainhit);
		void Take(int chain, int chainhit, double e);
		
		void BulkInsert(std::vector<Particle3D *> &p3d);
	
		void Reset();
		void dump();



	private:
	
		std::vector<hit> hits;
		int find(int chain, int chainhit);


};

#endif

