#ifndef PARTICLETRUTHHELPER_H
#define PARTICLETRUTHHELPER_H


#include <vector>
#include "NueAna/ParticlePID/ParticleFinder/Truth/ParticleTruthObject.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include <map>

class ParticleTruthHelper
{
	public:

		ParticleTruthHelper();	
		~ParticleTruthHelper();

		void Reset(){main_map.clear(); test_map.clear();};
		
		void AddStrip(int mypid, int trkid, int myplane, int mystrip, double myenergy,int view);

		void Process(ParticleObjectHolder & p);

	private:


	std::map<int, std::map<int, std::map<int, std::map<int, double> > > > main_map;

	std::map<int, std::map<int, std::map<int, std::map<int, int> > > > test_map;

	std::map<int, int> pid_map;

};

#endif
