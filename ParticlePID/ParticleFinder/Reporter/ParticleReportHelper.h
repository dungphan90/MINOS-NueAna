#ifndef PARTICLEREPORTHELPER_H
#define PARTICLEREPORTHELPER_H


#include <vector>
#include "NueAna/ParticlePID/ParticleFinder/Truth/ParticleTruthObject.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObjectHolder.h"
#include "NueAna/ParticlePID/ParticleFinder/ParticleObject.h"
#include "NueAna/ParticlePID/ParticleFinder/Reporter/ParticleReportObject.h"
#include <map>

class ParticleReportHelper
{
	public:

		ParticleReportHelper();	
		~ParticleReportHelper();

		void Reset(){found_map.clear(); truth_map.clear();};

		void addtruth(ParticleTruthObject *h);
		void addfound(ParticleObject *h);
		
		

		void Process(ParticleObjectHolder & p);

	private:

		std::map<int, std::map<int, std::map<int, double> > > found_map;
		std::map<int, std::map<int, std::map<int, double> > > truth_map;






};

#endif
