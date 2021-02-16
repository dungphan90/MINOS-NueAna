#ifndef ParticleConventions_H
#define ParticleConventions_H

namespace ParticleConventions
{
	double EnergyCorrection(double meu, int detector=-1, int release=-1, bool isMC=true);
}

namespace ClassType = ParticleConventions;

#endif

