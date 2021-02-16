#include "NueAna/ParticlePID/ParticleAna/Tools/ParticleConventions.h"
#include "Conventions/ReleaseType.h"

#include <cstdlib>





double ParticleConventions::EnergyCorrection(double meu, int detector, int type, bool /*isMC*/)
{
       if(ReleaseType::IsDogwood(type)&& (ReleaseType::IsDaikon(type)||ReleaseType::IsData(type)))
        {
                  if(detector==2) //far
                  {
                    float offset =   0.489987;
                    float slope  =  0.0387301;
                    return meu*slope+offset;
                }else if(detector==1) //near
                {
                        float offset = 0.4803261;
                        float slope = 0.03799819;
                        return meu*slope + offset;
                }
        }

	printf("UNKNOWN TYPE in ParticleConventions::EnergyCorrection!\n");
	exit(1);
}

