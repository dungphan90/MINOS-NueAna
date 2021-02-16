#ifndef NUERWHELPERS_H
#define NUERWHELPERS_H

class NtpMCTruth;
class ANtpTruthInfoBeam;
class ANtpTruthInfoBeamNue;

namespace NueRWHelpers
{
  float Oscillate(NtpMCTruth *mcth,float L, float dm2, float theta23, float UE32);
  float Oscillate(ANtpTruthInfoBeam *ib,float L, float dm2, float theta23, float UE32);
  float Oscillate(int nuFlavor, int nonOscNuFlavor, float Energy, 
                  float L, float dm2, float theta23, float U);
  float Oscillate(ANtpTruthInfoBeamNue *ib);

  
  //hierarchy = 1 => NO; -1 => RO  
  float OscillateMatter(int nuFlavor, int nonOscNuFlavor, float Energy,
			float L, float dm2, float theta23, float UE32,
			float delta=0,int hierarchy = 1);
  float OscillateMatter(NtpMCTruth *mcth,float L, float dm2, 
			float theta23, float UE32,
			float delta=0,int hierarchy = 1);
  float OscillateMatter(ANtpTruthInfoBeam *ib, float L, 
			float dm2, float theta23, float UE32,
			float delta=0,int hierarchy = 1);
}


#endif //NUERWHELPERS_H
