#ifndef NUEMINIANA_H
#define NUEMINIANA_H                                                               

#include <vector>
#include <map>
#include "Conventions/Detector.h"
#include "Conventions/BeamType.h"
#include "Conventions/ReleaseType.h"

#include "NueAna/NueAnaTools/Selection.h"

#include "TObject.h"

using namespace std;                                   

class NueRecord;
class NueMini;

class NueMiniAna
{
   public:
     NueMiniAna() {};
     virtual ~NueMiniAna() {};

     void FillMini(NueRecord *nr, NueMini *nm);
     void FillRecord(NueRecord *nr, NueMini *nm);

};

#endif 
