#ifndef NueMiniAnaPID_H
#define NueMiniAnaPID_H                                                               

#include <vector>
#include <map>
#include "Conventions/Detector.h"
#include "Conventions/BeamType.h"
#include "Conventions/ReleaseType.h"

#include "NueAna/NueAnaTools/Selection.h"

#include "TObject.h"

using namespace std;                                   

class NueRecord;
class NueMiniPID;

class NueMiniAnaPID
{
   public:
     NueMiniAnaPID() {};
     virtual ~NueMiniAnaPID() {};

	 //need to speficy domrcc to determine which precord branch to use
     void FillMini(NueRecord *nr, NueMiniPID *nm, int domrcc=0);
     void FillRecord(NueRecord *nr, NueMiniPID *nm, int domrcc=0);

};

#endif 
