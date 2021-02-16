#ifndef NUEPID_H
#define NUEPID_H

#include "Record/RecRecordImp.h"
#include "Record/RecHeader.h"
#include "NueAna/NuePIDHeader.h"

class NuePID : public RecRecordImp<NuePIDHeader>
{
public:
   NuePID();
   NuePID(const NuePIDHeader &pidhead);
   virtual ~NuePID();

   Int_t IsNue; //1 yes, -1 no, 0 ?
   Float_t likelihood;
private:
   ClassDef(NuePID,1)
      
};

#endif //NUEPID_H
