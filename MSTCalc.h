#ifndef MSTCALC_H
#define MSTCALC_H
#include "TObject.h"


class MSTCalc: public TObject
{

public:
   MSTCalc();
   MSTCalc(const MSTCalc &mst);
   virtual ~MSTCalc();

   virtual void Print(Option_t* option="") const;
   void Reset();
   void Clear(Option_t* option = "");

   int enmsg;
   float ew1;
   int enn1;
   float esm1;
   float ewtot;
   int enntot;
   float esmtot;
   float e4w;
   float e4sm;
   int e4nn;
   float eb1;
   float eb25;
   int enbranch;

   int onmsg;
   float ow1;
   int onn1;
   float osm1;
   float owtot;
   int onntot;
   float osmtot;
   float o4w;
   float o4sm;
   int o4nn;
   float ob1;
   float ob25;
   int onbranch;

   float *eallw1;//[enn1]
   float *oallw1;//[onn1]

   float *eallm1;//[enn1]
   float *oallm1;//[onn1]

   double eeprob;
   double oeprob;
   double ealpha;
   double oalpha;

   double ebeta;
   double obeta;


private:
   ClassDef(MSTCalc,2)
};

#endif //MSTCALC_H
