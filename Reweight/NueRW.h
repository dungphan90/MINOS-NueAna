#ifndef NUERW_H
#define NUERW_H

#include "TObject.h"
#include "Conventions/Detector.h"

class NueRW : public TObject
{
public:
   NueRW();
   virtual ~NueRW();

   typedef enum EFileType {
      kUnknown = 0,
      kBEAM = 1,
      kNUE = 2,
      kTAU = 3,
      kAGG = 4
   } FileType_t;

   static const char *AsString(FileType_t t);

   virtual std::ostream &Print(std::ostream &os) const;
   virtual void Print(Option_t *option="") const;
   void Reset();

   int FindEBin(float E);

   int fRun;
   int fSubRun;
   Detector::Detector_t fDet;
   FileType_t fFileType;
   int ftgt;

   float qel_ma;
   float res_ma;
   float coh_ma;
   float qel_fa0;
   float qel_eta;
   float res_omega;
   float res_z;
   float coh_r0;
   float coh_rei;
   float kno_a1;
   float kno_a2;
   float kno_a3;
   float kno_a4;
   float kno_b;
   float kno_r112;
   float kno_r122;
   float kno_r132;
   float kno_r142;
   float kno_r113;
   float kno_r123;
   float kno_r133;
   float kno_r143;

   float dm2;
   float ss2th;
   float UE32;

   int randrow;
   int nfiles;
   int nsnarls;
   int nevents;
   int neventswpid;
   int nacc;

   float nsig;
   float nbg;
   float nnueb;
   float nnumu;
   float nnutau;
   float nnc;

   int EBINS;
   float EBINW;

   float *nsigE;//[EBINS]
   float *nbgE;//[EBINS]
   float *nnuebE;//[EBINS]
   float *nnumuE;//[EBINS]
   float *nnutauE;//[EBINS]
   float *nncE;//[EBINS]


   
   const NueRW operator + (const NueRW &rw2) const {
      NueRW a;
      if(fDet!=rw2.fDet||randrow!=rw2.randrow){
	 a.randrow=-99;
	 return a;
      }
      a.fDet=fDet;
      a.randrow=randrow;
      a.fFileType=kAGG;

      a.qel_ma=qel_ma;
      a.res_ma=res_ma;
      a.coh_ma=coh_ma;
      a.qel_fa0=qel_fa0;
      a.qel_eta=qel_eta;
      a.res_omega=res_omega;
      a.res_z=res_z;
      a.coh_r0=coh_r0;
      a.coh_rei=coh_rei;
      a.kno_a1=kno_a1;
      a.kno_a2=kno_a2;
      a.kno_a3=kno_a3;
      a.kno_a4=kno_a4;
      a.kno_b=kno_b;
      a.kno_r112=kno_r112;
      a.kno_r122=kno_r122;
      a.kno_r132=kno_r132;
      a.kno_r142=kno_r142;
      a.kno_r113=kno_r113;
      a.kno_r123=kno_r123;
      a.kno_r133=kno_r133;
      a.kno_r143=kno_r143;
      a.dm2=dm2;
      a.ss2th=ss2th;
      a.UE32=UE32;

      a.nfiles=nfiles+rw2.nfiles;
      a.nsnarls=nsnarls+rw2.nsnarls;      
      a.nevents=nevents+rw2.nevents;
      a.neventswpid=neventswpid+rw2.neventswpid;
      a.nacc=nacc+rw2.nacc;
      a.nsig=nsig+rw2.nsig;
      a.nbg=nbg+rw2.nbg;
      a.nnueb=nnueb+rw2.nnueb;
      a.nnumu=nnumu+rw2.nnumu;
      a.nnutau=nnutau+rw2.nnutau;
      a.nnc=nnc+rw2.nnc;
      
      for(int i=0;i<EBINS;i++){
	a.nsigE[i]=nsigE[i]+rw2.nsigE[i];
	a.nbgE[i]=nbgE[i]+rw2.nbgE[i];
	a.nnuebE[i]=nnuebE[i]+rw2.nnuebE[i];
	a.nnumuE[i]=nnumuE[i]+rw2.nnumuE[i];
	a.nnutauE[i]=nnutauE[i]+rw2.nnutauE[i];
	a.nncE[i]=nncE[i]+rw2.nncE[i];
      }
      return a;
   }


   const NueRW operator / (const float s) const {
      NueRW a;
      a.fRun=fRun;
      a.fSubRun=fSubRun;
      a.fDet=fDet;
      a.fFileType=fFileType;
      a.ftgt=ftgt;
      a.qel_ma=qel_ma;
      a.res_ma=res_ma;
      a.coh_ma=coh_ma;
      a.qel_fa0=qel_fa0;
      a.qel_eta=qel_eta;
      a.res_omega=res_omega;
      a.res_z=res_z;
      a.coh_r0=coh_r0;
      a.coh_rei=coh_rei;
      a.kno_a1=kno_a1;
      a.kno_a2=kno_a2;
      a.kno_a3=kno_a3;
      a.kno_a4=kno_a4;
      a.kno_b=kno_b;
      a.kno_r112=kno_r112;
      a.kno_r122=kno_r122;
      a.kno_r132=kno_r132;
      a.kno_r142=kno_r142;
      a.kno_r113=kno_r113;
      a.kno_r123=kno_r123;
      a.kno_r133=kno_r133;
      a.kno_r143=kno_r143;      
      a.dm2=dm2;
      a.ss2th=ss2th;
      a.UE32=UE32;

      a.randrow=randrow;
      a.nfiles=nfiles;
      a.nsnarls=nsnarls;
      a.nevents=nevents;
      a.neventswpid=neventswpid;
      a.nacc=nacc;

      a.EBINS=EBINS;
      a.EBINW=EBINW;
      
      a.nsig=nsig/s;
      a.nbg=nbg/s;
      a.nnueb=nnueb/s;
      a.nnumu=nnumu/s;
      a.nnutau=nnutau/s;
      a.nnc=nnc/s;

      for(int i=0;i<EBINS;i++){
	a.nsigE[i]=nsigE[i]/s;
	a.nbgE[i]=nbgE[i]/s;
	a.nnuebE[i]=nnuebE[i]/s;
	a.nnumuE[i]=nnumuE[i]/s;
	a.nnutauE[i]=nnutauE[i]/s;
	a.nncE[i]=nncE[i]/s;
      }
      return a;
   }

   
      
private:
   ClassDef(NueRW,4)
      
};

#endif //NUERW_H
