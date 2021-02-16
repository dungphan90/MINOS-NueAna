#ifndef MANAGEDHIT_H
#define MANAGEDHIT_H

#include "TObject.h"

namespace Managed
{

class ManagedHit : public TObject
{
	public:
		ManagedHit(int view,int plane, int strip, double z,double t,double e);
		ManagedHit(){id=idcounter++;};
		virtual ~ManagedHit();

		static int idcounter;
		int id;
		static void ResetIDCounter();
		void AdvanceID();

		//to take energy from a hit... returns how much is actually removed (if we ask for too much, we return only what is available)
		double TakeEnergy(double e){ if(e<e_remaining)e_remaining-=e;else{e=e_remaining;e_remaining=0;}return e;};
		
		double SetEnergy(double e){ if(e>e_remaining)return e_remaining; else e_remaining=e;return e_remaining;};
		double GetT(){return t;};
		double GetZ(){return z;};
		int GetView(){return view;};
		int   GetStrip(){return strip;};
		int	  GetPlane(){return plane;};
		double GetEOriginal(){return e_original;};
		double GetERemaining(){return e_remaining;};
		int GetID(){return id;};
		
	private:
		double t;
		double z;
		int view;
		int plane;
		int strip;
		double e_original;
		double e_remaining;
		
	ClassDef(ManagedHit,1);

};


}

#endif

