#ifndef StripHit_H
#define StripHit_H



#include "TObject.h"

class StripHit : public TObject
{

	public:
	
	StripHit();
	~StripHit();

	int view;
	int plane;
	int strip;
	double t;
	double z;
	double e;


		

     private:
        ClassDef(StripHit,1)		
};


#endif

