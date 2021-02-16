#ifndef EVENTQUALITY_H
#define EVENTQUALITY_H


#include "TObject.h"


#include <vector>




class EventQuality : public TObject

{

	public:
	
		EventQuality();
		virtual ~EventQuality();

		void Reset();
		
		//general quality check
		bool IsQuality();


		// if there is a long muon found in only one view, then the event is bad
		//this can be caused if the event is too close to the coil or the edge of the detector
		int single_view_long_muon;
		int single_view_primary_shower;

		int foundlongmuon;
		int foundprimaryshower;

     private:
        ClassDef(EventQuality,1)


};

#endif

