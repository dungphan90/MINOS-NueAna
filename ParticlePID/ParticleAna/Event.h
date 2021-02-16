#ifndef EVENT_H
#define EVENT_H

#include "TObject.h"

class Event : public TObject
{

	public:
	
		Event();
		virtual ~Event();
		virtual void Clear(Option_t* option = "");

		double vtx_u;
		double vtx_v;
		double vtx_z;
		
		double sr_vtx_u; //!
		double sr_vtx_v; //!
		double sr_vtx_z; //!
		
		double min_u;
		double min_v;
		double min_z;
		double max_u;
		double max_v;
		double max_z;
		
		double visenergy;		


		double large_minz;
		double large_minu;
		double large_minv;
		double large_maxz;
		double large_maxu;
		double large_maxv;	
		
		double unused_e; //!
		double unused_e_avg; //!
		double unused_e_rms; //!	
		double pidA;
		double pidB;
		double pidC;
		double pidD;
		double pidE;
		double pidF;
		
		int unused_strips; //!
		int inFiducial; 
		int contained;
		int foundlongmuon;
		int foundprimaryshower;
		int nstrips;
		int nclusters;

		

	private:
		void Init();

		
		ClassDef(Event,3)

};

#endif

