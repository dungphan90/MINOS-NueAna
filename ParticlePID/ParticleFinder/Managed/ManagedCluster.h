#ifndef MANAGEDCLUSTER_H
#define MANAGEDCLUSTER_H

#include "TObject.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedHit.h"

#include <vector>
#include <map>

namespace Managed
{

class ManagedCluster : public TObject
{
	public:

		ManagedCluster();
		//ManagedCluster(size_t a){};
		virtual ~ManagedCluster();


		static int idcounter;
		int id;

		double z;
		double dz;
		double t;
		double dt;
		int view;
	
		double tmin;
		double tmax;
		double zmin;
		double zmax;
	
		int inuse;
	
		double rms_t;
		double e;
	
		std::vector<int> hitplane;
		std::vector<int> hitstrip;
		std::vector<double> hite;
		std::vector<double> hitz;	
		std::vector<double> hitt;
		std::vector<int> hit_id;

		void Reset();
		void AdvanceID();
		void Insert(double z,double t,double energy,int plane,int strip,int hit_id);
		void Finalize();

		void ResetIDCounter();
		
		std::map<double,int> tsortmap;

		int GetStatus(){return status;};
		void SetStatus(int status){this->status=status;};

	private:
		int status;
	ClassDef(ManagedCluster,1)

};

}

#endif

