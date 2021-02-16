#ifndef HoughLine_h
#define HoughLine_h

#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedCluster.h"

#include <vector>
#include <math.h>
#include "TObject.h"
class HoughLine : public TObject
{
	public:
		HoughLine();
		~HoughLine();

		HoughLine(double theta, double r, double offset_t, double offset_z);
		void SetHoughParams(double theta, double r, double offset_t, double offset_z);
		
		double GetExpectedT(double z);

		void Reset();
		
		void ResetHits(int keep_bounds=0);
		
		void AddCluster(Managed::ManagedCluster *cl);
	
		std::vector<int> cluster_id;
	
		double sum_e;
		double start_z;
		double start_t;
		double end_z;
		double end_t;
		int ncluster;
		
		double theta;
		double r;
		double offset_t;
		double offset_z;
		double phi;
		double chi2;
 
 		int primary;
 
 	private:
 	
 		ClassDef(HoughLine,1);
};



/*


bool CompareChi2(const HoughLine & l,const HoughLine &r)
{
	return l.chi2/l.ncluster > r.chi2/r.ncluster;
}	


bool CompareT(const HoughLine & l,const HoughLine &r)
{
	if(fabs(l.start_t-r.start_t)<0.0001)
	{
		if(fabs(l.end_t-r.end_t)<0.0001)
		{
			return CompareChi2(l,r);
		}
		return l.end_t<r.end_t;
	}
	
	return l.start_t<r.start_t;
	
}


bool CompareLength(const HoughLine & l,const HoughLine &r)
{
	if(fabs(l.start_z-r.start_z)<0.0001)
	{
		if(fabs(l.end_z-r.end_z)<0.0001)
		{
			return CompareT(l,r);
		}

		return l.end_z<r.end_z;
	}
	
	return l.start_z<r.start_z;
	
}

bool CompareTotalEnergy(const HoughLine & l,const HoughLine &r)
{
	return l.sum_e > r.sum_e;
}

bool CompareForwardAndClusters(const HoughLine & l,const HoughLine &r)
{
	if(l.ncluster == r.ncluster)
	{
		//are they the same number of planes?
		if(fabs(l.end_z-l.start_z + r.start_z -r.end_z)/0.06 < 1)
		{		
			//return the most forward one
			return l.start_z>r.start_z;
		}
		//select the shortest one
		return l.end_z-l.start_z > r.end_z-r.start_z;
		
	}
	
	return l.ncluster < r.ncluster;
	
}
*/

extern "C" {
	bool CompareLength(const HoughLine & l,const HoughLine &r);
	bool CompareForwardAndClusters(const HoughLine & l,const HoughLine &r);
	bool CompareTotalEnergy(const HoughLine & l,const HoughLine &r);
	bool CompareChi2(const HoughLine & l,const HoughLine &r);
	bool CompareT(const HoughLine & l,const HoughLine &r);
}


#endif

