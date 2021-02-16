#include "NueAna/ParticlePID/ParticleFinder/HoughLine.h"
#include <math.h>

ClassImp(HoughLine)

bool operator < (const HoughLine & left, const HoughLine & right)
{
	if(left.ncluster < right.ncluster)return true;
	return false;
}


HoughLine::HoughLine()
{
	Reset();
}

HoughLine::~HoughLine()
{}

HoughLine::HoughLine(double theta, double r, double offset_t, double offset_z)
{
	SetHoughParams(theta, r, offset_t, offset_z);
}


void HoughLine::AddCluster(Managed::ManagedCluster *cl)
{

	//printf("z %f t %f\n",cl->z,cl->t);

	if(!cl)return;
	cluster_id.push_back(cl->id);
	
//	printf("adding cluster %d to hl %f %f  %f %f  ",cl->id,start_z,start_t,end_z,end_t);

	sum_e+=cl->e;
	if(ncluster==0 && start_z==0) //need to check start_z in case we did a resethits(1)
	{
		start_z=end_z=cl->z;
		start_t=end_t=cl->t;
	}

	if(start_z > cl->z)
	{
		start_z = cl->z;
		start_t = cl->t;
	}
	
	if(end_z < cl->z)
	{
		end_z = cl->z;
		end_t = cl->t;
	}
	
//	printf(" now at %f %f  %f %f\n",start_z,start_t,end_z,end_t);

	ncluster++;
	
	double exp_t = GetExpectedT(cl->z);
	chi2+=(cl->t-exp_t)*(cl->t-exp_t);
}


void HoughLine::SetHoughParams(double theta, double r, double offset_t, double offset_z)
{
	Reset();
	this->theta=theta;
	this->r=r;
	this->offset_t=offset_t;
	this->offset_z=offset_z;
	
	
	//determine the angle off of z "phi
	double d=1;
	double zpos1 = offset_z + (offset_t + r/sin(theta)-d)*sin(theta)/cos(theta);
	double zpos2 = offset_z + (offset_t + r/sin(theta)+d)*sin(theta)/cos(theta);
	
	phi = atan(2*d/(zpos2-zpos1));
	
}


void HoughLine::Reset()
{

		
	theta=0;
	r=0;
	offset_t=0;
	offset_z=0;
	chi2=0;
	phi=0;
	
	ResetHits(0);

}

void HoughLine::ResetHits(int keep_bounds)
{

	cluster_id.clear();
	
	if(!keep_bounds)
	{
		start_z=0;
		start_t=0;
		end_z=0;
		end_t=0;
	}
	
	
	ncluster=0;
	
	sum_e=0;
	primary=0;	
}

		
double HoughLine::GetExpectedT(double z)
{


	return (- cos(theta)/sin(theta))*(z-offset_z)+r/sin(theta) + offset_t;

}



extern bool CompareChi2(const HoughLine & l,const HoughLine &r)
{
        return l.chi2/l.ncluster > r.chi2/r.ncluster;
}       


extern bool CompareT(const HoughLine & l,const HoughLine &r)
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

extern bool CompareLength(const HoughLine & l,const HoughLine &r)
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

extern bool CompareTotalEnergy(const HoughLine & l,const HoughLine &r)
{
        return l.sum_e > r.sum_e;
}

extern bool CompareForwardAndClusters(const HoughLine & l,const HoughLine &r)
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







