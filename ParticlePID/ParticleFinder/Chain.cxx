#include "NueAna/ParticlePID/ParticleFinder/Chain.h"
#include "MessageService/MsgService.h"

#include <math.h>

CVSID("$Id: Chain.cxx,v 1.3 2014/02/18 03:34:30 rhatcher Exp $");

ClassImp(Chain)

int Chain::lastid = 0;

Chain::Chain()
{

	start_t=0.0;
	end_t=0.0;
	start_z=0.0;
	end_z=0.0;
	particletype=0;
	sum_e=0.0;
	weighted_t=0.0;
	entries=0;
	avg_slope=0.0;
	avg_offset=0.0;
	last_slope=0.0;
	front_slope=0.0;
	back_slope=0.0;
	front_offset=0.0;
	back_offset=0.0;
	parentChain=-1;
	myId=-1;
	level=0;
	muonfrac=0.0;
	interior_muonfrac=0.0;
	emlike=0.0;
	num_peaks=0;
	strict_decreasing=1;
	strict_increasing=1;
        lastpeake=0.0;
	lastpeakprob=0;
	sum_z=0.0;
	sum_t=0.0;
	sum_z_z=0.0;
	sum_z_t=0.0;
	muonlike=0;
	a=0.0;
	b=0.0;




	t.clear();
	z.clear();
	e.clear();
	cluster_id.clear();
	children.clear();

	myId=lastid++;
	
	
	muon_threshold_max=2;
	muon_threshold_min=0.5;
	muon_min_hits=2; 
	
	available=1;
	
}

Chain::~Chain()
{
	t.clear();
	z.clear();
	e.clear();
	cluster_id.clear();
	children.clear();
}

void Chain::ResetCounter()
{
	lastid=0;
}


//interpolate the front of the chain....
double Chain::interpolate(double test_z)
{
	double int_t=0;
	if(z.size()<1)return 0;
/*
	printf("nentries %d\n",z.size());
	for(unsigned int i=0;i<z.size();i++)
	{
		printf("zt %f %f\n",z[i],t[i]);
	
	}
*/	
	
	//hack to remove entries with the same z position.... we should really be merging clusters.....
	
	std::vector<double>z;
	std::vector<double>t;
	std::vector<double>e;


/*
	printf("raw nentries %d\n",this->z.size());
	for(unsigned int i=0;i<this->z.size();i++)
	{
		printf("zt %f %f\n",this->z[i],this->t[i]);
	
	}	*/
	
	
	z.push_back(this->z[0]);
	t.push_back(this->t[0]*this->e[0]);
	e.push_back(this->e[0]);
	int cnt=0;
	for(unsigned int i=1;i<this->z.size();i++)
	{
		if(fabs(this->z[i]-this->z[i-1])>0.001)
		{
			cnt++;		
			z.push_back(this->z[i]);
			t.push_back(this->t[i]*this->e[i]);
			e.push_back(this->e[i]);	
		}else{
			t[cnt]+=this->t[i]*this->e[i];
			e[cnt]+=this->e[i];
		}	
	}	
	
	for(unsigned int i=0;i<z.size();i++)
	{	
		t[i]/=e[i];
	}

/*	
	printf("clean nentries %d\n",z.size());
	for(unsigned int i=0;i<z.size();i++)
	{
		printf("zt %f %f\n",z[i],t[i]);
	
	}	*/

	

	if(z.size()>2)
	{	
	
		if(test_z>=end_z)
		{
			int i=z.size()-1;
			if(i>1)
				int_t = interpolate(z[i],t[i],z[i-1],t[i-1],z[i-2],t[i-2],test_z);
		}else if(test_z<=start_z)
		{
			if(z.size()>2)
				int_t = interpolate(z[0],t[0],z[1],t[1],z[2],t[2],test_z);
		}else{
			int i=z.size()-1;
			if(i>1)
				int_t = interpolate(z[i],t[i],z[i-1],t[i-1],z[i-2],t[i-2],test_z);
			if(z.size()>2)
				int_t += interpolate(z[0],t[0],z[1],t[1],z[2],t[2],test_z);
				
			int_t /=2;	
		}

	}else if(z.size()>1)
	{
		int end = z.size()-1;
		double dz = z[end]-z[end-1];
		double dt = t[end]-t[end-1];
		if(dz!=0)
		{
			double slope = dt/dz;
			int_t = (test_z-z[end-1])*slope+t[end-1];
		}
	}else if(z.size()==1)
	{
		int_t = t[0];
	}


	

	return int_t;
}


double Chain::interpolate(double z0, double t0, double z1, double t1, double z2, double t2, double z)
{
	double int_t=0;
  
    //for quadratic	
	double L0 = (z-z1)*(z-z2)/( (z0-z1)*(z0-z2));
	double L1 = (z-z0)*(z-z2)/( (z1-z0)*(z1-z2));	
	double L2 = (z-z0)*(z-z1)/( (z2-z0)*(z2-z1));

	int_t = t0*L0 + t1*L1 + t2*L2;

//	printf("quad inter L0 %f L1 %f L2 %f   z0 %f z1 %f z2 %f  t0 %f t1 %f t2 %f\n",L0,L1,L2,z0,z1,z2,t0,t1,t2);

	return int_t;
}











void Chain::Recalc()
{

		std::vector<double>mt=t;
		std::vector<double>mz=z;
		std::vector<double>me=e;
		std::vector<int>mc=cluster_id;
		
		ClearHits();
		
		for(unsigned int i=0;i<mt.size();i++)
			if(me[i]>0)insert(mt[i],mz[i],me[i],cluster_id[i]);



}


//cause the chain to reverse itself... 
//care must be taken with parent, children, which are not effected by this function
void Chain::Reverse()
{

	//printf("reversing\n");

		std::vector<double>mt;
		for(unsigned int i=0;i<t.size();i++)mt.push_back(t[i]);
		std::vector<double>mz;
		for(unsigned int i=0;i<z.size();i++)mz.push_back(z[i]);
		std::vector<double>me;
		for(unsigned int i=0;i<e.size();i++)me.push_back(e[i]);
		
		std::vector<int>mc;
		for(unsigned int i=0;i<cluster_id.size();i++)mc.push_back(cluster_id[i]);
		
		//printf("sv\n");
		
		ClearHits();
		
		//printf("cleared hits\n");
		
		for(int i=mt.size()-1;i>-1;i--)
			add_to_back(mt[i],mz[i],me[i],cluster_id[i]);


		//printf("done\n");

}

//for clearing information about hits....
//does not reset id, parent, or children information
void Chain::ClearHits()
{
 	start_t=0.0;
 	end_t=0.0;
 	start_z=0.0;
 	end_z=0.0;
 	sum_e=0.0;
	weighted_t=0.0; 
 	entries=0;
 	avg_slope=0.0;
 	avg_offset=0.0;
 	last_slope=0.0;
 	front_slope=0.0;
 	back_slope=0.0;
 	muonfrac=0.0;
 	interior_muonfrac=0.0;
 	emlike=0.0;
 	num_peaks=0;
 	strict_decreasing=1;
 	strict_increasing=1;
 	lastpeake=0.0;
 	lastpeakprob=0;
 	sum_z=0.0;
 	sum_t=0.0;
 	sum_z_z=0.0;
 	sum_z_t=0.0;
	muonlike=0;
	a=0.0;
    b=0.0;


	z.clear();
	e.clear();
	t.clear();
	cluster_id.clear();

}


void Chain::PrintChain()
{

	if(!MsgService::Instance()->IsActive("Chain",Msg::kDebug))return;

	printf("Chain %d \n",myId);
	printf("%d entries    start_t %f start_z %f  end_t %f end_z %f\n",entries,start_t,start_z,end_t,end_z);
	printf("sume %f   avg_slope %f lastslope %f  frontslope %f backslope %f\n",sum_e,avg_slope,last_slope,front_slope,back_slope);
	printf("muonfrac %f, intmf %f, emlike %f, peaks %d, sdec %d, sinc %d\n",muonfrac, interior_muonfrac, emlike, num_peaks,strict_decreasing, strict_increasing);
	
	printf(" hits (t,z,e) : ");
	for(unsigned int i=0;i<t.size();i++)
		printf(" (%f, %f, %f - %d)", t[i],z[i],e[i],cluster_id[i]);
	printf("\n");


}


void Chain::updateMuonFrac(double /*it*/, double /*iz*/, double ie)
{
	if(ie < muon_threshold_max && ie > muon_threshold_min)
	{
	
		if(entries>2)
		{
			int cntfirst=0;
			if ( e[0] <  muon_threshold_max && e[0] > muon_threshold_min ) cntfirst=1;
			interior_muonfrac = (double)(muonlike - cntfirst) / (double)(entries-2);
		}
		
		muonlike++;
		if (muon_min_hits <= entries ) muonfrac = (double)muonlike / (double)entries;
	}

}

void Chain::updateNumPeaks(double /*it*/, double /*iz*/, double ie)
{
  if (entries>1) {
	if (ie < lastpeake && lastpeakprob) {
		num_peaks++;
		lastpeakprob = 0;
	} else if (ie > lastpeake) {
      lastpeakprob = 1;
	}

    if(ie <= lastpeake) strict_increasing=0;
    if(ie >= lastpeake) strict_decreasing=0;
  }
	
  lastpeake = ie;

}

void Chain::updateEMLike(double /*it*/, double /*iz*/, double /*ie*/)
{
	//for now its emlike if it has 1 peak.... make it better later...
	emlike = num_peaks==1;
}


int Chain::good_slope()
{
	//tell us if we have enough points to believe the slope value
	
	int good=0;
	
	if(t.size()>1)good=1;
	
	return good;
}


//insert a cluster in the proper position, keeping the beginning upstream
void Chain::insert(double it, double iz, double ie, int my_cluster_id)
{

	double lastt = 0.0;
	double lastz = 0.0;
	double laste = 0.0;

	if(t.size()>0)
	{
		lastt=t.back();
		lastz=z.back();
		laste=e.back();
	}
	
	int insert_idx=-1;
	for(unsigned int i=0;i<t.size();i++)
	{
		if(iz<z[i])
		{
			insert_idx=i;
			break;
		}
	}
	
	if(t.size()==0)
	{
		start_t=it;
		start_z=iz;
		end_t=it;
		end_z=iz;	
	}	
/*	
	if(insert_idx==0)
	{
		start_t=it;
		start_z=iz;
	}
	
	if(insert_idx==t.size()-1)
	{
		end_t=it;
		end_z=iz;
	}
*/

	if(start_z>iz)
	{
		start_z=iz;
		start_t=it;
	}
	
	if(end_z<iz)
	{
		end_z=iz;
		end_t=it;
	}
	
	weighted_t*=sum_e;
	weighted_t +=it*ie;
	
	sum_e+=ie;
	
	weighted_t /=sum_e;
	
	entries++;
	if(insert_idx>-1)
	{
		std::vector<double>::iterator itr_f=t.begin();
		itr_f+=insert_idx;
		t.insert(itr_f,it);
		itr_f=z.begin();
		itr_f+=insert_idx;
		z.insert(itr_f,iz);
		itr_f=e.begin();
		itr_f+=insert_idx;
		e.insert(itr_f,ie);
		std::vector<int>::iterator itr_i=cluster_id.begin();
		itr_i+=insert_idx;
		cluster_id.insert(itr_i,my_cluster_id);
	}else{
		t.push_back(it);
		z.push_back(iz);
		e.push_back(ie);
		cluster_id.push_back(my_cluster_id);
	}
	//printf("added %f %f %f %d\n",it,iz,ie,my_cluster_id);
		
	sum_z+=iz;
	sum_t+=it;
	sum_z_z+=iz*iz;
	sum_z_t+=iz*it;
	
		
	double n=(double)entries;
	if(n>1 && (n*sum_z_z-(sum_z)*(sum_z))!=0) 
	{
		a=(sum_t*sum_z_z-sum_z*sum_z_t)/(n*sum_z_z-(sum_z)*(sum_z));
		b=(n*sum_z_t-sum_z*sum_t)/(n*sum_z_z-(sum_z)*(sum_z));	
	}
	
	avg_slope = b;
	avg_offset = a;
	
	if(entries>1 && lastz-iz)last_slope =(lastt-it)/(lastz-iz);
	
	
	if(entries>1 && entries<=4)
	{
		front_slope=b;
		front_offset=a;
		back_slope=b;
		back_offset=a;
	}
	if(entries>4)
	{
		//back slope
		int n=4;
		double sum_z_t=0;
		double sum_z=0;
		double sum_t=0;
		double sum_z_z=0;
		
		for(unsigned int i=e.size()-1;i>e.size()-5;i--)
		{
			sum_z+=z[i];
			sum_t+=t[i];
			sum_z_t+=z[i]*t[i];
			sum_z_z+=z[i]*z[i];
		}
		
	
		if(fabs(n*sum_z_z-(sum_z)*(sum_z))>1e-10)
		{	
			back_slope=(n*sum_z_t-sum_z*sum_t)/(n*sum_z_z-(sum_z)*(sum_z));
			back_offset=(sum_t*sum_z_z-sum_z*sum_z_t)/(n*sum_z_z-(sum_z)*(sum_z));
		}
		
		//front slope
		n=4;
		if(entries>12)n=entries/3;
		if(n>10)n=10;
		sum_z_t=0;
		sum_z=0;
		sum_t=0;
		sum_z_z=0;
		
		for(int i=0;i<n;i++)
		{
			sum_z+=z[i];
			sum_t+=t[i];
			sum_z_t+=z[i]*t[i];
			sum_z_z+=z[i]*z[i];
		}
		
		if(fabs(n*sum_z_z-(sum_z)*(sum_z))>1e-10)
		{	
			front_slope=(n*sum_z_t-sum_z*sum_t)/(n*sum_z_z-(sum_z)*(sum_z));
			front_offset=(sum_t*sum_z_z-sum_z*sum_z_t)/(n*sum_z_z-(sum_z)*(sum_z));		
		}
	}
	
		
	updateMuonFrac(it, iz, ie);
	updateNumPeaks(it, iz, ie);
	updateEMLike(it, iz, ie);

}


void Chain::add_to_back(double it, double iz, double ie, int my_cluster_id)
{

	double lastt = 0.0;
	double lastz = 0.0;
	double laste = 0.0;

	if(t.size()>0)
	{
		lastt=t.back();
		lastz=z.back();
		laste=e.back();
	}else{
		start_t=it;
		start_z=iz;
	}
	
	weighted_t*=sum_e;
	weighted_t +=it*ie;
	
	sum_e+=ie;
	
	weighted_t /=sum_e;
	
	entries++;

	t.push_back(it);
	z.push_back(iz);
	e.push_back(ie);
	cluster_id.push_back(my_cluster_id);
	
	//printf("added %f %f %f %d\n",it,iz,ie,my_cluster_id);
		
	sum_z+=iz;
	sum_t+=it;
	sum_z_z+=iz*iz;
	sum_z_t+=iz*it;
	
	end_t=it;
	end_z=iz;
	
	
	double n=(double)entries;
	if(n>1 && fabs(n*sum_z_z-(sum_z)*(sum_z))>1e-10 )
	{
		a=(sum_t*sum_z_z-sum_z*sum_z_t)/(n*sum_z_z-(sum_z)*(sum_z));
		b=(n*sum_z_t-sum_z*sum_t)/(n*sum_z_z-(sum_z)*(sum_z));	
	}
	
	avg_slope = b;
	avg_offset = a;
	
	if(entries>1 && lastz-iz)last_slope =(lastt-it)/(lastz-iz);
	
	
	if(entries>1 && entries<=4)
	{
		front_slope=b;
		front_offset=a;
		back_slope=b;
		back_offset=a;
	}
	if(entries>4)
	{
		int n=4;
		double sum_z_t=0;
		double sum_z=0;
		double sum_t=0;
		double sum_z_z=0;
		
		for(unsigned int i=e.size()-1;i>e.size()-5;i--)
		{
			sum_z+=z[i];
			sum_t+=t[i];
			sum_z_t+=z[i]*t[i];
			sum_z_z+=z[i]*z[i];
		}
		
		if(fabs(n*sum_z_z-(sum_z)*(sum_z))>1e-10)
		{	
			back_slope=(n*sum_z_t-sum_z*sum_t)/(n*sum_z_z-(sum_z)*(sum_z));
			back_offset=(sum_t*sum_z_z-sum_z*sum_z_t)/(n*sum_z_z-(sum_z)*(sum_z));
		}
	}
	
	
	updateMuonFrac(it, iz, ie);
	updateNumPeaks(it, iz, ie);
	updateEMLike(it, iz, ie);

}
