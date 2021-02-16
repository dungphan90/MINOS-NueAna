#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"

ClassImp(Particle3D)


Particle3D::Particle3D()
{

	start_u=0.0;
	end_u=0.0;
	start_v=0.0;
	end_v=0.0;
	start_z=10000.0;
	end_z=0.0;
	particletype=Particle3D::other;
	sum_e=0.0;
	muonfrac=0.0;
	entries=0;
	numshared=0;
	theta=0;
	phi=0;
	emfit_a=0;
	emfit_b=0;
	emfit_e0=0;
	emfit_a_err=0;
	emfit_b_err=0;
	emfit_e0_err=0;
	emfit_prob=0;
	emfit_chisq=0;
	emfit_ndf=0;
	calibrated_energy=0;
	muonlike=0;
	avg_rms_t=0;





	u.clear();
	v.clear();
	z.clear();
	e.clear();
	chain.clear();
	chainhit.clear();
	shared.clear();
	lockshared.clear();
	types.clear();
	rms_t.clear();
	
	muon_threshold_max=2.5;
	muon_threshold_min=0.01;
	muon_min_hits=2; 
	
	cmp_chisq=0;
	cmp_ndf = 0;
	peakdiff=0;	
		
	pred_e_a=0;
	pred_g_a=0;
	pred_b=0;
	pred_e0=0;
	pred_e_chisq=0;
	pred_e_ndf=0;
	pred_g_chisq=0;
	pred_g_ndf=0;
	
	post_over=0;
	post_under=0;
	pre_over=0;
	pre_under=0;	
	
	
	pp_chisq=0;
	pp_ndf=0;
	pp_igood=0;
	pp_p=0;	
	
	
}

Particle3D::~Particle3D()
{
	u.clear();
	v.clear();
	z.clear();
	e.clear();
	chain.clear();
	chainhit.clear();
	shared.clear();
	lockshared.clear();
	types.clear();
	rms_t.clear();
}


void Particle3D::Clean()
{

//copy vector information
		std::vector<double>iu=u;
		std::vector<double>iv=v;		
		std::vector<double>iz=z;
		std::vector<double>ie=e;
		std::vector<int>ichain=chain;		
		std::vector<int>ichainhit=chainhit;	
		std::vector<int>iview=view;
		std::vector<double>irms_t=rms_t;
		
		std::vector<int>ishared=shared;  
		std::vector<int>ilockshared=lockshared; 

//clear existing information
	ResetHits();


//add back entries with >0 e
	
	for(unsigned int i=0;i<iu.size();i++)
	{
		if(ie[i]<0.01)continue;
		add_to_back( iu[i],  iv[i],  iz[i],  ie[i], ichain[i],  ichainhit[i], iview[i],  irms_t[i]);
		shared[shared.size()-1]=ishared[i];
		lockshared[lockshared.size()-1]=ilockshared[i];
		if(ishared[i])numshared++;
	
	}
	finalize();
	
//	printf("cleaned particle now has %d entries\n",entries);
}

void Particle3D::ResetHits()
{
	u.clear();
	v.clear();
	z.clear();
	e.clear();
	chain.clear();
	chainhit.clear();
	shared.clear();
	lockshared.clear();
	rms_t.clear();
	view.clear();
	
 	start_u=0;
 	end_u=0;
 	start_v =0;
 	end_v =0;
 	start_z = 10000.0;
 	end_z = 0;
 
  	sum_e=0.0;
  	muonfrac=0.0;
  	entries=0;
  	muonlike=0;
  	avg_rms_t=0;
  	numshared=0;
  	
	
	pred_e_a=0;
	pred_g_a=0;
	pred_b=0;
	pred_e0=0;
	pred_e_chisq=0;
	pred_e_ndf=0;
	pred_g_chisq=0;
	pred_g_ndf=0;	
	
	cmp_chisq=0;
	cmp_ndf = 0;
	peakdiff=0;	
		
		
		
	emfit_a=0;
	emfit_b=0;
	emfit_e0=0;
	emfit_a_err=0;
	emfit_b_err=0;
	emfit_e0_err=0;
	emfit_prob=0;
	emfit_chisq=0;
	emfit_ndf=0;
		
			  	
}

void Particle3D::SetShared(int ichain, int ichainhit)
{
	//find it
	int idx = -1;
	for(unsigned int i=0;i<z.size();i++)
	{
		if(chain[i]==ichain && chainhit[i]==ichainhit)
		{
			idx=i;
			break;
		}
	
	}
	if(idx<0)return;
	if(shared[idx]==0)numshared++;
	shared[idx]=1;

}

void Particle3D::UnsetShared(int ichain, int ichainhit)
{
	//find it
	int idx = -1;
	for(unsigned int i=0;i<z.size();i++)
	{
		if(chain[i]==ichain && chainhit[i]==ichainhit)
		{
			idx=i;
			break;
		}
	
	}
	if(idx<0)return;
	if(shared[idx]==1)numshared--;
	shared[idx]=0;


}




double Particle3D::GetNextEnergy(int ichain,int ichainhit,int inview)
{
	double energy =0;
	
	//find it
	int idx = -1;
	for(unsigned int i=0;i<z.size();i++)
	{
		if(chain[i]==ichain && chainhit[i]==ichainhit)
		{
			idx=i;
			break;
		}
	
	}
	
	//first look in other view
	for(unsigned int i=idx+1;i<z.size();i++)
	{
		if(view[i]==inview)
		{
			energy=e[i];
			break;
		}
	}
	
	inview=inview == 2 ? 3 :2;
	//then look in this view
	if(energy==0)
	for(unsigned int i=idx+1;i<z.size();i++)
	{
		if(view[i]==inview)
		{
			energy=e[i];
			break;
		}
	}
	
	
	
	return energy;
}

double Particle3D::GetPreviousEnergy(int ichain,int ichainhit,int inview)
{
	double energy =0;
	
	//find it
	int idx = -1;
	for(unsigned int i=0;i<z.size();i++)
	{
		if(chain[i]==ichain && chainhit[i]==ichainhit)
		{
			idx=i;
			break;
		}
	
	}
	
	for(int i=idx-1;i>-1;i--)
	{
		if(view[i]==inview)
		{
			energy=e[i];
			break;
		}
	}
	
	inview=inview == 2 ? 3 :2;
	//then look in this view
	if(energy==0)
		for(int i=idx-1;i>-1;i--)
	{
		if(view[i]==inview)
		{
			energy=e[i];
			break;
		}
	}
	
	return energy;
}



void Particle3D::SetEnergy(int ichain,int ichainhit,int inview, double energy)
{
	
	//find it
	for(unsigned int i=0;i<z.size();i++)
	{
		if(chain[i]==ichain && chainhit[i]==ichainhit && view[i]==inview)
		{
			e[i]=energy;
			return;
		}
	
	}
}

int Particle3D::ShareLocked(int ichain, int ichainhit)
{
	for(unsigned int i=0;i<z.size();i++)
	{
		if(chain[i]==ichain && chainhit[i]==ichainhit)
		{
			return lockshared[i];
			
		}
	
	}
	return -1;
}


void Particle3D::finalize()
{
	double etot=0;
	double rms_t_e=0;
	
	muonlike=0;
	for(int i=0;i<entries;i++)
	{
		rms_t_e+=e[i]*rms_t[i];
		etot+=e[i];
		
		if(e[i] < muon_threshold_max && e[i] > muon_threshold_min)
		{
			muonlike++;
		}
	
		
	}
	
	if (muon_min_hits <= entries ) muonfrac = (double)muonlike / (double)entries;
	
	if(etot>0)
	avg_rms_t = rms_t_e/etot;
	sum_e=etot;

}

void Particle3D::add_to_back(double iu, double iv, double iz, double ie,int ichain, int ichainhit,int iview, double irms_t)
{

	sum_e+=ie;
	
	
	entries++;

	u.push_back(iu);
	v.push_back(iv);
	z.push_back(iz);
	e.push_back(ie);
	shared.push_back(0);
	lockshared.push_back(0);
	chain.push_back(ichain);
	chainhit.push_back(ichainhit);
	view.push_back(iview);
	rms_t.push_back(irms_t);
	
	if(ie < muon_threshold_max && ie > muon_threshold_min)
	{
		muonlike++;
		if (muon_min_hits <= entries ) muonfrac = (double)muonlike / (double)entries;
	}
	

	
	if(start_z > iz)
	{
		start_z=iz;
		start_u=iu;
		start_v=iv;

	}
	
	if(end_z < iz)
	{
		end_z=iz;
		end_u=iu;
		end_v=iv;

	}


}
