#include "NueAna/ParticlePID/ParticleFinder/LongMuonFinder.h"
#include <cassert>
#include <map>
#include <math.h>
#include <algorithm>

#include "MessageService/MsgService.h"
#include "NueAna/ParticlePID/ParticleFinder/Chain.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedCluster.h"


#include <sstream>


CVSID("$Id: LongMuonFinder.cxx,v 1.3 2009/09/11 05:00:40 gmieg Exp $");

Particle3D * LongMuonFinder::foundparticle3d=0;
Chain * LongMuonFinder::chain_u=0;
Chain * LongMuonFinder::chain_v=0;


struct point
{
	int view;
	double z;
	double t;
	double e;
	int chain;
	int chainhit;
	double rms_t;
};
	
bool pointgreaterlmf(point p1, point p2){return p1.z < p2.z;}

//returns 0 if in a fully instrumented region
//returns 1 if in partial due to z
//returns 2 if in partial due to t (you need to look at other views in this case to make a determination)
int LongMuonFinder::IsPartiallyInstrumented(double t, double z, int view)
{
	if(detector==Detector::kFar)return 0;

	if(detector==Detector::kNear)
	{
		if(view==2)
		{
			if(z<0 ||z > 7)return 1;
			if(t<-0.2 || t>2.3)return 2;
			return 0;
		
		}else if(view==3)
		{
			if(z<0 ||z > 7)return 1;
			if(t<-2.3 || t>0.2)return 2;
			return 0;			
		}else
//		printf("unknown view in LongMuonFinder::IsPartiallyInstrumented\n");
		abort();
	}


//	printf("unknown detector in LongMuonFinder::IsPartiallyInstrumented\n");
	abort();
}

std::pair<int,int> LongMuonFinder::CountInPartiallyInstrumentedRegion(Chain *viewU, Chain *viewV)
{
	int ucount=0;
	int vcount=0;
	
	if(detector!=Detector::kNear)return std::pair<int,int>(ucount,vcount);

	std::map<double,std::pair<double,int> > zorder;// z, <t, view>
	
	for(unsigned int i=0;i<viewU->z.size();i++)
		zorder.insert(std::pair<double,std::pair<double,int> >(viewU->z[i],std::pair<double,int>(viewU->t[i],2)));
		
	for(unsigned int i=0;i<viewV->z.size();i++)
		zorder.insert(std::pair<double,std::pair<double,int> >(viewV->z[i],std::pair<double,int>(viewV->t[i],3)));

	

	//count the number in each plane that are in the partially instrumented region
	
	std::map<double,std::pair<double,int> >::iterator it;
	
	for(it=zorder.begin();it!=zorder.end();it++)
	{
		//values in the other views for extrapolation!
		double prevz=0;
		double prevt=0;
		double nextz=0;
		double nextt=0;

		double prevz2=0;
		double prevt2=0;
		double nextz2=0;
		double nextt2=0;
		
		int thisview = it->second.second;
		double thist = it->second.first;
		double thisz = it->first;
		
		if(thisview !=2 && thisview !=3)
		{
			printf("unknown view\n");
			abort();
		}
		
		std::map<double,std::pair<double,int> >::iterator itforward;
		std::map<double,std::pair<double,int> >::iterator itbackward;
		
		for(itforward=it;itforward!=zorder.end();itforward++)
		{
			if(itforward->second.first!=thisview )
			{
				if(!nextz)
				{
					nextz=itforward->first;
					nextt=itforward->second.second;
				}else if(fabs(itforward->first - nextz)>0.04){  //make sure we are in a different plane
					nextz2=itforward->first;
					nextt2=itforward->second.second;				
					break;
				}
			}	
		}
		for(itbackward=it;itbackward!=zorder.begin();itbackward--)
		{
			if(itbackward->second.first!=thisview)
			{
				if(!nextz)
				{
					prevz=itbackward->first;
					prevt=itbackward->second.second;
				}else if(fabs(itbackward->first - prevz)>0.04){  //make sure we are in a different plane
					prevz2=itbackward->first;
					prevt2=itbackward->second.second;				
					break;
				}
			}
		}
		
	//	printf("n %f %f p %f %f\n",nextz,nextz2,prevz,prevz2);
		
		if((nextz?1:0)+ (nextz2?1:0)+ (prevz?1:0)+ (prevz2?1:0) <2)break;//we can't do anything with information in only one view! need at least two points in the other view!
		
		//we need an estimation of the other view t position at the point in question...
		
		double est_other_t=0;
		
		double slope=0;
		double offset=0;
		
		if(nextz && prevz)
		{
			slope = (nextt-prevt)/(nextz-prevz);
			offset = prevt-slope*prevz;
		}else if(nextz && nextz2)
		{
			slope = (nextt-nextt2)/(nextz-nextz2);
			offset = nextt2-slope*nextz2;
		}else if(prevz && prevz2)
		{
			slope = (prevt-prevt2)/(prevz-prevz2);
			offset = prevt2-slope*prevz2;
		}
		
		est_other_t = thisz*slope+offset;
		
	//	printf("checking view %d z %f thist %f othert %f\n",thisview,thisz,thist,est_other_t);
		
		int partialviewthis = IsPartiallyInstrumented(thist,thisz,thisview);
		int partialviewother = IsPartiallyInstrumented(est_other_t,thisz,thisview==2?3:2);
		
		if(partialviewthis || partialviewother)
		{
			if(thisview==2)ucount++;
			else vcount++;
		}
		
	}
	
	
	
	
	return std::pair<int,int>(ucount,vcount);
	
				

}



LongMuonFinder::LongMuonFinder(Detector::Detector_t d)
{
	detector=d;
	Reset();
}

LongMuonFinder::~LongMuonFinder()
{}

void LongMuonFinder::Reset()
{
	if(foundparticle3d){delete foundparticle3d; foundparticle3d=0;}
	foundparticle=0;
	if(chain_u){delete chain_u; chain_u=0;}
	if(chain_v){delete chain_v; chain_v=0;}
	single_view_long_muon=0;
}

int LongMuonFinder::FindLongMuon(Managed::ClusterManager *cm)
{

	MSG("LongMuonFinder",Msg::kDebug)<<"looking for long muons\n";
	Reset();

	cluster_manager=cm;
	
	//try to find a long muon chain in each view
	if(chain_u){delete chain_u;chain_u=0;}
	if(chain_v){delete chain_v;chain_v=0;}
	chain_u = FindMuonChain(cm,2);
	chain_v = FindMuonChain(cm,3);
	
	int qual_u=0;
	int qual_v=0;
	
	
	
	if(chain_u)chain_u->Recalc();
	if(chain_v)chain_v->Recalc();
	
	std::pair<int,int> partialcount = CountInPartiallyInstrumentedRegion(chain_u, chain_v);
	
	//printf("partial count %d %d\n",partialcount.first,partialcount.second);
	
	if(chain_u)qual_u=CheckChainQuality(chain_u,2,partialcount.first);
	if(chain_v)qual_v=CheckChainQuality(chain_v,3,partialcount.second);
	

	if(qual_u&&!qual_v)
	{
		if(chain_u->end_z - chain_u->start_z > 2)single_view_long_muon=1;
	}else if(qual_v&&!qual_u)
	{
		if(chain_v->end_z - chain_v->start_z > 2)single_view_long_muon=1;
	}
	
	//if we have a quality chain in each view... then we will be making a particle
	if(qual_u && qual_v)
	{
	
		MSG("LongMuonFinder",Msg::kDebug)<<"qqqqq "<<qual_u <<" "<< qual_v<<"\n";
		
	
		//look in both views to find at least 3u and 3v consecutive planes with only 1 strip... that is where we say the muon exits the rest of the shower/partices/etc
		double isolation_z = FindIsolationZ();
		

		int qual=CheckChainOverlap(chain_u, chain_v, isolation_z);
		if(qual)
		{
		
		ClearFrontVertex(chain_u, chain_v);

		if(isolation_z>0)
		{		
			//need to remove energy from clusters forward of that point which are too high for a muon track
			
			RemoveNonMuonEnergy(chain_u,2,isolation_z);
			RemoveNonMuonEnergy(chain_v,3,isolation_z);
			
		
			//need to absorb clusters along the chain which came off the muon (brehms, etc)
			AbsorbMuonClusters(chain_u,2,isolation_z);
			AbsorbMuonClusters(chain_v,3,isolation_z);
		}
	

		//save clusters involved...

		for(int i=0;i<chain_u->entries;i++)
		{
			int newid = cluster_manager->SaveCluster(chain_u->cluster_id[i],chain_u->e[i],1);
			Managed::ManagedCluster * newcluster = cluster_manager->GetSavedCluster(newid);
			if(chain_u->e[i]!=newcluster->e)
				MSG("LongMuonFinder",Msg::kDebug)<<"saving with different e "<<chain_u->e[i]<<" "<<newcluster->e<<"\n";
			chain_u->e[i]=newcluster->e;
			chain_u->cluster_id[i]=newid;
		}

		for(int i=0;i<chain_v->entries;i++)
		{
			int newid = cluster_manager->SaveCluster(chain_v->cluster_id[i],chain_v->e[i],1);
			Managed::ManagedCluster * newcluster = cluster_manager->GetSavedCluster(newid);
			if(chain_v->e[i]!=newcluster->e)
				MSG("LongMuonFinder",Msg::kDebug)<<"saving with different e "<<chain_v->e[i]<<" "<<newcluster->e<<"\n";
			chain_v->e[i]=newcluster->e;
			chain_v->cluster_id[i]=newid;
		}		
		
		
		
		
	
	
		MSG("LongMuonFinder",Msg::kDebug)<<"making particle\n";
		MakeParticle3D();
		if(foundparticle3d)
		{
			foundparticle=1;
		
			foundparticle3d->particletype=Particle3D::muon;
		
			MSG("LongMuonFinder",Msg::kDebug)<<"particle made\n";
			DumpParticle3D();
			
			chain_u->available=0;
			chain_v->available=0;
		}
		}
	}


	MSG("LongMuonFinder",Msg::kDebug)<<"printing long muon chains\nu\n";
	if(chain_u)chain_u->PrintChain();MSG("LongMuonFinder",Msg::kDebug)<<"v\n";
	if(chain_v)chain_v->PrintChain();
	MSG("LongMuonFinder",Msg::kDebug)<<"done\n";
	

	MSG("LongMuonFinder",Msg::kDebug)<<"done looking for long muons\n";
	return foundparticle;
}


//make sure that the chains both start within one plane of each other... 
//for now, simply do this by finding the plane below the most upstream hit
void LongMuonFinder::ClearFrontVertex(Chain * chain_u, Chain * chain_v)
{
	double start_z_u = chain_u->start_z;
	double start_z_v = chain_v->start_z;

	if(fabs(start_z_u-start_z_v)<0.04)return; //close enough
	
	double start_z = start_z_u < start_z_v ? start_z_v : start_z_u - 0.05;
	
	Chain * toclear =  start_z_u < start_z_v ? chain_u : chain_v;
	//now delete all hits up to start_z
	
	Chain tmp = *toclear;
	
	toclear->ClearHits();
	
	for(int i=0;i<tmp.entries;i++)
	{
		if(tmp.z[i]>start_z)
			toclear->insert(tmp.t[i], tmp.z[i], tmp.e[i], tmp.cluster_id[i]);
	}
	
}

//make sure that the chains actually overlap in z (so from the same particle)
//and that the chain is of sufficient length
int LongMuonFinder::CheckChainOverlap(Chain * chain_u, Chain * chain_v, double clear_z)
{

	MSG("LongMuonFinder",Msg::kDebug)<<"checking chain overlap clear_z at "<<clear_z<<"\n";
	MSG("LongMuonFinder",Msg::kDebug)<<"isolation distance u "<<chain_u->end_z-clear_z<<"\n";
	MSG("LongMuonFinder",Msg::kDebug)<<"isolation distance v "<<chain_v->end_z-clear_z<<"\n";


	//first check that the muon is sufficiently long past the isolation point!
	double require_isolation_distance=1.0;
	if(clear_z <1e-9)return 0; //no valid z clearance measured!
	if(chain_u->end_z-clear_z<require_isolation_distance && 
		fabs(chain_u->end_z-chain_u->start_z)*0.7 > clear_z)return 0;
	if(chain_v->end_z-clear_z<require_isolation_distance && 
		fabs(chain_v->end_z-chain_v->start_z)*0.7 > clear_z)return 0;
	
	//now make sure that there is sufficient overlap
	double require_overlap_distance=1.0;
	
	double start_z = chain_u->start_z > chain_v->start_z ? chain_u->start_z : chain_v->start_z;
	double end_z = chain_u->end_z < chain_v->end_z ? chain_u->end_z : chain_v->end_z;

	MSG("LongMuonFinder",Msg::kDebug)<<"overlap "<<end_z-start_z<<"\n";
	
	if(end_z-start_z < require_overlap_distance)return 0;
	return 1;
	
	

}

void LongMuonFinder::RemoveNonMuonEnergy(Chain *ch,int view,double past_z)
{


	if(!cluster_manager)return;
	
	std::map<double, std::map<double, int> > * cluster_map = cluster_manager->GetClusterMap(view);
    std::map<double, std::map<double, int> >::iterator p_iterr;
    std::map<double, int >::iterator s_iterr;

	//first determine the average energy deposition
	double sum_e=0;
	int cnt_e=0;
	for(p_iterr=cluster_map->begin();p_iterr!=cluster_map->end(); p_iterr++)
    {	
    	if(p_iterr->first-past_z<0.01)continue;
    
    	for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
        {
        	double e =cluster_manager->GetCluster(s_iterr->second)->e;
        	if(e<0.1)continue;
        	cnt_e++;
        	sum_e+=e;
        }
    }

	if(cnt_e<1)return;
	double avg_e=sum_e/(double)cnt_e;
	

	for(p_iterr=cluster_map->begin();p_iterr!=cluster_map->end(); p_iterr++)
    {	
    	if(p_iterr->first-past_z>0.01)return;
    	for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
        {
        	Managed::ManagedCluster *c = cluster_manager->GetCluster(s_iterr->second);
        	if(c->e < avg_e * 1.2)continue;
        	for(unsigned int i=0;i<ch->cluster_id.size();i++)
        	{
        		if(ch->cluster_id[i]==c->id)
        		{	
        			MSG("LongMuonFinder",Msg::kDebug)<<"changing cluster energy in chain from "<<ch->e[i]<<" to "<<avg_e<<"\n";
        			ch->e[i]=avg_e;

        			break;
        		}
        	}
			
		}
	}
	

}

void LongMuonFinder::AbsorbMuonClusters(Chain *ch,int view,double past_z)
{
	//for each cluster in the chain, merge it with all other clusters in the plane past z distance 

	if(!cluster_manager)return;
	
	std::map<double, std::map<double, int> > * cluster_map = cluster_manager->GetClusterMap(view);

    std::map<double, std::map<double, int> >::iterator p_iterr;
    std::map<double, int >::iterator s_iterr;

	std::vector<int> clusters;

	double last_z=0;
	for(p_iterr=cluster_map->begin();p_iterr!=cluster_map->end(); p_iterr++)
    {	
    	if(p_iterr->first-past_z<0.02)continue;
    		
    	if(fabs(last_z-p_iterr->first)>0.01)
		{
			if(clusters.size()>0)MergeChainClusters(ch,&clusters);
			clusters.clear();
			last_z=p_iterr->first;
		}
	
	
    	
    	for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
        {
			clusters.push_back(s_iterr->second);
		}
	}
	
	if(clusters.size()>0)MergeChainClusters(ch,&clusters);
	clusters.clear();
}


void LongMuonFinder::MergeChainClusters(Chain * ch, std::vector<int> *clusters)
{
	if(clusters->size()<1)return;

/*
	//find the appropriate spot in the chain for these clusters... 
	Managed::ManagedCluster *current=0;
	
	unsigned int idx;
	for(idx=0;idx<ch->z.size();idx++)
	{
		if(fabs(ch->z[idx]-(*clusters)[0]->z)<0.01)
		{
			current = cluster_manager->GetCluster(ch->cluster_id[idx]);
			break;
		}
	}
		
	int add_to_chain=0;	
	if(!current)
	{	
		add_to_chain=1;
		current=(*clusters)[0];
	}
*/

	
	int mergeid=(*clusters)[0];
	for(unsigned int i=1;i<(*clusters).size();i++)
	{
		mergeid=cluster_manager->MergeClusters(mergeid,(*clusters)[i]);
	}


	Managed::ManagedCluster *current = cluster_manager->GetCluster(mergeid);

	//do we already have a cluster for that plane?
	int add_to_chain=1;	
	int idx=-1;
	for(idx=0;idx<(int)ch->z.size();idx++)
	{
		if(fabs(ch->z[idx]-current->z)<0.01)
		{
			add_to_chain=0;
			break;
		}
	}


	if(add_to_chain)
	{
		ch->insert(current->t, current->z, current->e, current->id);
	}else{
		ch->e[idx]=current->e;
		//for our purposes... we just want to account from brehm like energy... we dont want to move the point of the muon track hit
	//	ch->t[idx]=current->t;
	//	ch->z[idx]=current->z;
		ch->cluster_id[idx]=current->id;
	}

}



double LongMuonFinder::FindIsolationZ()
{
	if(!cluster_manager)return 0;
	
	//require 3 consecutive u and v planes to signify start of muon isolation
	std::map<double, std::map<double, int> > * cluster_map = cluster_manager->GetClusterMap();


    std::map<double, std::map<double, int> >::iterator p_iterr;
    std::map<double, int >::iterator s_iterr;

	double start_z[6];
	for(int i=0;i<6;i++)start_z[i]=0;
	int ucount=0;
	int vcount=0;
	int cnt=0;
	int lastview=0;
	
	std::vector<double> goodz;
	
    for(p_iterr=cluster_map->begin();p_iterr!=cluster_map->end(); p_iterr++)
    {	
    	
		if(fabs(start_z[5]-p_iterr->first)>0.03)
		{
		
			
			for(int i=0;i<5;i++)start_z[i]=start_z[i+1];
			start_z[5]=p_iterr->first;
			
			if(cnt<=1)
			{
				if(lastview==2)ucount++;
				else vcount++;
				
				goodz.push_back(p_iterr->first);
			}else{
				ucount=0;
				vcount=0;
				goodz.clear();
			}
			
			MSG("LongMuonFinder",Msg::kDebug)<<"z "<<start_z[5] <<" lastview "<<lastview <<", ucount "<< ucount<<" vcount "<<vcount <<" pcnt "<< cnt<<"\n";
				
			cnt=0;
		}

    	if(ucount>=3 && vcount>=3)break;

    	for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
        {
        	Managed::ManagedCluster *this_cluster = cluster_manager->GetCluster(s_iterr->second);
			int view=this_cluster->view;
			
        	//did we have an empty plane?
        	//or do we have a very large, unmuon type hit?
        	if( /*(cnt==0 && view==lastview) ||*/ this_cluster->e >3)
        	{
				ucount=0;
				vcount=0;
				goodz.clear();
        	}
        	
        	
        
        	lastview=view;
        	cnt++;
		}
	
	}

	MSG("LongMuonFinder",Msg::kDebug)<<"!!!!! u "<<ucount<<" v "<<vcount<<" z "<<start_z[0]<<" good z "<<goodz[0]<<"\n";
	if(ucount>=3 && vcount>=3) return goodz[0]; //start_z[0]
	
	return 0;	

}



// check to see if this is a muon-like chain with a sufficient number of hits 
int LongMuonFinder::CheckChainQuality(Chain *c, int /*view*/, int partialcount)
{
	if(!c)return 0;
	if(c->entries<3)return 0;	
	
	int qual=1;
	
	//look at muon-like strip fraction
	//and sparsity (# planes with hits/#planes from front to back of muon)
	
	int muonlike=0;
	double hitplane=0.;
	
	double last_hitplane=-100.;
	
	double min = 1000000.;
	double max = 0.;
	
	double partial=(double)partialcount;
	
	for(int i=0;i<c->entries;i++)
	{
		Managed::ManagedCluster *clus = cluster_manager->GetCluster(c->cluster_id[i]);
	
		if(c->e[i]>0.5 && c->e[i]<2.5)muonlike++;
		if(fabs(last_hitplane - c->z[i])>0.03)
		{
			hitplane+=1;
			last_hitplane=c->z[i];
		}
		
		min = clus->hitplane[0]< min? clus->hitplane[0]:min;
		max = clus->hitplane[0]> max? clus->hitplane[0]:max;
	
/*		int part=IsPartiallyInstrumented(c->t[i], c->z[i], view);
		if(part)
		{
			//if(part==1) //its in z, so no question!
				partial+=1;
			//if(part==2)
			//!!! we have no knowledge of the other view in this function!
			//so we should just err on the side of caution and count all failures...
			//this will result in sparsities > 1 !
			//this will also result in under counted sparsities!
		}
*/	}
	
	double muonfrac = (double)muonlike / (double)c->entries;
	
	double nplanes = (max-min)/2.0 ; //get the total number of planes in this view that the track could go through// / 0.0708  ;
	
	//double detscale=1;
	//if(partial>0)detscale=4.*partial/hitplane + 1.*(1-partial/hitplane);
	
	
	//printf("partially instrumented planes %f max %f min %f estplanes %f hit planes %f adj planes %f\n",partial,max,min, (max-min)/0.0708, hitplane,hitplane+partial*8.0);
	
	//adjust number of hitplanes for partially instrumented regions
	hitplane += 8*partial;
	
	double sparsity = nplanes<=0? 1 : (double)hitplane / nplanes ; //length/size of u+v plane
	
	double sparsity_cut=0.8;
	
	//we don't know about the other view... so if we have some partial hits....make a large guess...count the rest as partial as well!
	if(partial/(hitplane-8*partial)>0.2 && sparsity<1)sparsity=(hitplane+8*(hitplane-partial))/nplanes;
	
	
	if (muonfrac < 0.5) qual=0;
	if (sparsity < sparsity_cut) qual=0;
	
	MSG("LongMuonFinder",Msg::kDebug)<<"muon chain check  muonfrac "<< muonfrac<<" sparsity "<< sparsity<<"  qual "<< qual<<" hitplane "<<hitplane<< " nplane "<<nplanes <<" partialhits "<<partial<<"\n";

	return qual;
}

void LongMuonFinder::MakeParticle3D()
{

	//verify that the chains in both views have sufficient overlap
	
	
	//make the 3d particle from these chains
	if(foundparticle3d){delete foundparticle3d;foundparticle3d=0;};
	foundparticle3d= new Particle3D();


	
	std::vector<point> points;
	
	double start =0;
	double end  =1000;
	
	double endu=0;
	double endv=0;
		
	
	
	
	int type=0;


		Chain *c = chain_u;
		
		if(c->particletype)
			type = c->particletype;
		
		for(int j=0;j<c->entries;j++)
		{
			if(c->parentChain==-1 && j==0)
				start = start < c->z[j]? c->z[j] : start;
			if( j == c->entries-1)
			{
				end = end < c->z[j] ? end : c->z[j];
			}
			endu=endu<c->z[j]?c->z[j]:endu;

			
			point a;
			a.z=c->z[j];
			a.t=c->t[j];
			a.e=c->e[j];
			a.chain=c->myId;
			a.chainhit=j;
			a.view=2;
			
			Managed::ManagedCluster * clu = cluster_manager->GetClusterSaver()->GetCluster(c->cluster_id[j]);
			if(clu)
			{
				a.rms_t=clu->rms_t;
				points.push_back(a);
			}
		}
	
	

		c = chain_v;
		
		if(c->particletype)
			type = c->particletype;
		for(int j=0;j<c->entries;j++)
		{
			if(c->parentChain==-1 && j==0)
				start = start < c->z[j]? c->z[j] : start;
			if( j == c->entries-1)
			{
				end = end < c->z[j] ? end : c->z[j];
			}
			endv=endv<c->z[j]?c->z[j]:endv;
			

			point a;
			a.z=c->z[j];
			a.t=c->t[j];
			a.e=c->e[j];
			a.chain=c->myId;
			a.chainhit=j;
			a.view=3;
			
			Managed::ManagedCluster * clu =cluster_manager->GetClusterSaver()->GetCluster(c->cluster_id[j]);
			if(clu)
			{
				a.rms_t=clu->rms_t;
				points.push_back(a);
			}
		}
	
	
	
	//we don't want the particle to extrapolate too far.... ...but now go all the way to the end
	double stopend = endu<endv?endv:endu;

		
	sort(points.begin(), points.end(),pointgreaterlmf);
	
	for(unsigned int i=0;i<points.size();i++)
	{
//		if( start - points[i].z > 0.045 )continue;
//		if( points[i].z - end > 0.045)continue; //within 1 planes, we want the end of the chain
	
//	printf("point %f %f %f %d\n",points[i].z,points[i].t,points[i].e,points[i].view);
	
		if(points[i].z>stopend)continue;
	
		int myview = points[i].view;
		int lower=-1;
		int upper=-1;
		for(int j=i-1;j>-1;j--)
		{
			if(points[j].view!=myview)
			{
				lower=j;
				break;
			}
		}
		for(unsigned int j=i+1;j<points.size();j++)
		{
			if(points[j].view!=myview)
			{
				upper=j;
				break;
			}
		}	
		
		
		
		double u = points[i].view==2 ? points[i].t : 0;
		double v = points[i].view==3 ? points[i].t : 0;
		double e = points[i].e;
		double z = points[i].z;
		int view = points[i].view;
		
		double rms_t = points[i].rms_t;
		

		
		if(lower>-1 && upper > -1 )// good we can extrapolate!
		{
			double s = (points[upper].t - points[lower].t) /  ( points[upper].z - points[lower].z);
			
			double t = s * (points[i].z-points[lower].z) + points[lower].t;
			
			u = myview == 2 ? u : t;
			v = myview == 3 ? v : t;
			

		}else if(lower>-1 && upper < 0)  //just user the closest other view t value
		{
		

			//do we have another lower value?
			int lower2=-1;
			for(int j=lower-1;j>-1;j--)
			{
				if(points[j].view!=myview && fabs(points[lower].z-points[j].z)>0.04)
				{
					lower2=j;
					break;
				}
			}
			
			if(lower2>-1 && fabs( points[lower].z - points[lower2].z) >0)
			{
				double s = (points[lower].t - points[lower2].t) /  ( points[lower].z - points[lower2].z);
			
				double t = s * (points[i].z-points[lower2].z) + points[lower2].t;
				u = myview == 2 ? u : t;
				v = myview == 3 ? v : t;
			
			}else{
				u = myview == 2 ? u : points[lower].t;
				v = myview == 3 ? v : points[lower].t;
			}
		}
		else if(upper>-1 && lower < 0)   //just user the closest other view t value
		{
		

		
			//do we have another upper value?
			int upper2=-1;
			for(unsigned int j=upper+1;j<points.size();j++)
			{
				if(points[j].view!=myview && fabs(points[upper].z-points[j].z)>0.04)
				{
					upper2=j;
					break;
				}
			}
			
			if(upper2>-1 && fabs( points[upper2].z - points[upper].z)>0)
			{
				double s = (points[upper2].t - points[upper].t) /  ( points[upper2].z - points[upper].z);
			
				double t = s * (points[i].z-points[upper].z) + points[upper].t;
				u = myview == 2 ? u : t;
				v = myview == 3 ? v : t;
			
			}else{
				u = myview == 2 ? u : points[upper].t;
				v = myview == 3 ? v : points[upper].t;
			}



			//lets use the vertex!

/*
			if(points[upper].view != points[i].view)
				upper=upper2;
				
			if(points[upper].view == points[i].view)
			{
*/
//dont yet have a vertex!
/*				double vz =vtx_z;
				double vt = 0;
				vt = points[upper].view == 2 ? vtx_u : vt;
				vt = points[upper].view == 3 ? vtx_v : vt;
			
				double t =vt;
				
				//if the point at z is not the same as z vtx, we need to extrapolate 
				if(fabs ( points[upper].z - vz)>0.001)
				{			
			
					double s = (points[upper].t - vt) /  ( points[upper].z - vz);
					t = s * (points[i].z-vz) + vt;			
			
				}
				
		//		printf("---view %d u %f v %f t %f\n",myview,u,v,t);
		//		printf("---vtx z %f t%f upper z %f t %f\n",vz,vt,points[upper].z,points[upper].t);
			

				u = myview == 2 ? u : t;
				v = myview == 3 ? v : t;
*/
				u = myview == 2 ? u : points[upper].t;
				v = myview == 3 ? v : points[upper].t;

//			}	
				

		}
		else if(upper==-1 && lower==-1) //we have an empty view!!!
		{

			u = myview == 2 ? u : 0;
			v = myview == 3 ? v : 0;
		}
		
		foundparticle3d->add_to_back(u,v,z,e,points[i].chain,points[i].chainhit,view,rms_t);
		
	//	printf("adding 3d point to particle %f %f %f --- %f\n",u,v,z,e);
	}
	
	
	if(type)
		foundparticle3d->particletype=(Particle3D::EParticle3DType)type;
		
	foundparticle3d->finalize();
	
	if(foundparticle3d->entries<1){delete foundparticle3d;foundparticle3d=0;}


}



Chain * LongMuonFinder::FindMuonChain(Managed::ClusterManager *cl, int view)
{

	MSG("LongMuonFinder",Msg::kDebug)<<"\n\nLooking for long muon Chain in view "<<view<<"\n";


	std::map<double, std::map<double, int> > * cluster_map = cl->GetClusterMap(view);


    std::map<double, std::map<double, int> >::reverse_iterator p_iterr;
    std::map<double, int >::iterator s_iterr;


	Chain * c = new Chain();

	double last_t=0;
	double last_d=100000;
	Managed::ManagedCluster*closest=0;
	double last_z=100000;
	double last_plane=100000;
	Managed::ManagedCluster*last_closest=0;
	int notadded=0;
	int planecount=0;
	int lastplanecount=0;
	Managed::ManagedCluster*last_added=0;
	
	Managed::ManagedCluster * closest_proj=0;
	double last_d_proj=100000;
	
	double proj_t=0;
	double lin_proj_t=0;
	
	int first=1;
	
	
/*	
	//for the first one......want to make sure we are choosing the hit closest to the track
	//go 5 planes and see where the track actually is...

	int npln=0;
	int ncnt=0;
	double mean_t=0;
    for(p_iterr=cluster_map->rbegin();p_iterr!=cluster_map->rend(); p_iterr++)
    {	
		

		if(fabs(last_plane-p_iterr->first)>0.03)npln++;
		if(npln>5)break;

         for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
         {
		 		Managed::ManagedCluster * clus = cl->GetCluster(s_iterr->second);
				ncnt++;
				mean_t += clus->t;
      	
      			if(npln==1)last_z=clus->z;
      		
      	 }
		last_plane=p_iterr->first+0.04;
		
		
	}
	
	if(ncnt>0)mean_t/=ncnt;
	
	printf("last 5 planes mean t %f\n",mean_t);

	//find the closest cluster to the mean t within 2 strips

	first=1;
	last_t=mean_t;
*/



	
	
	//for the first one......want to make sure we are choosing the hit closest to the track
	//go 5 planes and see where the track actually is...
	//look for sigle hits in each view and extrapolate to end plane

	int npln=0;
	int ncnt=0;


		double sum_z_t=0;
		double sum_z=0;
		double sum_t=0;
		double sum_z_z=0;
		
		
	double lasts_t=0;
	double lasts_z=0;
	int n=0;

	double end_z=0;

    for(p_iterr=cluster_map->rbegin();p_iterr!=cluster_map->rend(); p_iterr++)
    {	
		

		if(fabs(last_plane-p_iterr->first)>0.03)
		{
			MSG("LongMuonFinder",Msg::kDebug)<<"newplane "<<npln<<" last plane count "<<ncnt<<"\n";
			if(ncnt>0)
			{
				lasts_z/=ncnt;
				lasts_t/=ncnt;
			
				sum_z+=lasts_z;
				sum_t+=lasts_t;
				sum_z_t+=lasts_z*lasts_t;
				sum_z_z+=lasts_z*lasts_z;
				n++;
			}
		
			npln++;
			last_plane=p_iterr->first;
			ncnt=0;
			lasts_z=0;
			lasts_t=0;
		}
		if(npln>5)break;

         for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
         {
		 		Managed::ManagedCluster * clus = cl->GetCluster(s_iterr->second);
		 		if(!clus)continue;
				ncnt++;
		//	printf("cluster at z t e %f %f %f \n",clus->z,clus->t,clus->e);
      			lasts_z+=clus->z;
      			lasts_t+=clus->t;
      			if(clus->z>end_z)end_z=clus->z;
      	 }
		
		
		
	}

		
	double back_slope=0;
	double back_offset=0;
	last_t=0;

	if(n>1)
	{
		back_slope=(n*sum_z_t-sum_z*sum_t)/(n*sum_z_z-(sum_z)*(sum_z));
		back_offset=(sum_t*sum_z_z-sum_z*sum_z_t)/(n*sum_z_z-(sum_z)*(sum_z));
		last_t=back_offset + back_slope*(end_z + 0.07); //want the proj hit in the plane past what we have in the data
		last_z=last_z+0.07;
		MSG("LongMuonFinder",Msg::kDebug)<<"!!! found back slope "<<back_slope<<" offset "<<back_offset<<" endz "<<end_z<<" proj "<<last_t<<" \n";
	
	}
	
	

	MSG("LongMuonFinder",Msg::kDebug)<<"last 5 planes proj t "<<last_t<<"\n";

	//find the closest cluster to the mean t within 2 strips

	first=1;
	last_plane=1000000;
last_d=100000;


/*	
	
	for(p_iterr=cluster_map->rbegin();p_iterr!=cluster_map->rend(); p_iterr++)
    {	
    	if(fabs(last_plane-p_iterr->first)>0.03)npln++;
		if(npln>5)break;
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
        {
		 		Managed::ManagedCluster * clus = cl->GetCluster(s_iterr->second);
				if(fabs(clus->t-mean_t)<0.05)
				{
					last_d=0;
					closest=clus;
				
				
					last_d_proj=0;
					closest_proj=clus;
					last_z=clus->z;
				
					printf("in first plane check %f %f  dist %f lastdist %f  lastdistproj %f\n",clus->z, clus->t,fabs(clus->t-lin_proj_t),last_d,last_d_proj);
							
      	 			planecount++;
				
					first=0;
				}
      	}
		last_plane=p_iterr->first;
		
		if(first)
		{
			p_iterr++;
			break;
		}
	}
*/	
	
	
	
	//end check
	
			double last_notadded_z=0;
		double last_notadded_t=0;
	
	
	int punch_through=0;
	double punch_through_t=0;
	double punch_through_z=0;
	double closest_punch_through=10000;
	
	//start in the plane after where we left off from starting hit check
    for(p_iterr=cluster_map->rbegin();p_iterr!=cluster_map->rend(); p_iterr++)
    {

		 	

	//std::map<double, std::map<double, int> >::reverse_iterator checker = p_iterr;
	//checker ++;
	
	int nextplane=fabs(last_plane-p_iterr->first)>0.03;

	//if(checker !=cluster_map->rend()) if(fabs(checker->first -p_iterr->first)<0.03)nextplane=0;



    	//are we now in a different plane and do we have something to add?
    	if(nextplane && closest)
    	{

			//see if we didn't add the last guy and it is the only one in the plane and lies upon 
			//the line made by the last added point and the current point being considered
			
			MSG("LongMuonFinder",Msg::kDebug)<<"notadded "<<notadded<<", planecount "<<planecount<<", lastplanecount "<<lastplanecount<<", last_added "<<(last_added!=0)<<", last_closest "<< (last_closest!=0)<<"\n";
			if(notadded&&planecount==1&&lastplanecount==1&&last_added&&last_closest)
			{
				double slope = (last_added->t-closest->t)/(last_added->z-closest->z);
				double p  = slope * (last_closest->z-closest->z)+closest->t;
				
				MSG("LongMuonFinder",Msg::kDebug)<<"expecting not added point at "<<p<<" and it is at "<<last_closest->t<<" with z "<<last_closest->z<<"\n";
						
				if(fabs(p-last_closest->t) < 0.05 +(c->front_slope==0?0.05:0))
				{
					c->insert(last_closest->t,last_closest->z,last_closest->e,last_closest->id);
					last_closest=0;
					notadded=0;
				}
			}



    		double prev_dt=0;
    		if(c->entries>1)prev_dt=fabs(c->t[0]-c->t[1]);
//    		if(closest &&( (prev_dt<0.001?0.025:prev_dt)*2+0.05 +0.02*(fabs(last_z-closest->z)/0.035)> last_d || first))
    		if(closest &&( /*(prev_dt<0.025?0.025:prev_dt)*2*/ 0.025 + 0.01*(fabs(last_z-closest->z)/0.035) + (lastplanecount==1?0.002:0)> last_d || first))
			{
				c->insert(closest->t,closest->z,closest->e,closest->id);
				first=0;
				last_t=closest->t;
    			last_z=closest->z;
    			MSG("LongMuonFinder",Msg::kDebug)<<"adding "<<last_z<<" "<<last_t<<"\n";
    			notadded=0;
    			last_added=closest;
    			lastplanecount=planecount;
			}else{
				MSG("LongMuonFinder",Msg::kDebug)<<"failed projection to z "<<closest->z<<" t "<<closest->t<<" gap "<<last_d<<"\n";
				
				if(closest_proj)
					MSG("LongMuonFinder",Msg::kDebug)<<"checking other extrap prev_dt "<<prev_dt<<" last_z "<<last_z<<" proj_z "<<closest_proj->z<<" lpc "<<lastplanecount<<" "<<(prev_dt<0.025?0.025:prev_dt)*2 +0.2+0.002*(fabs(last_z-closest_proj->z)/0.035) + (lastplanecount==1?0.02:0)<<" > "<< last_d_proj<<"\n";
				if(closest_proj &&( /*(prev_dt<0.025?0.025:prev_dt)*2*/ 0.025 +0.01*(fabs(last_z-closest_proj->z)/0.035) + (lastplanecount==1?0.002:0)> last_d_proj || first))
				{
					c->insert(closest_proj->t,closest_proj->z,closest_proj->e,closest_proj->id);
					first=0;
					last_t=closest_proj->t;
    				last_z=closest_proj->z;
					MSG("LongMuonFinder",Msg::kDebug)<<"adding "<<last_z<<" "<<last_t<<"\n";
					notadded=0;
    				last_added=closest_proj;
    				lastplanecount=planecount;
				}else{
				
					lastplanecount=planecount;
					last_closest=closest;
					closest=0;
					closest_proj=0;
					notadded=1;
				}
			}
    		


			
    	
    		closest=0;
    		last_d=100000;
			closest_proj=0;
			last_d_proj=100000;
		}
		
		if(nextplane)
		{
			MSG("LongMuonFinder",Msg::kDebug)<<"\n\n-------\nnext plane\n";
    		MSG("LongMuonFinder",Msg::kDebug)<<"last "<<last_plane<<" cur "<<p_iterr->first<<"\n";
    	}
    
		

		if(nextplane)
		{
			planecount=0;
		
			proj_t=last_t;
			if(c->entries>3) proj_t = c->interpolate(p_iterr->first);
		 	
			lin_proj_t = c->front_slope*c->start_z+c->front_offset;
			if(lin_proj_t==0)lin_proj_t=last_t;
			MSG("LongMuonFinder",Msg::kDebug)<<"quad proj "<<proj_t<<" front proj "<<lin_proj_t<<"\n";
		 	

			//punch-through....
			//check next plane --- only  1 hit?
			//extrap to next plane
			//is the hit in the next plane sufficiently close? <1 strip?
			//set a flag....
			//then the closest in this current plane is the closest to the punchthrough line

			punch_through=0;
			punch_through_t=0;
			punch_through_z=0;
			closest_punch_through=10000;	
			std::map<double, std::map<double, int> >::reverse_iterator z_it = p_iterr;
			z_it++;
			std::map<double, int>::iterator t_it;
			
			//max number of planes to look forward
			int maxcheck=4;
			
//			printf("next plane\n");
			while(maxcheck>0 && z_it !=cluster_map->rend())
			{

			
			punch_through=0;
			punch_through_t=0;
			punch_through_z=0;
			double last_z = z_it->first;
			int cnt=0;
		//	printf("\n\nlast_z %f\n",last_z);
			for(;z_it!=cluster_map->rend();z_it++)
			{
				if(fabs(last_z - z_it->first)>0.03)break;
				for(t_it=z_it->second.begin();t_it!=z_it->second.end();t_it++)
				{
					cnt++;
					if(cnt>1)break;
					punch_through_z=z_it->first;
					punch_through_t=t_it->first;
		//			printf("%f %f   \n",punch_through_z,punch_through_t);
				}
				if(cnt>1)break;
			}
		//	printf("cnt %d\n",cnt);
			if(cnt==1)
			{
				//is that point close to where we expect it?
				double exp  = c->front_slope*(punch_through_z) +c->front_offset ;
				
		//		printf("fs %f fo %f cz %f ct %f pz %f pt %f\n",c->front_slope,c->front_offset,c->start_z,c->start_t,punch_through_z,punch_through_t);
				
		//		printf("ext %f  diff %f\n",exp,fabs(exp-punch_through_t));
				if(fabs(exp-punch_through_t)<0.02)
				{
					punch_through=1;
					//save the expected point in the plane currently under consideration!
					double curz=p_iterr->first;
					double dz = (c->start_z - punch_through_z);
					double dt = (c->start_t - punch_through_t);
					
			//		printf("pt %f pz %f dz %f dt %f sl %f\n",punch_through_t,punch_through_z,dz,dt,dt/dz);
					

					punch_through_t=dz ? dt/dz*(curz-punch_through_z)+punch_through_t : punch_through_t;
					punch_through_z=curz;
				}
			}

		//	if(punch_through)printf("punch through found at z %f t%f\n",punch_through_z,punch_through_t);
			if(punch_through)break;
			
			z_it++;
			maxcheck--;			
			}

		}	
		

         for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
         {
		 		Managed::ManagedCluster * clus = cl->GetCluster(s_iterr->second);


				if(fabs(clus->t-lin_proj_t)<last_d)
				{
					last_d=fabs(clus->t-lin_proj_t);
					closest=clus;
				}
				
				if(fabs(clus->t-proj_t)<last_d_proj)
				{
					last_d_proj=fabs(clus->t-proj_t);
					closest_proj=clus;
				}
				
				//also consider a hit at the same dt
				if(last_d < 0.05 && fabs(clus->t-last_t)<last_d)
				{
					last_d=fabs(clus->t-last_t);
					closest=clus;
				}	
				
				//also consider if the last hit was not added... if we added that last signle hit, would this hit then fit better?
				if(notadded && lastplanecount==1)
				{
					double slope = (last_notadded_t-last_t)/(last_notadded_z-last_z);
					double proj = last_t - (last_z-clus->z) * slope;
					MSG("LongMuonFinder",Msg::kDebug)<<"notadded fit projects to z "<<clus->z<<" t "<<proj<<"\n";
					MSG("LongMuonFinder",Msg::kDebug)<<"slope from past points z t "<<last_notadded_z<<" "<<last_notadded_t<<"    "<<last_z<<" "<<last_t<<"\n";
					if(fabs(proj - clus->t)<0.15)
					{
						closest=clus;
						last_d =  fabs(proj - clus->t);
						MSG("LongMuonFinder",Msg::kDebug)<<"closest!\n";
					} 
				}
				
				//give priority to the punchthrough
				if(punch_through)
				{
					if(fabs(clus->t-punch_through_t)<closest_punch_through)
					{
						closest_punch_through=fabs(clus->t-punch_through_t);
						closest=clus;
					}	
				}
				
				MSG("LongMuonFinder",Msg::kDebug)<<clus->z<<" "<< clus->t<< " dist "<<fabs(clus->t-lin_proj_t)<<" lastdist "<<last_d<<" lastdistproj "<<last_d_proj<<"\n";
							
      	 	planecount++;
      	 }
		last_plane=p_iterr->first;
		
		//pick up the last hit details... only use if you know there is only one hit in this plane!
		for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
        {
		 	Managed::ManagedCluster * clus = cl->GetCluster(s_iterr->second);
			last_notadded_z=clus->z;
			last_notadded_t=clus->t;
		}
	}


    		double prev_dt=0;
    		if(c->entries>1)prev_dt=fabs(c->t[0]-c->t[1]);
//    		if(closest &&( (prev_dt<0.001?0.025:prev_dt)*2+0.05 +0.02*(fabs(last_z-closest->z)/0.035)> last_d || first))
    		if(closest &&( /*(prev_dt<0.025?0.025:prev_dt)*2*/ 0.025 + 0.01*(fabs(last_z-closest->z)/0.035) + (lastplanecount==1?0.002:0)> last_d || first))
			{
				c->insert(closest->t,closest->z,closest->e,closest->id);
				first=0;
				last_t=closest->t;
    			last_z=closest->z;
    			MSG("LongMuonFinder",Msg::kDebug)<<"adding "<<last_z<<" "<<last_t<<"\n";
    			notadded=0;
    			last_added=closest;
    			lastplanecount=planecount;
			}else{
				if(closest)MSG("LongMuonFinder",Msg::kDebug)<<"failed projection to z "<<closest->z<<" t "<<closest->t<<" gap is "<<last_d<<"\n";
				else MSG("LongMuonFinder",Msg::kDebug)<<"no closest\n";
				
				if(closest_proj)
					MSG("LongMuonFinder",Msg::kDebug)<<"checking other extrap prev_dt "<<prev_dt<<" last_z "<<last_z<<" proj_z "<<closest_proj->z<<" lpc "<<lastplanecount<<" "<<(prev_dt<0.025?0.025:prev_dt)*2 +0.2+0.002*(fabs(last_z-closest_proj->z)/0.035) + (lastplanecount==1?0.02:0)<<" > "<< last_d_proj<<"\n";
				if(closest_proj &&( /*(prev_dt<0.025?0.025:prev_dt)*2*/ 0.025 +0.01*(fabs(last_z-closest_proj->z)/0.035) + (lastplanecount==1?0.002:0)> last_d_proj || first))
				{
					c->insert(closest_proj->t,closest_proj->z,closest_proj->e,closest_proj->id);
					first=0;
					last_t=closest_proj->t;
    				last_z=closest_proj->z;
					MSG("LongMuonFinder",Msg::kDebug)<<"adding "<<last_z<<" "<<last_t<<"\n";
					notadded=0;
    				last_added=closest_proj;
    				lastplanecount=planecount;
				}else{
				
					lastplanecount=planecount;
					last_closest=closest;
					closest=0;
					closest_proj=0;
					notadded=1;
				}
			}
    		


			
    	
   
			
			
return c;			


}


void LongMuonFinder::DumpParticle3D()
{
	ostringstream s;
	
	s << "Long Muon Particle3D found "<<endl;

		Particle3D * p = foundparticle3d;
		if (p==0)return;
	
	

	
		s << "\n---Particle3d "  << "---\nstart, stop (u,v,z) (" << p->start_u << ", "<< p->start_v << ", " << p->start_z <<") (" <<p->end_u<<", "<< p->end_v << ", " << p->end_z <<")"<<endl;
		s << "entries "<< p->entries << " muonlike " << p->muonfrac <<endl;
		

		
		s<<"types: ";
		for(unsigned int j=0;j<p->types.size();j++)
		{
			switch( p->types[j].type)
			{
				case 	ParticleType::em:
					s<<"em ";
					break;
				case ParticleType::emshort:
					s<<"emshort  ";
					break;				
				case ParticleType::muon:
					s<<"muon ";
					break;
				case ParticleType::prot:
					s<<"prot ";
					break;		
				case ParticleType::pi0:
					s<<"pi0 ";
					break;
				case ParticleType::uniquemuon:
					s<<"uniquemuon ";
					break;	
			
			}
				
		}
		s<<endl;
		
		s<<"particletype : "<<p->particletype<<endl;
		
		
		s<<"par a " << p->emfit_a << " par b "<< p->emfit_b << " par e0 "<< p->emfit_e0 << " cale "<<p->calibrated_energy<<" chisq " << p->emfit_chisq<<" ndf " << p->emfit_ndf<<endl;
		s<<"emprob " << p->emfit_prob <<" avg rms_t "<<p->avg_rms_t<<endl;
				
		s << "points (u,v,z,e - chain, chainhit, chainview - shared - rms_t - view) : ";
		for(int j=0;j<p->entries;j++)
			s << "(" << p->u[j]<<", "<<p->v[j]<<", "<<p->z[j]<<", "<<p->e[j]<<" - " << p->chain[j] <<", "<<p->chainhit[j]<<", "<<p->view[j]<<" - "<<p->shared[j]<<" - "<< p->rms_t[j]<< " - " << p->view[j]<<") ";
	
	/*
	s << "points (e - shared) : ";
		for(int j=0;j<p->entries;j++)
			s << "(" <<p->e[j]<<" - "<<p->shared[j]<<") ";
	
	*/
	
	

			
			

			
			
		s<<endl<<endl;
		
		
	
	
	

	MSG("LongMuonFinder",Msg::kDebug) <<s.str();

	
}

