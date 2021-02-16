#include "NueAna/ParticlePID/ParticleFinder/PrimaryShowerFinder.h"
#include <map>
#include <math.h>
#include <algorithm>

#include "MessageService/MsgService.h"

#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedCluster.h"
#include <sstream>
#include "NueAna/ParticlePID/ParticleFinder/ShwFit.h"


ClassImp(PrimaryShowerFinder)
CVSID("$Id: PrimaryShowerFinder.cxx,v 1.3 2009/09/11 05:00:40 gmieg Exp $");


TH2D * PrimaryShowerFinder::houghmapU =0;
TH2D * PrimaryShowerFinder::houghmapV =0;

TH2D * PrimaryShowerFinder::intU=0;
TH2D * PrimaryShowerFinder::intV=0;	

ChainHelper *PrimaryShowerFinder::chu=0;
ChainHelper *PrimaryShowerFinder::chv=0;

Particle3D * PrimaryShowerFinder::foundparticle3d=0;

Chain * PrimaryShowerFinder::chain_u=0;
Chain * PrimaryShowerFinder::chain_v=0;
		

struct spoint
{
	int view;
	double z;
	double t;
	double e;
	int chain;
	int chainhit;
	double rms_t;
};
	
bool spointgreaterlmf(spoint p1, spoint p2){return p1.z < p2.z;}

PrimaryShowerFinder::PrimaryShowerFinder()
{

	Reset();
}

PrimaryShowerFinder::~PrimaryShowerFinder()
{
	Reset();
}

void PrimaryShowerFinder::Reset(int reset_vertex)
{
	ran=0;
	if(foundparticle3d){delete foundparticle3d;foundparticle3d=0;}
	foundparticle=0;
	chain_u=0; 
	chain_v=0;
	single_view_long_shower=0;


//	if(houghmapU)delete houghmapU;
//		houghmapU=0;
//	if(houghmapV)delete houghmapV;
//		houghmapV=0;


	if(houghmapU)houghmapU->Reset();
	if(houghmapV)houghmapV->Reset();

	houghlinesU.clear();
	houghlinesV.clear();
	houghlineMatch.clear();

	if(intU)delete intU;
		intU=0;
	if(intV)delete intV;
		intV=0;

	if(reset_vertex)
	{
		vtx_u=0;
		vtx_v=0;
		vtx_z=0;
		foundvertex=0;
		foundvertexU=0;
		foundvertexV=0;
	}
}

void PrimaryShowerFinder::MakeChains(int view)
{

	ChainHelper *ch = view==2?chu:chv;
	if(!ch)return;
	

	

	
	std::vector<HoughLine> * hl = view==2?&houghlinesU:&houghlinesV;
	
	if(hl->size()<1)
	{
		//does the view have only one cluster? if so... make that a chain!
		std::map<double, std::map<double, int> > * cluster_map = cluster_manager->GetClusterMap(view);
  		std::map<double, std::map<double, int> >::reverse_iterator p_iterr;
    	std::map<double, int >::iterator s_iterr;

		Managed::ManagedCluster *cluster=0;
		int cluster_count=0;
	 	for(p_iterr=cluster_map->rbegin();p_iterr!=cluster_map->rend(); p_iterr++)
    	{	
 	
        	for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
			{
				cluster = cluster_manager->GetCluster(s_iterr->second);
	   			cluster_count++;
	   		}
    	}
    	if(cluster_count==1 && cluster)
    	{
			int nc = ch->NewChain();
			Chain *c = ch->GetChain(nc);
			
		
		//	int newid = cluster_manager->SaveCluster(cluster->id,cluster->e,2);
		//	cluster = cluster_manager->GetCluster(newid);
			
			if(cluster)
			{
				c->insert(cluster->t, cluster->z, cluster->e, cluster->id);
			}    	
			MSG("PrimaryShowerFinder",Msg::kDebug)<<"single cluster chain made in view "<<view<<"\n";
    	
    	}
		
	
		return;
	}
	
	
	//make a copy....
	
	std::vector<HoughLine> hlcopy = *hl;

	hl=&hlcopy;

	Managed::ClusterManager *mycm=cluster_manager;
	

	//probably vertex 
	//to check if chains after the first are vertex pointing or other chain pointing
	double prob_vtx_z=0;
	double prob_vtx_t=0;

	if(view ==2 && foundvertexU || foundvertex)
	{
		prob_vtx_z=vtx_z;
		prob_vtx_t=vtx_u;
		
		//printf("USING VERTEX U %f %f\n",prob_vtx_z,prob_vtx_t);
	}

	if(view ==3 && foundvertexV || foundvertex)
	{
		prob_vtx_z=vtx_z;
		prob_vtx_t=vtx_u;
		
		//printf("USING VERTEX V %f %f\n",prob_vtx_z,prob_vtx_t);
	}	

/////////
	//uncomment to do a temp	
	Managed::ClusterManager cm=*cluster_manager;
	cm.SetClusterSaver(0);
	cm.SetHitManager(0);
	
	mycm=&cm;
	mycm->ClearInUse();

////////

	
	std::vector<int>foundChains;//store the found chains here as well... so we can split up hough lines as needed
	
	
	//need to load any existing chains from the long muon finder!!!!
	/////////
	for(unsigned int i=0;i<ch->finished.size();i++)
	{
		foundChains.push_back(ch->finished[i].myId);
	
	}
	
	
	
	/////////
	

	//do a full reload...
	for(unsigned int k=0;k<hl->size();k++)
	{
//		printf("%d\n",k);
		(*hl)[k].ResetHits(0);
		LoadCloseHits(&(*hl)[k],mycm,view,0.0412*1.5);
//		printf("%f %f   %f %f\n",(*hl)[k].start_z,(*hl)[k].start_t,(*hl)[k].end_z,(*hl)[k].end_t);	
	}
	std::vector<HoughLine> hlcopy2 = *hl;
		
//	
//	printf("!!!!!!!!!!!!!!!!!!!!!!!\nlooking for chains\n");
	for(int ccc=0;ccc<10;ccc++)
	{
	//	printf("try %d with %d houghlines\n",ccc,hl->size());
		for(unsigned int k=0;k<hl->size();k++)
		{
	
	//		printf("%d\n",k);
			(*hl)[k].ResetHits(1); //keep the bounds
			
		//	printf("reset hl %d\n",k);
			LoadCloseHits(&(*hl)[k],mycm,view,0.0412*1.5,1);
	//		printf("current %d,%f %f %f %f\n",(&(*hl)[k])->ncluster,(&(*hl)[k])->start_z,(&(*hl)[k])->start_t,(&(*hl)[k])->end_z,(&(*hl)[k])->end_t);
			
			
//			printf("loaded %d hits\n",(*hl)[k].ncluster);
			
		
			//does this hough line intersect any found chains?  if so, we should split it!
			for(unsigned int c=0;c<foundChains.size();c++)
			{
//				printf("c %d of %d\n",c,foundChains.size());
//				printf("current %d,%f %f %f %f\n",(&(*hl)[k])->ncluster,(&(*hl)[k])->start_z,(&(*hl)[k])->start_t,(&(*hl)[k])->end_z,(&(*hl)[k])->end_t);
				
				Chain *chain = ch->GetChain(foundChains[c]);
				if(!chain)continue;
				if(chain->entries<2)continue;
				
				HoughLine *lnew = SplitHoughLine(chain,&(*hl)[k],mycm);
//				printf("split  %d,%f %f %f %f\n",(&(*hl)[k])->ncluster,(&(*hl)[k])->start_z,(&(*hl)[k])->start_t,(&(*hl)[k])->end_z,(&(*hl)[k])->end_t);
				
							
				if(lnew && lnew->ncluster>0)
				{
///					printf("new from split\n");
					hl->push_back(*lnew);
//					printf("%d %f %f %f %f\n",lnew->ncluster,lnew->start_z,lnew->start_t,lnew->end_z,lnew->end_t);
				}
				
				if(lnew)delete lnew;
				//delete it if it is empty!
				if((*hl)[k].ncluster<1)
				{
//					printf("hl size %d hl at position %d is empty...\n",hl->size(),k);
					
				//	std::vector<HoughLine> hh;
				//	for(unsigned int aw=0;aw<hl->size();aw++)
				//		if(aw!=k)hh.push_back((*hl)[aw]);
					
				//	*hl=hh;
					
					hl->erase(hl->begin()+k,hl->begin()+k+1);
					k--;
//					printf("%d %d\n",k,hl->size());
					break;//restart the check now that this hl is split
				}

			}
//			printf("e %d\n",k);
		}
	
		CleanHoughLines(0,hl);
		
		//so we know where to attach secondary chains
		std::vector<int>closest_chain;
		std::vector<double>closest_chain_z;
		for(unsigned int k=0;k<hl->size();k++)
		{
			closest_chain.push_back(-1);
			closest_chain_z.push_back(0);
		}


		std::vector<int>toskip;

		if( ( view ==2 && foundvertexU == 1 ) || ( view==3 && foundvertexV==1)  || foundvertex) //if we don't have a vertex, we don't have any chains... so consider all of them and skip this part....
		for(unsigned int k=0;k<hl->size();k++)
		{
			HoughLine * l = &(*hl)[k];

			int keep=0;
			//check for chain pointing	
			double lastdist=10000;
			for(unsigned int i=0;i<foundChains.size();i++)
			{	
				Chain *c = ch->GetChain(foundChains[i]);
				if(!c)continue;
				if(c->entries<2)continue; //can come from when a chain is split at the first cluster
		


				//measure the distance of intersection
				double vz = vtx_z;
				double vt = view==2? vtx_u:vtx_v;
				
				
				double c_start_z, c_start_t, c_end_z, c_end_t =0;
				
				if( fabs((c->start_z-vz)*(c->start_z-vz)+(c->start_t-vt)*(c->start_t-vt)) < fabs((c->end_z-vz)*(c->end_z-vz)+(c->end_t-vt)*(c->end_t-vt)) )
				{
					c_start_z = c->start_z;
					c_start_t = c->start_t;
					c_end_z = c->end_z;
					c_end_t = c->end_t;	
				}else{
					c_start_z = c->end_z;
					c_start_t = c->end_t;
					c_end_z = c->start_z;
					c_end_t = c->start_t;				
				}
				
			         if(!(c_end_z-c_start_z))continue;

				double l_start_z, l_start_t, l_end_z, l_end_t =0;
				
				if( fabs((l->start_z-c_start_z)*(l->start_z-c_start_z)+(l->start_t-c_start_t)*(l->start_t-c_start_t)) < fabs((l->end_z-c_start_z)*(l->end_z-c_start_z)+(l->end_t-c_start_t)*(l->end_t-c_start_t)) )
				{
					l_start_z = l->start_z;
					l_start_t = l->start_t;
					l_end_z = l->end_z;
					l_end_t = l->end_t;	
				}else{
					l_start_z = l->end_z;
					l_start_t = l->end_t;
					l_end_z = l->start_z;
					l_end_t = l->start_t;				
				}
				

				//find intersection
				
				double lslope = (l_end_t-l_start_t) / (l_end_z-l_start_z) ;
				double loffset = l_end_t - l_end_z*lslope;
				
				double cslope = (c_end_t-c_start_t) / (c_end_z-c_start_z) ;
				double coffset = c_end_t - c_end_z*cslope;
				
				double int_z = cslope-lslope ? (loffset - coffset ) / ( cslope - lslope):l_start_z;
				double int_t = cslope * int_z + coffset;
												
				//require intersection to be behind the start of the chain

				//if(fabs((c_start_z-vz)*(c_start_z-vz)+(c_start_t-vt)*(c_start_t-vt)) > fabs((int_z-vz)*(int_z-vz)+(int_t-vt)*(int_t-vt)) )continue;
				
				//calculate distance along chain path... where start is 0 and going towards end increases
				//if its <0  its bad....
				double r_c_end_z=c_end_z-c_start_z;
				double r_c_end_t=c_end_t-c_start_t;
				double r_l_end_z=l_end_z-c_start_z;
				double r_l_end_t=l_end_t-c_start_t;
				double r_l_start_z=l_start_z-c_start_z;
				double r_l_start_t=l_start_t-c_start_t;
				
				double theta = atan(r_c_end_t/r_c_end_z);				
				double rr_int_z = (int_z-c_start_z) * cos(theta) + (int_t-c_start_t) * sin(theta);				
				double rr_int_t = (int_t-c_start_t) * cos(theta) - (int_z-c_start_z) * sin(theta);					
				if(rr_int_z<0.03)continue;// points to far fowards of chain
				
				
				
				
				double ltheta = atan( ( l->end_t-l->start_t) / (l->end_z-l->start_z) );				
				double rr_l_int_z = (int_z-l->start_z) * cos(ltheta) + (int_t-l->start_t) * sin(ltheta);				
	
				//printf("in hl coords, int is at z %f\n",rr_l_int_z);
	
				//require the intersection to be forwards of the hl!
				if(rr_l_int_z>0)continue;				
				
						
				
				//its the intersection past the end of the chain?  if so... check the relative angle between the hl and the chain to keep it reasonable
				double rr_c_end_z = r_c_end_z * cos(theta) + r_c_end_t * sin(theta);					
				//double rr_c_end_t =  -r_c_end_z * sin(theta) + r_c_end_t * cos(theta); 				


				//rotate the coordinates so the chain is along the x axis

				double rr_l_end_z = r_l_end_z * cos(theta) + r_l_end_t * sin(theta);
				double rr_l_end_t =  -r_l_end_z * sin(theta) + r_l_end_t * cos(theta); 
				double rr_l_start_z = r_l_start_z * cos(theta) + r_l_start_t * sin(theta);
				double rr_l_start_t =  -r_l_start_z * sin(theta) + r_l_start_t * cos(theta);				
	
				
				if(rr_int_z > rr_c_end_z)
				{
					//count the angle from the intersection!
					double dy = rr_l_start_t - rr_int_t;
					double dx = rr_l_start_z - rr_int_z;
					
					double intangle = dx ? atan(dy/dx) : 0;
					if(intangle>3.141592/6 || dx==0)continue; //to steep!
				
				}
				
				//is it the closest?  
				
				double dist = sqrt((l_start_z - int_z)*(l_start_z - int_z)+(l_start_t - int_t)*(l_start_t - int_t));
				
				//if the intersection is past the end of the chain... need to include the distance from the intersection to the end of the chain
				if(fabs((c_end_z-vz)*(c_end_z-vz)+(c_end_t-vt)*(c_end_t-vt)) < fabs((int_z-vz)*(int_z-vz)+(int_t-vt)*(int_t-vt)) )
				{
					dist+= sqrt((c_end_z - int_z)*(c_end_z - int_z)+(c_end_t - int_t)*(c_end_t - int_t));
				}

				//is the closest ?
				if(dist > lastdist)continue;
				//printf("closer to chain %f %f %f %f at z %f\n",c_start_z,c_start_t,c_end_z,c_end_t,int_z);
				//if so, check for orientation before recording
			
				if(fabs(rr_l_end_t) <= fabs(rr_l_start_t) && rr_l_start_z < rr_l_end_z ||
					fabs(rr_l_end_t) >= fabs(rr_l_start_t) && rr_l_start_z > rr_l_end_z 
				 )
				{
					//its closer... but no good... we want to invalidate the last guess...
					closest_chain[k]=-1;
					closest_chain_z[k]=int_z;
					lastdist=dist;	
					keep=0;				
					continue;
				}

				//printf("in chain coords line has start %f %f end %f %f\n",rr_l_start_z,rr_l_start_t,rr_l_end_z,rr_l_end_t);
				//printf("hl start %f %f end %f %f\n",l_start_z,l_start_t,l_end_z,l_end_t);
				//record it
				closest_chain[k]=foundChains[i];
				closest_chain_z[k]=int_z;
				lastdist=dist;
				keep=1;
				
				//printf("keeping %d it points to chain %f %f %f %f at z %f\n",k,c_start_z,c_start_t,c_end_z,c_end_t,int_z);
				
			}
			
			
			//dont have anything yet?
			//see if we are vertex pointing
			if(closest_chain[k]<0)
			{
				double exp_vtx_t = l->GetExpectedT(vtx_z);
				if((view==2&&foundvertexU || view==3&&foundvertexV) && fabs(exp_vtx_t-(view==2?vtx_u:vtx_v))<2*0.0451) //smaller than 2 strips away?
				{
				
					//need to make sure that the vertex is not within the hl!
	
				
					double theta = atan( ( l->end_t-l->start_t) / (l->end_z-l->start_z) );				
					double rr_end_z = (l->end_z-l->start_z) * cos(theta) + (l->end_t-l->start_t) * sin(theta);				
					//double rr_end_t = (l->end_t-l->start_t) * cos(theta) - (l->end_z-l->start_z) * sin(theta);					

					double rr_vtx_z = (vtx_z-l->start_z) * cos(theta) + (exp_vtx_t-l->start_t) * sin(theta);				
					//double rr_vtx_t = (exp_vtx_t-l->start_t) * cos(theta) - (vtx_z-l->start_z) * sin(theta);	
	
					//printf("in hl coords:  hl 0 0 to %f %f  vtx %f %f\n",rr_end_z,rr_end_t,rr_vtx_z,rr_vtx_t);
	
					if(rr_vtx_z< 0 || rr_vtx_z > rr_end_z)
					{				
						//printf("keeping %d vp estimate vtx %f %f actual %f %f\n",k,vtx_z,l->GetExpectedT(vtx_z),vtx_z,view==2?vtx_u:vtx_v);
						keep=1; //we'll keep it
					}
				}	
			}
		

			if(keep)
			{
				//printf("keeping %d\n",k);
				continue;
			}
			//if we get here... we don't want this hl!
			//printf("skipping %d\n",k);
			toskip.push_back(k);
		
		}


/*


		std::vector<int>toskip;
		//make sure that all hough line are vertex pointing or pointing to a another chain and on the correct angle
		if(foundChains.size()>0)
		{
			for(unsigned int k=0;k<hl->size();k++)
			{
				HoughLine * l = &(*hl)[k];
	//					printf("current %d,%f %f %f %f\n",(&(*hl)[k])->ncluster,(&(*hl)[k])->start_z,(&(*hl)[k])->start_t,(&(*hl)[k])->end_z,(&(*hl)[k])->end_t);
						
							
				int keep=0;				
				//is it vertex pointing?
				if(foundvertex && fabs(l->GetExpectedT(prob_vtx_z)-prob_vtx_t)<2*0.0451) //smaller than 2 strips away?
				{
					printf("keeping %d vp\n",k);
					keep=1; //we'll keep it
					//we still want to see if it is pointing to another chain so we know if it needs to be a secondary
				}

				//now we need to see if it is pointing to a current chain ... and find the closest one!
				int closest_pointing=-1;
				double closest_dist=0;
				double dist=0;		
				for(unsigned int i=0;i<foundChains.size();i++)
				{
					Chain *c = ch->GetChain(foundChains[i]);
					double exp_t_start = l->GetExpectedT(c->start_z);
					double exp_t_end = l->GetExpectedT(c->end_z);

			
					if(( exp_t_start < c->start_t && exp_t_end > c->end_t ) || ( exp_t_start > c->start_t && exp_t_end < c->end_t))
					{
					//		//check for orientation
					//		
					//		if(fabs(c->end_z - l->end_z) > fabs(c->end_z - l->start_z))continue;
					//		if(fabs(c->start_z - l->start_z) > fabs(c->start_z - l->end_z))continue;
							
							MSG("PrimaryShowerFinder",Msg::kDebug)<<"Passing orientation on "<<c->start_z<<" "<<c->start_t<<" "<<c->end_z<<" "<<c->end_t<<"\n";

							//check the distance
							//the way we found the intersections requires the intersection to be within the chain....
							for(unsigned int kk=0;kk<c->z.size()-1;kk++)
							{
								double exp_t = l->GetExpectedT(c->z[kk]);
								double exp_t_n = l->GetExpectedT(c->z[kk+1]);
								
								if(( exp_t>=c->t[kk]&&exp_t_n<=c->t[kk+1] ) || ( exp_t<=c->t[kk]&&exp_t_n>=c->t[kk+1] ) )
								{//intersection is bewteen these two points!
							
										//for now, put it close to the more forward hit...
										//later we can actually extrapolate it
										
										double t = c->t[kk];
										dist = fabs(t-exp_t);
										closest_chain_z[k]=c->z[kk];
									
								}
							
							}
							
							
						if(closest_pointing<0)//the first one
						{
							closest_pointing=foundChains[i];
							closest_dist = 	dist;		
						}else
						{
							if(dist<closest_dist)
							{
								closest_pointing=foundChains[i];
								closest_dist = 	dist;								
							}
						}
				
						


					}
				}
				
				for(unsigned int i=0;i<foundChains.size();i++)
				{		
									Chain *c = ch->GetChain(foundChains[i]);	
				
	//printf("end chain %f %f  start hl %f %f\n",c->end_z,c->end_t,l->start_z,l->start_t);

				
				
				if(c->end_z < l->start_z)
				{
				
							//also check to see if this chain is at the end of another chain
							if(
								fabs(c->end_t - l->GetExpectedT(c->end_z))<0.1  
								&& l->start_z - c->end_z  <0.07*3)
							{
							
								closest_chain_z[k]=c->end_z;
								dist = fabs(c->end_t - l->GetExpectedT(c->end_z));

							}

						//	printf("end chain %f %f  start hl %f %f\n",c->end_z,c->end_t,l->start_z,l->start_t);
				
				}
				
				
				
						if(closest_pointing<0)//the first one
						{
							closest_pointing=foundChains[i];
							closest_dist = 	dist;		
						}else
						{
							if(dist<closest_dist)
							{
								closest_pointing=foundChains[i];
								closest_dist = 	dist;								
							}
						}
				
				
				
				}
				
				
				if(closest_pointing>-1)
				{
						Chain *c = ch->GetChain(closest_pointing);
						//its pointing... so now check the orientation
						double close_z;
						double close_t;
						printf("prob vtx z %f t %f\n",prob_vtx_z,prob_vtx_t);
						printf("chain start z %f t %f, end z %f t %f\n",c->start_z,c->start_t,c->end_z,c->end_t);
						if(fabs((c->start_t-prob_vtx_t)*(c->start_t-prob_vtx_t)+(c->start_z-prob_vtx_z)*(c->start_z-prob_vtx_z))
							<
							fabs((c->end_t-prob_vtx_t)*(c->end_t-prob_vtx_t)+(c->end_z-prob_vtx_z)*(c->end_z-prob_vtx_z)))
						{	
							close_z=c->start_z;
							close_t=c->start_t;
						}
						else
						{
							close_z=c->end_z;
							close_t=c->end_t;
						}
						
						printf("close z %f t %f\n",close_z,close_t);
						double ch_close_z;
						double ch_close_t;
						double ch_far_z;
						double ch_far_t;
						
			//			if(fabs(l->start_z - close_z) < fabs(l->end_z - close_z))
			//			{
						if(fabs((l->start_t-prob_vtx_t)*(l->start_t-prob_vtx_t)+(l->start_z-prob_vtx_z)*(l->start_z-prob_vtx_z))
							<
							fabs((l->end_t-prob_vtx_t)*(l->end_t-prob_vtx_t)+(l->end_z-prob_vtx_z)*(l->end_z-prob_vtx_z)))
						{	
			
			
							ch_close_z = l->start_z;
							ch_close_t = l->start_t;
							ch_far_z = l->end_z;
							ch_far_t = l->end_t;
						}else{
							ch_close_z = l->end_z;
							ch_close_t = l->end_t;
							ch_far_z = l->start_z;
							ch_far_t = l->start_t;
						}
						
						
						//this will cause issues with vertical chains!
						double dst = fabs(c->interpolate(ch_close_z)-ch_close_t);
						double det = fabs(c->interpolate(ch_far_z)-ch_far_t);
						printf("!!! %f\n",c->interpolate(ch_far_z));
						printf(" close z %f t %f  far z %f t %f  dst %f det %f\n",ch_close_z,ch_close_t,ch_far_z,ch_far_t,dst,det);
						
						if(dst<det)
						{
							keep=1;//we'll keep it
							printf("keeping %d cp\n",k);
							closest_chain[k]=closest_pointing;
						}
							
				
				}
				
				if(keep)continue;
				//if we get here... we don't want this hl!
				printf("skipping %d\n",k);
				toskip.push_back(k);
	
			}
		}
	


*/







	for(unsigned int i=0;i<hl->size();i++)
	{
		MSG("PrimaryShowerFinder",Msg::kDebug)<<"hl "<<i<<" "<<(*hl)[i].start_z<<" "<<(*hl)[i].start_t<<" "<<(*hl)[i].end_z<<" "<<(*hl)[i].end_t<<"\n";
	}
	
	
	
	std::map<double,int> orderer;
//	printf("\n");
	for(int i=0;i<(int)hl->size();i++)
	{

		int save=1;
		for(int j=0;j<(int)toskip.size();j++)
		{
			if(toskip[j]==i)
			{
				save =0;
				break;
			}
		}
		if(!save)continue;

		HoughLine * l = &(*hl)[i];
	
		double weight = l->ncluster *cos(l->phi)*l->sum_e*l->sum_e;

		orderer.insert(make_pair(weight,i));
		
//		printf("%d %f %f | %f\n",l->ncluster,l->phi,l->sum_e,weight);
	}		
//	printf("\n");

	std::map<double,int>::reverse_iterator it_orderer;

	HoughLine *h = 0;	

	//int cnt=0;
	//for(it_orderer = orderer.rbegin();it_orderer!=orderer.rend();it_orderer++)
	//{
	
	it_orderer = orderer.rbegin();
	while(it_orderer!=orderer.rend())
	{
		if((*hl)[it_orderer->second].ncluster<2)it_orderer++;
		else break;
	}
	
	if(it_orderer==orderer.rend())break;
	
		int i= it_orderer->second;
		
		HoughLine * l = &(*hl)[i];
		h=l;
		
		MSG("PrimaryShowerFinder",Msg::kDebug)<<"chose "<<l->ncluster<<" "<<l->phi<<" "<<l->sum_e<<" | "<<it_orderer->first<<"\n";

	//	cnt++;
	//	if(cnt>3)break;


		int nc = ch->NewChain();
		Chain *c = ch->GetChain(nc);
		foundChains.push_back(nc);

		MSG("PrimaryShowerFinder",Msg::kDebug)<<"making chain from houghline "<<i<<" "<<l->start_z<<" "<<l->start_t<<" "<<l->end_z<<" "<<l->end_t<<"\n";
		MSG("PrimaryShowerFinder",Msg::kDebug)<<"STARTING WITH "<< ch->finished.size()<<" CHAINS\n";
	
		for(unsigned int i=0;i<h->cluster_id.size();i++)
		{
			Managed::ManagedCluster *cl = mycm->GetCluster(h->cluster_id[i]);
			
			
		//	int newid = mycm->SaveCluster(cl->id,cl->e,2);
		//	cl = mycm->GetCluster(newid);
			
			if(cl)
			{
			
				
				
				
				c->insert(cl->t, cl->z, cl->e, cl->id);
				MSG("PrimaryShowerFinder",Msg::kDebug)<<"("<<cl->t<<" "<< cl->z<<" "<< cl->e<<" "<< cl->id<<") s"<<cl->GetStatus()<<" v"<<cl->view<<"\n";
			}

			mycm->MarkUsed(cl->id);
		
		}
		//attach the chain if needed...
		
		if(closest_chain[i]>-1)
		{

			Chain *mom = ch->GetChain(closest_chain[i]);
			Chain newdaughter = ch->AttachAt(mom,c,closest_chain_z[i]);
			if(newdaughter.entries>0)ch->AddFinishedChain(newdaughter);
			MSG("PrimaryShowerFinder",Msg::kDebug)<<"kept by closest_chain, attempting to attach to chain "<<mom->start_z<<" "<<mom->start_t<<" "<<mom->end_z<<" "<<mom->end_t<<" at z "<<closest_chain_z[i]<<"\n";
		}
		
		if(foundChains.size()==1 && !foundvertex)
		{
			prob_vtx_z=c->start_z;
			prob_vtx_t=c->start_t;
			
			if(view==2)foundvertexU=1;
			if(view==3)foundvertexV=1;
			if(foundvertexU && foundvertexV) foundvertex=1;
			vtx_z=prob_vtx_z;
			if(view==2)vtx_u=prob_vtx_t;
			if(view==3)vtx_v=prob_vtx_t;
		}
	
		MSG("PrimaryShowerFinder",Msg::kDebug)<<"ENDING WITH "<<ch->finished.size()<<" CHAINS\n";
		MSG("PrimaryShowerFinder",Msg::kDebug)<<"done\n";

	//}
	
	
	
	}

	MSG("PrimaryShowerFinder",Msg::kDebug)<<"done making chains\n";
	
	mycm->ClearInUse();
}


HoughLine  * PrimaryShowerFinder::SplitHoughLine(Chain *c,HoughLine *hl, Managed::ClusterManager *mycm)
{
	//here we want to split the hough line if the chain bisects it... event by the chain projectionion
	
	if(!mycm)mycm=cluster_manager;
	
	if(!c || !hl)return 0;
	
	//check for intersection?
	int intersection=0;
	double intersection_z=0;
	double intersection_t=0;
	
	
	
	//consider the various intersections:
	
//	printf("looking for intersection chain %f %f %f %f hl %f %f %f %f\n",c->start_z,c->start_t,c->end_z,c->end_t,hl->start_z,hl->start_t,hl->end_z,hl->end_t);
	
	//cross inside
/*	if(c->z.size()<1)return 0;
	for(unsigned int i=0;i<c->z.size()-1;i++)
	{
		double et1 = hl->GetExpectedT(c->z[i]);
		double et2 = hl->GetExpectedT(c->z[i+1]);
		
		double max_t = et1>et2?et1:et2;
		double min_t = et1<et2?et1:et2;
		
		double min_c_t = c->t[i]<c->t[i+1]?c->t[i]:c->t[i+1];
		double max_c_t = c->t[i]>c->t[i+1]?c->t[i]:c->t[i+1];
		
		printf("minct %f mint %f maxct %f maxt %f\n",min_c_t,min_t,max_c_t,max_t);
		
		if( (min_c_t<min_t && max_c_t>max_t) || (max_t>max_c_t && min_t < min_c_t))
		//intersection
		{
			intersection=1;
			intersection_z = (c->z[i]+c->z[i+1])/2;
			intersection_t = (min_c_t+max_c_t)/2;
		
		}
	}
*/	
	if(c->z.size()<1)return 0;
	
	double last_int_z=c->z[0];
	double last_int_dt=c->t[0]-hl->GetExpectedT(c->z[0]);
	
	for(unsigned int i=1;i<c->z.size();i++)
	{
		double et1 = hl->GetExpectedT(c->z[i]);

		double dt=c->t[i]-et1;	
			
//		printf("last z %f dt %f  now z %f dt %f\n",last_int_z,last_int_dt,c->z[i],dt);
		
		//intersection
		if( ( last_int_dt<0 && dt >0 ) || ( last_int_dt >0 && dt < 0) )
		{
			intersection=1;

/*
			intersection_z = c->z[i];
			intersection_t = c->t[i];
*/

			//find the actual intersection
			
			double cslope=(c->z[i]-c->z[i-1]) ? (c->t[i]-c->t[i-1])/(c->z[i]-c->z[i-1]) : 0;
			double coffset=c->t[i]-cslope*c->z[i];
			
			double ht2 = hl->GetExpectedT(c->z[i]);
			double ht1 = hl->GetExpectedT(c->z[i-1]);

			double hslope = (c->z[i]-c->z[i-1]) ? (ht2-ht1)/(c->z[i]-c->z[i-1]) : 0;
			double hoffset = ht2 - hslope*c->z[i];
			
			intersection_t = (hslope-cslope) ? (hslope*coffset - cslope*hoffset)/(hslope-cslope) : 0;
			if(cslope)intersection_z = (intersection_t - coffset ) / cslope;
			else intersection_z=0;
			

			break;
		}
		
		last_int_z=c->z[i];
		last_int_dt=dt;
		
	}
	
	
	//cross with projection
	
	
	
	if(!intersection)return 0;

//	printf("intersection found\n");

	double hl_min_z = hl->start_z<hl->end_z?hl->start_z:hl->end_z;
	double hl_max_z = hl->start_z>hl->end_z?hl->start_z:hl->end_z;
	
	if(intersection_z < hl_min_z || intersection_z > hl_max_z)return 0;
	
	HoughLine *newHL=new HoughLine();
	*newHL=*hl;
	newHL->ResetHits();
	
	HoughLine oldHL=*hl;
	oldHL.ResetHits();
	
	//otherwise, split the hough line at the intersection point....
	for(unsigned int i=0;i<hl->cluster_id.size();i++)
	{
		Managed::ManagedCluster *cl = mycm->GetCluster(hl->cluster_id[i]);
		
		if(cl->z-intersection_z<-0.01)  //require a buffer here in case the chain has two hits in the same z plane which don't have exactly same z values
		{
			oldHL.AddCluster(cl);
//			printf("adding to old hl  %f %f\n",cl->z,cl->t);
		}
		else
		{
			newHL->AddCluster(cl);
//			printf("adding to new hl  %f %f\n",cl->z,cl->t);
		}
	}
	
	

	
	//its possible that we made a new houghline which has the same hits as the old houghline, but the old line now is empty....
	//in such a case... just keep the old hough line
	//this is required to prevent infitie and useless recursion
	
	if(oldHL.ncluster<1 && newHL->ncluster>0)
	{
		*hl=*newHL;
		delete newHL;
		return 0;
	}
	
	*hl=oldHL;	
	return newHL;
}



int PrimaryShowerFinder::FindFirstIntersection(int view, double &t, double &z)
{
	std::vector<HoughLine> * hl;
	if(view==2)hl=&houghlinesU;
	else if(view==3)hl=&houghlinesV;
	else return 0;

	TH2D * inth = view==2?intU:intV;
	
	if(inth)delete inth;
	inth=0;
	
	double min_t = view==2?cluster_manager->minu:cluster_manager->minv;
	double max_t = view==2?cluster_manager->maxu:cluster_manager->maxv;
	inth = new TH2D("","",50,cluster_manager->minz,cluster_manager->maxz,50,min_t,max_t);
	inth->SetDirectory(0);

	double best_z = 10000;
	double best_t=0;

	//iterate over each line
	for(unsigned int i=0;i<(*hl).size();i++)
	{
		//nested loop
		for(unsigned int j=0;j<(*hl).size();j++)
		{
			if(i==j)continue;
		
			//computer intersection and save if better
			HoughLine *h1 = &(*hl)[i];
			HoughLine *h2 = &(*hl)[j];
			
			
			z=
			(
				h2->r/sin(h2->theta) + h2->offset_t 
				- h1->r/sin(h1->theta) - h1->offset_t 
				+(cos(h2->theta)/sin(h2->theta))*h2->offset_z 
				-(cos(h1->theta)/sin(h1->theta))*h1->offset_z
			) 
				/
			(
				cos(h2->theta)/sin(h2->theta)
				-
				cos(h1->theta)/sin(h1->theta)
			);
			
			if(z<best_z && z>0)
			{
				best_z=z;
				best_t=	h1->GetExpectedT(best_z);			
			}
			
			//require the intersection to be within 2 planes of the start/end of either houghline
			if(fabs(h1->start_z-z)>0.07 && fabs(h2->start_z-z)>0.07 && fabs(h1->end_z-z)>0.07 && fabs(h2->end_z-z)>0.07)continue;
			
				inth->Fill(z,h1->GetExpectedT(z),h1->ncluster+h2->ncluster);
	
		}
	}
	
	inth->Draw("colz");
	if(best_z >=10000)return 0;


	t=best_t;
	z=best_z;
	return 1;
}




void PrimaryShowerFinder::Bundle(int view)
{
	std::vector<HoughLine> * hl;
	if(view==2)hl=&houghlinesU;
	else if(view==3)hl=&houghlinesV;
	else return;


	std::vector<std::pair<int,int> >to_merge;
	
	//iterate over each line
	for(unsigned int i=0;i<(*hl).size();i++)
	{
	
//			printf("houghline with e %f clusters %d (t,z) (%f,%f) (%f,%f) chi2/ndf %f  offz %f expt %f phi %f\n", (*hl)[i].sum_e, (*hl)[i].ncluster, (*hl)[i].start_t, (*hl)[i].start_z, (*hl)[i].end_t, (*hl)[i].end_z, (*hl)[i].chi2/(*hl)[i].ncluster,(*hl)[i].offset_z,(*hl)[i].GetExpectedT((*hl)[i].offset_z),(*hl)[i].phi);
			
		//nested loop
		for(int j=i+1;j<(int)(*hl).size();j++)
		{

			//if(i==j)continue;
			
			//did we already match this chain?
			int done=0;
			for(int k=0;k<(int)to_merge.size();k++)
			{
				if(to_merge[k].first==j || to_merge[k].second==j)
				{
					done=1;
					break;
				}
			}
			if(done)continue;
			
		//
//			printf("\t %d ccc e %f clusters %d (t,z) (%f,%f) (%f,%f) chi2/ndf %f  offz %f expt %f phi %f\n", j,(*hl)[j].sum_e, (*hl)[j].ncluster, (*hl)[j].start_t, (*hl)[j].start_z, (*hl)[j].end_t, (*hl)[j].end_z, (*hl)[j].chi2/(*hl)[j].ncluster,(*hl)[j].offset_z,(*hl)[j].GetExpectedT((*hl)[j].offset_z),(*hl)[j].phi);
		
			if(fabs((*hl)[i].phi - (*hl)[j].phi)<0.05)
			{
				if(fabs( (*hl)[i].GetExpectedT((*hl)[i].offset_z) - (*hl)[j].GetExpectedT((*hl)[i].offset_z)) < 0.05)
				{
					to_merge.push_back(make_pair(i,j));
//					printf("tm %d %d %f %f %f\n",i,j,(*hl)[i].offset_z,(*hl)[i].GetExpectedT((*hl)[i].offset_z),(*hl)[j].GetExpectedT((*hl)[i].offset_z));
				
				}
			
			
			}
	
		}
	}
	
	if(to_merge.size()<1)return;
	
	std::vector<HoughLine> toKeep;

	//keep non merged lines
	for(unsigned int i=0;i<hl->size();i++)
	{
		int keep=1;
		for(unsigned int j=0;j<to_merge.size();j++)
		{
			if(to_merge[j].first ==(int)i || to_merge[j].second ==(int)i)
			{
				keep=0;
				break;
			}
		}
		if(keep){
			toKeep.push_back((*hl)[i]);
//			printf("not merging %d\n",i);
		}
	
	}


	for(unsigned int i=0;i<to_merge.size();i++)
	{
/*
		//simple for now... just pick one out of the group
		if(i==0 || (i>0 && (to_merge[i].first != to_merge[i-1].first) ) )
		{

			toKeep.push_back((*hl)[to_merge[i].first]);
			
			printf("%d \n",to_merge[i].first);
		}
*/
		std::vector<int> mlist;
		
		mlist.push_back(to_merge[i].first);
		for(unsigned int j=i+1;j<to_merge.size();j++)
		{
			if(to_merge[j].first !=to_merge[i].first)break;
			mlist.push_back(to_merge[j].second);
		}


		//now we have a list of houghlines to merge...
		//do a weighted average by energy
		
		double theta=0;
		double r=0;
		double offset_t=0;
		double offset_z=0;
		double sum_e=0;
		
		for(unsigned int j=0;j<mlist.size();j++)
		{
			HoughLine *h = &(*hl)[mlist[j]];
			theta += h->theta*h->sum_e;
			r += h->r*h->sum_e;
			offset_t += h->offset_t*h->sum_e;
			offset_z += h->offset_z*h->sum_e;
			sum_e+=h->sum_e;
		}
		
		if(sum_e<0.1)continue;
		theta/=sum_e;
		r/=sum_e;
		offset_t/=sum_e;
		offset_z/=sum_e;
		
		HoughLine l( theta,  r,  offset_t,  offset_z);
		LoadCloseHits(&l,cluster_manager,view,0.0412);
		toKeep.push_back(l);

	}



	
	(*hl)=toKeep;
	
//	printf("# of hough lines %d\n",hl->size());

}


void PrimaryShowerFinder::CleanHoughLines(int view,std::vector<HoughLine>*hhh=0)
{

	
	std::vector<HoughLine> * hl;
	if(!hhh)
	{
		if(view==2)hl=&houghlinesU;
		else if(view==3)hl=&houghlinesV;
		else return;
	}else{
		hl=hhh;
	}
	
	
	
	std::vector<int>todel;
	//remove all but first vertical hough lines in this view
	
	//first hit z is:
	double first_z = 1000;
	int keep_v=-1;
	for(unsigned int i=0;i<hl->size();i++)
	{
		if(fabs((*hl)[i].phi-3.141592/2 )<0.1){
			if((*hl)[i].start_z < first_z)
			{
				first_z = (*hl)[i].start_z;
				keep_v=i;
			}
		}
	}			
		
	for(int i=0;i<(int)hl->size();i++)
	{
		if(fabs((*hl)[i].phi-3.141592/2 )<0.1){
			if(keep_v!=i)todel.push_back(i);
		}
	}	
	
		
			
	//remove hough lines that have the closest distance between a pair of adjacent hits larger than some value
	double min_dist = 0.2;
	for(unsigned int i=0;i<hl->size();i++)
	{
	

		double small_dist=1000;
		int done=0;
		for(unsigned int j=0;j<(*hl)[i].cluster_id.size();j++)
		{
			
			for(unsigned int k=j+1;k<(*hl)[i].cluster_id.size();k++)
			{
				int id1 = (*hl)[i].cluster_id[j];
				int id2 = (*hl)[i].cluster_id[k];
				
				Managed::ManagedCluster * c1 = cluster_manager->GetCluster(id1);
				Managed::ManagedCluster * c2 = cluster_manager->GetCluster(id2);
				if(!c1 || !c2)continue;
				
			//	double dist = sqrt((c1->t-c2->t)*(c1->t-c2->t) + (c1->z-c2->z)*(c1->z-c2->z));
				
			//	//adjust the distance based on the angle of the houghline
			//	dist *=cos((*hl)[i].phi);
			
				double dist = fabs(c1->z-c2->z);
				
				if(dist < small_dist)small_dist=dist;		
				
		//		printf("dist %f\n",dist);		
				
				if(small_dist < min_dist)
				{
					done=1;
					break;
				}
			}
			if(done)break;	
		}
		if(done)continue;
		
		//if we are here.. we want to delete it
		todel.push_back(i);
		
	}


	//hit density
	for(unsigned int i=0;i<hl->size();i++)
	{
		//double len = sqrt(((*hl)[i].end_z-(*hl)[i].start_z)*((*hl)[i].end_z-(*hl)[i].start_z) + ((*hl)[i].end_t-(*hl)[i].start_t)*((*hl)[i].end_t-(*hl)[i].start_t) );
		
		double len=fabs((*hl)[i].end_z-(*hl)[i].start_z);
		double frac = len>0?(*hl)[i].ncluster/(len/0.07):0;
		
		if(frac <0.5)todel.push_back(i);
	}	



	//need to save?
	if(todel.size()<1)return;
	
	std::vector<HoughLine> temp;
	for(int i=0;i<(int)hl->size();i++)
	{
		int save=1;
		for(int j=0;j<(int)todel.size();j++)
		{
			if(todel[j]==i)
			{
				save=0;
				break;
			}
		}	
		if(!save)continue;
		temp.push_back((*hl)[i]);
	}

	*hl=temp;


}


Chain * PrimaryShowerFinder::ExtractMuon(Chain *c)
{
	if(!c)return 0;
	
	int required_consecutive_muon_hits=2;
	
	if(c->entries<required_consecutive_muon_hits)return 0;
	
	MSG("PrimaryShowerFinder",Msg::kDebug)<<"extract muon\n";
	
	//do we have a muon in this chain?
	
	double muoncount=0;
	double muonenergy=0;
	double sume=0;
	
	for(int i=0;i<c->entries;i++)
	{
		if(c->e[i]<2.5)
		{
			muoncount++;
			muonenergy+=c->e[i];
		}
		sume+=c->e[i];
	
	}

	double muonfrac = muoncount / c->entries;

	//double efrac=muonenergy/sume;
	
	MSG("PrimaryShowerFinder",Msg::kDebug)<<"mf "<<muonfrac<<"\n";
	if(muonfrac<0.2 ) return 0;//no muon
	
	
	
	//find where the muon exists the other particles... require n in a row muon hits
	
	int startmuon=0;
	int nmuon=0;
	
	for(int i=0;i<c->entries;i++)
	{
	
	//	printf("%d %f\n",i,c->e[i]);
		startmuon=i;
		if(c->e[i]<2.5)
		{
			nmuon++;
			if(nmuon>=required_consecutive_muon_hits)break;  //require n muon hits in a row
		
		}else
		{
			nmuon=0;
		}
	}
	
	if(nmuon<required_consecutive_muon_hits)return 0;//not enough muon like!
	//printf("sm %d\n",startmuon);
	startmuon=startmuon-nmuon+1;
	
	MSG("PrimaryShowerFinder",Msg::kDebug)<<"muon starts at "<<startmuon<<"\n";
	
	//determine average muon energy in this particle
	//using the muon hits, but ignoring large ones which could be muon->e shower
	
	double muone=0;
	muoncount=0;
	for(int i=startmuon;i<c->entries;i++)
	{
		if(c->e[i]<2.5)
		{
			muone+=c->e[i];
			muoncount++;
		}
	}
	
	double adjmuon = muone/muoncount;
	
	MSG("PrimaryShowerFinder",Msg::kDebug)<<"avg muon e "<<adjmuon<<"\n";
	
	Chain ctemp;
	
	std::vector<double> savedE;
	for(int i=0;i<startmuon;i++)
	{
		if(c->e[i]-adjmuon>0)
		{
			ctemp.insert(c->t[i],c->z[i],c->e[i]-adjmuon,c->cluster_id[i]);
			savedE.push_back(adjmuon);
		}else{
			savedE.push_back(0);
		}
		MSG("PrimaryShowerFinder",Msg::kDebug)<<"AAAAAAAAAAAA "<<c->t[i]<<" "<<c->z[i]<<" e "<<c->e[i]-adjmuon<<"\n";
	
	}
	

	Chain *muonchain = new Chain();
	muonchain->level=c->level; //so we get the parend id, etc
	muonchain->parentChain=c->parentChain;
	muonchain->children=c->children;
	
	
	for(int i=0;i<startmuon;i++)
	{
		if(savedE[i]>0)
			muonchain->insert(c->t[i],c->z[i],savedE[i],c->cluster_id[i]);
	}	
	
	for(int i=startmuon;i<c->entries;i++)
	{
		muonchain->insert(c->t[i],c->z[i],c->e[i],c->cluster_id[i]);
	}		



	*c=ctemp;
	
	return muonchain;

}






int PrimaryShowerFinder::FindPrimaryShower(Managed::ClusterManager *cm)
{

	MSG("PrimaryShowerFinder",Msg::kDebug)<<"looking for primary shower\n";
	Reset(0); //keep any vertex that we might have set

	ran=1;

	cluster_manager=cm;
	
	
	MakeHoughMap(2);
	MakeHoughMap(3);

//DumpHoughLines(2);
//DumpHoughLines(3);

	CleanHoughLines(2);
	CleanHoughLines(3);
	
	Bundle(2);
	Bundle(3);

	MakeChains(2);
	MakeChains(3);


//	FindHoughLineMatches();
	
	
	//if(houghlineMatch.size()<1)
		
	
	
	//try to find a long Shower chain in each view
	//chain_u = FindShowerChain(cm,2);
	//chain_v = FindShowerChain(cm,3);
	
	
//	chain_u = MakeShowerChain(2);
//	chain_v = MakeShowerChain(3);


	chain_u =0;
	chain_v=0;
	
	for(unsigned int i=0;i<chu->finished.size();i++)
	{
		if(chu->finished[i].available)
		{
			chain_u=&chu->finished[i];
			break;
		}
	}
	
	for(unsigned int i=0;i<chv->finished.size();i++)
	{
		if(chv->finished[i].available)
		{
			chain_v=&chv->finished[i];
			break;
		}
	}

	if(!chain_u || !chain_v)return 0;
	

//make temp chains for use in the finder which represent a max path
	


	//printf("building chain\n");
	chu->print_finished();	
	std::pair< std::vector<int>, double> path;
	path = chu->FindMaxPath(chain_u);	
	chain_u=new Chain();


	//printf("chain made from ");
//	for(int i=0 ;i <path.first.size();i++)
//		printf("%d ",path.first[i]);
//	printf("\n");
		


	for(int i=0 ;i <(int)path.first.size();i++)
	{
		Chain *c = chu->GetChain(path.first[i]);
		if(!c)continue;
		for(int j=0;j<c->entries;j++)
		{
			chain_u->insert(c->t[j],c->z[j],c->e[j],c->cluster_id[j]);
			MSG("PrimaryShowerFinder",Msg::kDebug)<<"u "<<c->z[j]<<" "<<c->t[j]<<" "<<c->e[j]<<"\n";
		}
	}
	std::vector<int> upath=path.first;

	//printf("building chain\n");
	chv->print_finished();
	path = chv->FindMaxPath(chain_v);	
	chain_v=new Chain();
	
	//printf("chain made from ");
	//	for(int i=0 ;i <path.first.size();i++)
	//		printf("%d ",path.first[i]);
	//	printf("\n");
		

	for(int i=0 ;i <(int)path.first.size();i++)
	{
		Chain *c = chv->GetChain(path.first[i]);
		if(!c)continue;
//		printf("chain %d\n",path.first[i]);
		for(int j=0;j<c->entries;j++)
		{
//			printf("id %d\n",c->cluster_id[j]);
			chain_v->insert(c->t[j],c->z[j],c->e[j],c->cluster_id[j]);
			MSG("PrimaryShowerFinder",Msg::kDebug)<<"v "<<c->z[j]<<" "<<c->t[j]<<" "<<c->e[j]<<"\n";
		}
	}
	std::vector<int> vpath=path.first;		
	
/*	
	printf("\n\n---------\nchain\n");
	for(int i=0;i<chain_u->entries;i++)
	{
		printf("%d %f %f %f\n",chain_u->cluster_id[i],chain_u->z[i],chain_u->t[i],chain_u->e[i]);
	}

	printf("\n\n---------\nchain\n");
	for(int i=0;i<chain_v->entries;i++)
	{
		printf("%d %f %f %f\n",chain_v->cluster_id[i],chain_v->z[i],chain_v->t[i],chain_v->e[i]);
	}	
*/	
	//muon extraction:
	/////////
	Chain * mchainU=0;
	Chain * mchainV=0;
	mchainU=ExtractMuon(chain_u);
	mchainV=ExtractMuon(chain_v);
	////////////
/*	
	printf("\n\n---------\nmuon chain\n");
	if(mchainU)for(int i=0;i<mchainU->entries;i++)
	{
		printf("%d %f %f %f\n",mchainU->cluster_id[i],mchainU->z[i],mchainU->t[i],mchainU->e[i]);
	}

	printf("\n\n---------\nmuon chain\n");
	if(mchainV)for(int i=0;i<mchainV->entries;i++)
	{
		printf("%d %f %f %f\n",mchainV->cluster_id[i],mchainV->z[i],mchainV->t[i],mchainV->e[i]);
	}	
	

	printf("\n\n---------\nchain\n");
	for(int i=0;i<chain_u->entries;i++)
	{
		printf("%d %f %f %f\n",chain_u->cluster_id[i],chain_u->z[i],chain_u->t[i],chain_u->e[i]);
	}

	printf("\n\n---------\nchain\n");
	for(int i=0;i<chain_v->entries;i++)
	{
		printf("%d %f %f %f\n",chain_v->cluster_id[i],chain_v->z[i],chain_v->t[i],chain_v->e[i]);
	}	
	*/
	
	//double ustart=chain_u->start_z;
	//double vstart=chain_v->start_z;
	//double uend=chain_u->end_z;
	//double vend=chain_v->end_z;
//	ExpandShowerChain(chain_u,2,vstart,vend);
//	ExpandShowerChain(chain_v,3,ustart,uend);
	
	
	int qual_u=0;
	int qual_v=0;
	
	if(chain_u)
	{
		chain_u->Recalc();
		qual_u=CheckChainQuality(chain_u);
	}
	if(chain_v)
	{
		chain_v->Recalc();
		qual_v=CheckChainQuality(chain_v);
	}

	if(qual_u&&!qual_v)
	{
		if(chain_u->end_z - chain_u->start_z > 2)single_view_long_shower=1;
	}else if(qual_v&&!qual_u)
	{
		if(chain_v->end_z - chain_v->start_z > 2)single_view_long_shower=1;
	}
	
	if(!qual_u || !qual_v)return 0;



	//should do something smarter... like try other chains...
	
	//if we have a quality chain in each view... then we will be making a particle
/*	if(qual_u && qual_v)
	{
	
		MSG("PrimaryShowerFinder",Msg::kDebug)<<"qqqqq "<<qual_u <<" "<< qual_v<<"\n";
		
	
		//look in both views to find at least 3u and 3v consecutive planes with only 1 strip... that is where we say the Shower exits the rest of the shower/partices/etc

		int qual=CheckChainOverlap(chain_u, chain_v);
		if(qual)
		{
		
		ClearFrontVertex(chain_u, chain_v);
*/

		
		
	
	
		MSG("PrimaryShowerFinder",Msg::kDebug)<<"making particle\n";
		MakeParticle3D();
		if(foundparticle3d)
		{
			foundparticle=1;
		
			foundparticle3d->particletype=Particle3D::electron;
		
			MSG("PrimaryShowerFinder",Msg::kDebug)<<"particle made\n";
			DumpParticle3D();
			
//			chain_u->available=0;
//			chain_v->available=0;

			for(int i=0 ;i <(int)upath.size();i++)
			{
				Chain *c = chu->GetChain(upath[i]);
				
				MSG("PrimaryShowerFinder",Msg::kDebug)<<"u chain "<<upath[i]<<"\n";
				
				c->available=0;
			}
			for(int i=0 ;i <(int)vpath.size();i++)
			{
				Chain *c = chv->GetChain(vpath[i]);
				
				MSG("PrimaryShowerFinder",Msg::kDebug)<<"v chain "<<vpath[i]<<"\n";
				
				c->available=0;
			}			
			
			//do a fit on the particle3d to see if it is like an electron..
			Particle3D *p=foundparticle3d;		
			ShwFit f;
			
			double startz=0;
			
			startz=p->z[0];
			
			for(int i=0;i<p->entries;i++)
			{
				//double planewidth=0.035;
				
		//		f.Insert(p->e[i],(int)((p->z[i]-startz)/planewidth));
			
		//		cout <<"inserting "<<p->e[i]<<" at plane "<<(int)((p->z[i]-startz)/planewidth)<<endl;
		
		//for now, do the simple thing and assume each hit is in the next plane...
			
				//			f.Insert(p->e[i],i+1);
			
			
				f.Insert3d(p->e[i],p->u[i],p->v[i],p->z[i]);
			
				//cout <<"inserting "<<p->e[i]<<" at plane "<<i+1<<endl;
		
			
			}
			
			//f.Fit();
		

			f.Fit3d(1);




			MSG("PrimaryShowerFinder",Msg::kDebug)<<"EM FIT! " << f.par_a << " " << f.par_b << " +- "<<f.par_b_err<<" "<<f.par_e0<<endl;
		
			p->emfit_a=f.par_a;
			p->emfit_b=f.par_b;
			p->emfit_e0=f.par_e0;
			p->emfit_a_err=f.par_a_err;
			p->emfit_b_err=f.par_b_err;
			p->emfit_e0_err=f.par_e0_err;
			p->emfit_prob=f.prob;
			p->calibrated_energy=f.par_e0;
			p->emfit_chisq=f.chisq;
			p->emfit_ndf=f.ndf;

			p->pred_e_a=f.pred_e_a;	
			p->pred_g_a=f.pred_g_a;
			p->pred_b=f.pred_b;
			p->pred_e0=f.pred_e0;
			p->pred_e_chisq=f.pred_e_chisq;
			p->pred_e_ndf=f.pred_e_ndf;
			p->pred_g_chisq=f.pred_g_chisq;
			p->pred_g_ndf=f.pred_g_ndf;				
			
			p->pre_over=f.pre_over;		
			p->pre_under=f.pre_under;		
			p->post_over=f.post_over;		
			p->post_under=f.post_under;		

			p->pp_chisq=f.pp_chisq;
			p->pp_ndf=f.pp_ndf;
			p->pp_igood=f.pp_igood;
			p->pp_p=f.pp_p;
	
			p->cmp_chisq = f.cmp_chisq;
			p->cmp_ndf = f.cmp_ndf;
			p->peakdiff = f.peakdiff;	
			
			int redochain=0; 
			//do we need to redo the chain because we are shifting the start past the first hits of the chain?
			if(f.zshift)
			{
				if(f.zshift < 0)redochain=1;
				p->start_z -=f.zshift;
				
				MSG("PrimaryShowerFinder",Msg::kDebug)<<"zshift "<<f.zshift<<"\n";
			}
			
			
			if(!f.conv)p->particletype=Particle3D::other;
 			else{
 			
 			
 			
 			
 			//redo the chain if needed
 			if(redochain)
 			{

				MSG("PrimaryShowerFinder",Msg::kDebug)<<"REDOING CHAINS\n";
 			
				Chain temp_u;
 				for(int i=0;i<chain_u->entries;i++)
				{
					if(chain_u->z[i]<p->start_z){
						MSG("PrimaryShowerFinder",Msg::kDebug)<<"throwing away chain hit with z "<<chain_u->z[i]<<"\n";
						continue;
					}
					temp_u.insert(chain_u->t[i],chain_u->z[i],chain_u->e[i],chain_u->cluster_id[i]);
				}
				temp_u.myId=chain_u->myId;
				*chain_u=temp_u;			

				Chain temp_v;
				for(int i=0;i<chain_v->entries;i++)
				{
					if(chain_v->z[i]<p->start_z){
						MSG("PrimaryShowerFinder",Msg::kDebug)<<"throwing away chain hit with z "<<chain_v->z[i]<<"\n";
						continue;
					}
					temp_v.insert(chain_v->t[i],chain_v->z[i],chain_v->e[i],chain_v->cluster_id[i]);
				}
				temp_v.myId=chain_v->myId;		
				*chain_v=temp_v;
 				chain_u->available=0;
				chain_v->available=0;

				delete foundparticle3d; 
				foundparticle=0;	
				foundparticle3d=0;
					
 				//remake particle
 				MakeParticle3D();

				if(!foundparticle3d)return 0;

				if(foundparticle3d)
				{
				foundparticle=1;	
				p=foundparticle3d;	
				
				f.Reset();
				
				double startz=0;
			
				startz=p->z[0];
			
				for(int i=0;i<p->entries;i++)
				{
					//double planewidth=0.035;
					f.Insert3d(p->e[i],p->u[i],p->v[i],p->z[i]);
				}
			

				
				foundparticle3d->particletype=Particle3D::electron;
				f.Fit3d();
			
				MSG("PrimaryShowerFinder",Msg::kDebug)<<"EM FIT! " << f.par_a << " " << f.par_b << " +- "<<f.par_b_err<<" "<<f.par_e0<<endl;
		
				p->emfit_a=f.par_a;
				p->emfit_b=f.par_b;
				p->emfit_e0=f.par_e0;
				p->emfit_a_err=f.par_a_err;
				p->emfit_b_err=f.par_b_err;
				p->emfit_e0_err=f.par_e0_err;
				p->emfit_prob=f.prob;
				p->calibrated_energy=f.par_e0;
				p->emfit_chisq=f.chisq;
				p->emfit_ndf=f.ndf;

				p->pred_e_a=f.pred_e_a;	
				p->pred_g_a=f.pred_g_a;
				p->pred_b=f.pred_b;
				p->pred_e0=f.pred_e0;
				p->pred_e_chisq=f.pred_e_chisq;
				p->pred_e_ndf=f.pred_e_ndf;
				p->pred_g_chisq=f.pred_g_chisq;
				p->pred_g_ndf=f.pred_g_ndf;				
			
				p->pre_over=f.pre_over;		
				p->pre_under=f.pre_under;		
				p->post_over=f.post_over;		
				p->post_under=f.post_under;	
				
				p->pp_chisq=f.pp_chisq;
				p->pp_ndf=f.pp_ndf;
				p->pp_igood=f.pp_igood;
				p->pp_p=f.pp_p;
	
				p->cmp_chisq = f.cmp_chisq;
				p->cmp_ndf = f.cmp_ndf;
				p->peakdiff = f.peakdiff;	
								
			
				if(!f.conv)p->particletype=Particle3D::other;
				if(f.zshift)
				{
					p->start_z -=f.zshift;
					MSG("PrimaryShowerFinder",Msg::kDebug)<<"zshift "<<f.zshift<<"\n";
				}
			
 				}
 			
 			}
 			
 			
 		}	

 			
 		//save clusters involved...

		//don't save until we know that we want it!

		for(int i=0;i<chain_u->entries;i++)
		{
			int newid = cluster_manager->SaveCluster(chain_u->cluster_id[i],chain_u->e[i],2);
			if(!newid)continue;//if we have an issue saving, such as there is no energy...
			Managed::ManagedCluster * newcluster = cluster_manager->GetCluster(newid);
			if(fabs(chain_u->e[i]-newcluster->e)>0.001)
				MSG("PrimaryShowerFinder",Msg::kDebug)<<"saving with different e "<<chain_u->e[i]<<" "<<newcluster->e<<"\n";
			chain_u->e[i]=newcluster->e;
			chain_u->cluster_id[i]=newid;
		}



		for(int i=0;i<chain_v->entries;i++)
		{
			int newid = cluster_manager->SaveCluster(chain_v->cluster_id[i],chain_v->e[i],2);
			if(!newid)continue;//if we have an issue saving, such as there is no energy...
			Managed::ManagedCluster * newcluster = cluster_manager->GetCluster(newid);
			
			//printf("saving cluster id %d e %f newcluster ok ? %d\n",chain_v->cluster_id[i],chain_v->e[i],newcluster!=0);
			//if(newcluster)printf("new cluster id %d  actual %d\n",newid,newcluster->id);
			if(fabs(chain_v->e[i]-newcluster->e)>0.001)
				MSG("PrimaryShowerFinder",Msg::kDebug)<<"saving with different e "<<chain_v->e[i]<<" "<<newcluster->e<<"\n";
			chain_v->e[i]=newcluster->e;
			chain_v->cluster_id[i]=newid;
		}		

 			
 		//save chains involved!
 		chain_u->available=0;
 		chain_v->available=0;

 		//printf("psf wants to add %d %d\n",chain_u->myId,chain_v->myId);
 		chu->AddFinishedChain(*chain_u);
 		chv->AddFinishedChain(*chain_v);
 		

 		
 		//did we have muon chains?
 		if(mchainU)chu->AddFinishedChain(*mchainU);	
		if(mchainV)chv->AddFinishedChain(*mchainV);		
 			
 		//we made the electron chain temporary....
 		//now remove the actual chains that were used to make the electron and muon extracted chains
 		
 		
 		for(unsigned int i=0;i<upath.size();i++)
 		{
 			chu->DeleteChain(upath[i]);
 		}
 
  		for(unsigned int i=0;i<vpath.size();i++)
 		{
 			chv->DeleteChain(vpath[i]);
 		}
 				

		
		//!!!! we also need to remake the particle here after adjusting chain energy!
 		//////////////
 		 				//remake particle
 
 				delete foundparticle3d; 
				foundparticle=0;	
				foundparticle3d=0;
					
 				//remake particle
 				MakeParticle3D();
				p=foundparticle3d;	
				if(!foundparticle3d)return 0;

				foundparticle=1;
								
				foundparticle3d->particletype=Particle3D::electron;
				f.Fit3d(1);
			
				MSG("PrimaryShowerFinder",Msg::kDebug)<<"EM FIT! " << f.par_a << " " << f.par_b << " +- "<<f.par_b_err<<" "<<f.par_e0<<endl;
		
				p->emfit_a=f.par_a;
				p->emfit_b=f.par_b;
				p->emfit_e0=f.par_e0;
				p->emfit_a_err=f.par_a_err;
				p->emfit_b_err=f.par_b_err;
				p->emfit_e0_err=f.par_e0_err;
				p->emfit_prob=f.prob;
				p->calibrated_energy=f.par_e0;
				p->emfit_chisq=f.chisq;
				p->emfit_ndf=f.ndf;

				p->pred_e_a=f.pred_e_a;	
				p->pred_g_a=f.pred_g_a;
				p->pred_b=f.pred_b;
				p->pred_e0=f.pred_e0;
				p->pred_e_chisq=f.pred_e_chisq;
				p->pred_e_ndf=f.pred_e_ndf;
				p->pred_g_chisq=f.pred_g_chisq;
				p->pred_g_ndf=f.pred_g_ndf;				
			
				p->pre_over=f.pre_over;		
				p->pre_under=f.pre_under;		
				p->post_over=f.post_over;		
				p->post_under=f.post_under;		

				p->pp_chisq=f.pp_chisq;
				p->pp_ndf=f.pp_ndf;
				p->pp_igood=f.pp_igood;
				p->pp_p=f.pp_p;
	
	
				p->cmp_chisq = f.cmp_chisq;
				p->cmp_ndf = f.cmp_ndf;
				p->peakdiff = f.peakdiff;	
				
			
				if(!f.conv)p->particletype=Particle3D::other;
 		
 		//////////////	
 			
 			
 			}
			
		
//		}
//	}



	MSG("PrimaryShowerFinder",Msg::kDebug)<<"printing long Shower chains\nu\n";
	if(chain_u)chain_u->PrintChain();MSG("PrimaryShowerFinder",Msg::kDebug)<<"v\n";
	if(chain_v)chain_v->PrintChain();
	MSG("PrimaryShowerFinder",Msg::kDebug)<<"done\n";
	

	MSG("PrimaryShowerFinder",Msg::kDebug)<<"done looking for long Showers\n";
	return foundparticle;
}



//once we make chains in both view, we need to see if one chain is too short compared to the other
//if so, try to extend it by finding near hits
void PrimaryShowerFinder::ExpandShowerChain(Chain * chain,int view, double start_z, double end_z)
{

	if(chain->start_z < start_z && chain->end_z > end_z)return;
	
//	printf("expanding shower chain from %f to %f\n",chain->start_z,chain->end_z);
				
	if(chain->start_z > start_z)
	{
		double z=chain->start_z-0.06;
		while(z>start_z-0.04)
		{
			//find hits in this z range and add the best one to the chain
			std::vector<int> cs = cluster_manager->FindClustersInZ(z,view);
	
	//		printf("found %d hits in z %f\n",cs.size(),z);
			
			int closest=-1;
			double closest_d=1000;
			//for now do a simple straight line from start t
			for(unsigned int i=0;i<cs.size();i++)	
			{
				Managed::ManagedCluster *c = cluster_manager->GetCluster(cs[i]);
	//			printf("id %d z %f t %f e %f\n",cs[i],c->z,c->t,c->e);
			
				if(fabs(c->t-chain->start_t)<closest_d)
				{
					closest_d=fabs(c->t-chain->start_t);
					closest=c->id;
				}
			}
			
			if(closest>-1 && closest_d < 0.1)
			{
				Managed::ManagedCluster *c = cluster_manager->GetCluster(closest);
				chain->insert(c->t,c->z,c->e,c->id);
			}
	
			z-= 0.06; //subtract off a plane
		}
	}


		
	if(chain->end_z < end_z)
	{
		double z=chain->end_z+0.06;
		while(z<end_z+0.04)
		{
			//find hits in this z range and add the best one to the chain
			std::vector<int> cs = cluster_manager->FindClustersInZ(z,view);
	
	//		printf("found %d hits in z %f\n",cs.size(),z);
			
			int closest=-1;
			double closest_d=1000;
			//for now do a simple straight line from start t
			for(unsigned int i=0;i<cs.size();i++)	
			{
				Managed::ManagedCluster *c = cluster_manager->GetCluster(cs[i]);
	//			printf("id %d z %f t %f e %f\n",cs[i],c->z,c->t,c->e);
			
				if(fabs(c->t-chain->end_t)<closest_d)
				{
					closest_d=fabs(c->t-chain->end_t);
					closest=c->id;
				}
			}
			
			if(closest>-1  && closest_d < 0.1)
			{
				Managed::ManagedCluster *c = cluster_manager->GetCluster(closest);
				chain->insert(c->t,c->z,c->e,c->id);
			}
	
			z+= 0.06; //subtract off a plane
		}
	}


}




//make sure that the chains both start within one plane of each other... 
//for now, simply do this by finding the plane below the most upstream hit
void PrimaryShowerFinder::ClearFrontVertex(Chain * chain_u, Chain * chain_v)
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
int PrimaryShowerFinder::CheckChainOverlap(Chain * chain_u, Chain * chain_v)
{

	

	
	//now make sure that there is sufficient overlap
	double require_overlap_distance=1.0;
	
	double start_z = chain_u->start_z > chain_v->start_z ? chain_u->start_z : chain_v->start_z;
	double end_z = chain_u->end_z < chain_v->end_z ? chain_u->end_z : chain_v->end_z;

	MSG("PrimaryShowerFinder",Msg::kDebug)<<"overlap "<<end_z-start_z<<"\n";
	
	if(end_z-start_z < require_overlap_distance)return 0;
	return 1;
	
	

}


// check to see if this is a Shower-like chain with a sufficient number of hits 
int PrimaryShowerFinder::CheckChainQuality(Chain *c)
{
	if(!c)return 0;
	if(c->entries<1)return 0;	
	
	int qual=1;
	
	//look at Shower-like strip fraction
	//and sparsity (# planes with hits/#planes from front to back of Shower)
	
	int Showerlike=0;
	int hitplane=0;
	
	double last_hitplane=-100;
	
	double min = 1000000;
	double max = 0;
	for(int i=0;i<c->entries;i++)
	{
		if(c->e[i]>0.5 && c->e[i]<2.5)Showerlike++;
		if(fabs(last_hitplane - c->z[i])>0.03)
		{
			hitplane++;
			last_hitplane=c->z[i];
		}
		
		min = c->z[i]< min? c->z[i]:min;
		max = c->z[i]> max? c->z[i]:max;
	}
	
	double Showerfrac = (double)Showerlike / (double)c->entries;
	double sparsity = max-min?(double)hitplane / (max-min) / 7.08 :1; //length/size of u+v plane
	

//definition of shwerfrac is wrong!	
//	if (Showerfrac < 0.5) qual=0;
	if (sparsity < 0.8) qual=0;
	
	if(c->muonfrac>0.5)qual=0;
	
	
	MSG("PrimaryShowerFinder",Msg::kDebug)<<"MUON FRAC "<<c->muonfrac <<"Shower chain check  Showerfrac "<< Showerfrac<<" sparsity "<< sparsity<<"  qual "<< qual<<"\n";

	return qual;
}


Chain * PrimaryShowerFinder::FindShowerChain(Managed::ClusterManager *cl, int view)
{

	MSG("PrimaryShowerFinder",Msg::kDebug)<<"\n\nLooking for long Shower Chain in view "<<view<<"\n";


	std::map<double, std::map<double, int> > * cluster_map = cl->GetClusterMap(view);


    std::map<double, std::map<double, int> >::reverse_iterator p_iterr;
    std::map<double, int >::iterator s_iterr;


	Chain * c = new Chain();

	Managed::ManagedCluster * largest=0;
	double last_plane=100000;
	
	std::vector<Managed::ManagedCluster *> maxplaneclusters;
	//first we want to find maximum hits in each plane
	//we will then check these hits to see if they can form a tightly aligned track
    for(p_iterr=cluster_map->rbegin();p_iterr!=cluster_map->rend(); p_iterr++)
    {	
    	//new plane and a largest cluster in the previous plane?
		if(fabs(last_plane-p_iterr->first)>0.03 && largest)
		{
			MSG("PrimaryShowerFinder",Msg::kDebug)<<"largest cluster found "<<largest->t<<" "<<largest->z<<" "<<largest->e<<"\n";
			maxplaneclusters.push_back(largest);
			largest=0;
		}
	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cl->GetCluster(s_iterr->second);
			if(!largest)
			{
				largest=clus;
				continue;
			}
			if(largest->e<clus->e)largest=clus;
	  	}
	}
	if(largest)
	{
		MSG("PrimaryShowerFinder",Msg::kDebug)<<"largest cluster found "<<largest->t<<" "<<largest->z<<" "<<largest->e<<"\n";
		maxplaneclusters.push_back(largest);
		largest=0;
	}

	
	for(unsigned int i=0;i<maxplaneclusters.size();i++)
		c->insert(maxplaneclusters[i]->t,maxplaneclusters[i]->z,maxplaneclusters[i]->e,maxplaneclusters[i]->id);		
			
	return c;			


}


void PrimaryShowerFinder::DumpParticle3D()
{
	ostringstream s;
	
	s << "Long Shower Particle3D found "<<endl;

		Particle3D * p = foundparticle3d;
		if (p==0)return;
	
	

	
		s << "\n---Particle3d "  << "---\nstart, stop (u,v,z) (" << p->start_u << ", "<< p->start_v << ", " << p->start_z <<") (" <<p->end_u<<", "<< p->end_v << ", " << p->end_z <<")"<<endl;
		s << "entries "<< p->entries << "  muonlike " << p->muonfrac <<endl;
		

		
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
					s<<"uniqueShower ";
					break;	
			
			}
				
		}
		s<<endl;
		
		s<<"particletype : "<<p->particletype<<endl;
		
		
		s<<"par a " << p->emfit_a << " par b "<< p->emfit_b << " par e0 "<< p->emfit_e0 << " cale "<<p->calibrated_energy<<" chisq " << p->emfit_chisq<<" ndf " << p->emfit_ndf<<endl;
		s<<"emprob " << p->emfit_prob <<" avg rms_t "<<p->avg_rms_t<<endl;
				
		s << "spoints (u,v,z,e - chain, chainhit, chainview - shared - rms_t - view) : ";
		for(int j=0;j<p->entries;j++)
			s << "(" << p->u[j]<<", "<<p->v[j]<<", "<<p->z[j]<<", "<<p->e[j]<<" - " << p->chain[j] <<", "<<p->chainhit[j]<<", "<<p->view[j]<<" - "<<p->shared[j]<<" - "<< p->rms_t[j]<< " - " << p->view[j]<<") ";
	
	/*
	s << "spoints (e - shared) : ";
		for(int j=0;j<p->entries;j++)
			s << "(" <<p->e[j]<<" - "<<p->shared[j]<<") ";
	
	*/
	
	

			
			

			
			
		s<<endl<<endl;
		
		
	
	
	

	MSG("PrimaryShowerFinder",Msg::kDebug) <<s.str();

	
}



void PrimaryShowerFinder::MakeParticle3D()
{

	//printf("making psf particle cu %d cv %d\n",chain_u->myId,chain_v->myId);

	//verify that the chains in both views have sufficient overlap
	
	
	//make the 3d particle from these chains
	if(foundparticle3d)delete foundparticle3d;
	foundparticle3d= new Particle3D();


	
	std::vector<spoint> spoints;
	
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

			
			spoint a;
			a.z=c->z[j];
			a.t=c->t[j];
			a.e=c->e[j];
			a.chain=c->myId;
			a.chainhit=j;
			a.view=2;
			

			Managed::ManagedCluster * clu = cluster_manager->GetCluster(c->cluster_id[j]);
		
			if(clu)
			{
				a.rms_t=clu->rms_t;
				spoints.push_back(a);
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
			

			spoint a;
			a.z=c->z[j];
			a.t=c->t[j];
			a.e=c->e[j];
			a.chain=c->myId;
			a.chainhit=j;
			a.view=3;
			

			Managed::ManagedCluster * clu = cluster_manager->GetCluster(c->cluster_id[j]);
			if(clu)
			{
				a.rms_t=clu->rms_t;
				spoints.push_back(a);
			}
		}
	
	MSG("PrimaryShowerFinder",Msg::kDebug)<<"MAKING PARTICLE WITH "<<spoints.size()<<" POINTS\n";
	
	//we don't want the particle to extrapolate too far.... ...but now go all the way to the end
	double stopend = endu<endv?endv:endu;

		
	sort(spoints.begin(), spoints.end(),spointgreaterlmf);
	
	for(unsigned int i=0;i<spoints.size();i++)
	{
//		if( start - spoints[i].z > 0.045 )continue;
//		if( spoints[i].z - end > 0.045)continue; //within 1 planes, we want the end of the chain
	
//	printf("spoint %f %f %f %d\n",spoints[i].z,spoints[i].t,spoints[i].e,spoints[i].view);
	
		if(spoints[i].z>stopend)continue;
	
		int myview = spoints[i].view;
		int lower=-1;
		int upper=-1;
		for(int j=i-1;j>-1;j--)
		{
			if(spoints[j].view!=myview)
			{
				lower=j;
				break;
			}
		}
		for(unsigned int j=i+1;j<spoints.size();j++)
		{
			if(spoints[j].view!=myview)
			{
				upper=j;
				break;
			}
		}	
		
		
		
		double u = spoints[i].view==2 ? spoints[i].t : 0;
		double v = spoints[i].view==3 ? spoints[i].t : 0;
		double e = spoints[i].e;
		double z = spoints[i].z;
		int view = spoints[i].view;
		
		double rms_t = spoints[i].rms_t;
		

		
		if(lower>-1 && upper > -1 )// good we can extrapolate!
		{
			double s = (spoints[upper].t - spoints[lower].t) /  ( spoints[upper].z - spoints[lower].z);
			
			double t = s * (spoints[i].z-spoints[lower].z) + spoints[lower].t;
			
			u = myview == 2 ? u : t;
			v = myview == 3 ? v : t;
			

		}else if(lower>-1 && upper < 0)  //just user the closest other view t value
		{
		

			//do we have another lower value?
			int lower2=-1;
			for(int j=lower-1;j>-1;j--)
			{
				if(spoints[j].view!=myview)
				{
					lower2=j;
					break;
				}
			}
			
			if(lower2>-1 && fabs( spoints[lower].z - spoints[lower2].z) >0)
			{
				double s = (spoints[lower].t - spoints[lower2].t) /  ( spoints[lower].z - spoints[lower2].z);
			
				double t = s * (spoints[i].z-spoints[lower2].z) + spoints[lower2].t;
				u = myview == 2 ? u : t;
				v = myview == 3 ? v : t;
			
			}else{
				u = myview == 2 ? u : spoints[lower].t;
				v = myview == 3 ? v : spoints[lower].t;
			}
		}
		else if(upper>-1 && lower < 0)   //just user the closest other view t value
		{
		

		
			//do we have another upper value?
			int upper2=-1;
			for(unsigned int j=upper+1;j<spoints.size();j++)
			{
				if(spoints[j].view!=myview)
				{
					upper2=j;
					break;
				}
			}
			
			if(upper2>-1 && fabs( spoints[upper2].z - spoints[upper].z)>0)
			{
				double s = (spoints[upper2].t - spoints[upper].t) /  ( spoints[upper2].z - spoints[upper].z);
			
				double t = s * (spoints[i].z-spoints[upper].z) + spoints[upper].t;
				u = myview == 2 ? u : t;
				v = myview == 3 ? v : t;
			
			}else{
				u = myview == 2 ? u : spoints[upper].t;
				v = myview == 3 ? v : spoints[upper].t;
			}



			//lets use the vertex!

/*
			if(spoints[upper].view != spoints[i].view)
				upper=upper2;
				
			if(spoints[upper].view == spoints[i].view)
			{
*/
//dont yet have a vertex!
/*				double vz =vtx_z;
				double vt = 0;
				vt = spoints[upper].view == 2 ? vtx_u : vt;
				vt = spoints[upper].view == 3 ? vtx_v : vt;
			
				double t =vt;
				
				//if the spoint at z is not the same as z vtx, we need to extrapolate 
				if(fabs ( spoints[upper].z - vz)>0.001)
				{			
			
					double s = (spoints[upper].t - vt) /  ( spoints[upper].z - vz);
					t = s * (spoints[i].z-vz) + vt;			
			
				}
				
		//		printf("---view %d u %f v %f t %f\n",myview,u,v,t);
		//		printf("---vtx z %f t%f upper z %f t %f\n",vz,vt,spoints[upper].z,spoints[upper].t);
			

				u = myview == 2 ? u : t;
				v = myview == 3 ? v : t;
*/
				u = myview == 2 ? u : spoints[upper].t;
				v = myview == 3 ? v : spoints[upper].t;

//			}	
				

		}
		else if(upper==-1 && lower==-1) //we have an empty view!!!
		{

			u = myview == 2 ? u : 0;
			v = myview == 3 ? v : 0;
		}
		
		foundparticle3d->add_to_back(u,v,z,e,spoints[i].chain,spoints[i].chainhit,view,rms_t);
		
	//	printf("adding 3d spoint to particle %f %f %f --- %f\n",u,v,z,e);
	}
	
	
	if(type)
		foundparticle3d->particletype=(Particle3D::EParticle3DType)type;
		
	foundparticle3d->finalize();
	
	if(foundparticle3d->entries<1){delete foundparticle3d;foundparticle3d=0;}


}



void PrimaryShowerFinder::LoadCloseHits(HoughLine * hl, Managed::ClusterManager *cm, int view, double max_tdist,int limit_to_current_size)
{
	if(!hl || !cm)return;

	std::map<double, std::map<double, int> > * cluster_map = cm->GetClusterMap(view);
    std::map<double, std::map<double, int> >::reverse_iterator p_iterr;
    std::map<double, int >::iterator s_iterr;


//printf("\nloading hit\n");
	double last_z=0;
	Managed::ManagedCluster *closest_cluster=0;
	double dist=1000;
    for(p_iterr=cluster_map->rbegin();p_iterr!=cluster_map->rend(); p_iterr++)
    {	
    	if(!last_z)last_z=p_iterr->first;
    	 	
  //  	printf("plane %f \n",p_iterr->first);
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
		{
				Managed::ManagedCluster *cl = cm->GetCluster(s_iterr->second);
				if(!cl)continue;
				if(cl->e<0.001)continue;
				
				//do we need to make sure that this cluster is contained withing the dimensions of the houghline?
				if(limit_to_current_size)
				{
					double min_z = hl->start_z < hl->end_z?hl->start_z:hl->end_z;
					double max_z = hl->start_z > hl->end_z?hl->start_z:hl->end_z;
					if(min_z-cl->z>0.07 || cl->z-max_z>0.07)
					{
			//			printf("not adding %f %f to hl with endpts %f %f  %f %f\n",cl->z,cl->t,hl->start_z,hl->start_t,hl->end_z,hl->end_t);
						continue;
					}
					
					//if this hl is contrained to a single plane... we need to keep hits within t!
					double min_t =  hl->start_t < hl->end_t ? hl->start_t:hl->end_t;
					double max_t =  hl->start_t > hl->end_t ? hl->start_t:hl->end_t;
					
					if(max_z-min_z<0.001 && (cl->t-min_t < 0.001 || cl->t - max_t > 0.001))continue;
					
					if(fabs(hl->phi-3.141592/2 )<0.1)
					{

						
						if(cl->t<min_t || cl->t>max_t)continue;
					}
				
				}				
				
				double exp_t = hl->GetExpectedT(cl->z);
//				printf(" cluster z %f\n",cl->z);
				//is it a vertical line?
				if(fabs(hl->phi-3.141592/2 )<0.1){
				
					//find where t=0
					double zpos = hl->offset_z + (hl->offset_t + hl->r/sin(hl->theta))*sin(hl->theta)/cos(hl->theta);
			//	printf("verticle match at %f against %f\n",zpos,cl->z);	
					//give more for the verticle hits	
					double tmt = max_tdist > 0.1 ? max_tdist:0.1;
					if(fabs(zpos-cl->z)<tmt)
					{
						hl->AddCluster(cl);
					}
				
				}else{
					//is it the closest?
					if(fabs(exp_t-cl->t)<max_tdist && fabs(exp_t-cl->t)<dist)
					{
						closest_cluster=cl;
						dist=fabs(exp_t-cl->t);
					}
				}	
					
		}

		std::map<double, std::map<double, int> >::reverse_iterator next_p_iterr=p_iterr;
		next_p_iterr++;
    	if(next_p_iterr!=cluster_map->rend())
    	if(fabs(p_iterr->first-next_p_iterr->first)>0.04 && closest_cluster)
    	{//new plane.. save the old plane
			hl->AddCluster(closest_cluster);
//			printf("next plane %f from %f\n",p_iterr->first,last_z);
  //  		printf("adding %f %f %f\n",closest_cluster->z,closest_cluster->t,closest_cluster->e);
    		closest_cluster=0;
    		dist=1000;
    		
    		last_z=p_iterr->first;
    	}


	}

	if(closest_cluster)
   	{//new plane.. save the old plane
 //  	printf("last_plane %f\n",last_z);
 //  			printf("adding %f %f %f\n",closest_cluster->z,closest_cluster->t,closest_cluster->e);
		hl->AddCluster(closest_cluster);
	}
}


void PrimaryShowerFinder::MakeHoughMap(int view)
{
	if(view !=2 && view !=3)return;
	std::map<double, std::map<double, int>  >::iterator p_iterr;
    std::map<double, int>::iterator s_iterr;
	
//	printf("----view %d----\n\n",view);
	
	if(houghmapU && view==2){houghmapU->Reset();}//delete houghmapU;houghmapU=0;}
	if(houghmapV && view==3){houghmapV->Reset();}//delete houghmapV;houghmapV=0;}
	
//	if(houghmapU==0)houghmapU= new TH2D("hhu","Hough U", 300, 0,3.141592,150,-1,1);
//	if(houghmapV==0)houghmapV= new TH2D("hhv","Hough V", 300, 0,3.141592,150,-1,1);

	if(houghmapU==0)houghmapU= new TH2D("hhu","Hough U", 100, 0,3.141592,50,-1,1);
	if(houghmapV==0)houghmapV= new TH2D("hhv","Hough V", 100, 0,3.141592,50,-1,1);


		
	TH2D * h = view==2?houghmapU:view==3?houghmapV:0;
	if(!h)return;

	//is the map empty?
	if(cluster_manager->GetClusterMap(view)->begin()==cluster_manager->GetClusterMap(view)->end())return; 

 //being on the first hit can caus some trouble... so move a bit away from it
	double off_z = cluster_manager->GetClusterMap(view)->begin()->first-0.1;
	double off_t = cluster_manager->GetClusterMap(view)->begin()->second.begin()->first-0.1;
	//Managed::ManagedCluster *largest=0;
  	//double last_plane=0;
  
  
   	TH2D copy;
	h->Copy(copy);    

	//find max value!
	double maxvale=0;
    for(p_iterr=cluster_manager->GetClusterMap(view)->begin();p_iterr!=cluster_manager->GetClusterMap(view)->end(); p_iterr++)
    {	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluster_manager->GetCluster(s_iterr->second);
			if(clus->e>maxvale)
			{
				maxvale=clus->e;
			}
		}
	}


    for(p_iterr=cluster_manager->GetClusterMap(view)->begin();p_iterr!=cluster_manager->GetClusterMap(view)->end(); p_iterr++)
    {	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluster_manager->GetCluster(s_iterr->second);
//printf("adding with e %f\n",clus->e);
			//don't care about low e hits...
			if(0.2 > clus->e)continue;
			for(double i=0;i<3.141592;i+=3.141592/(h->GetNbinsX()-1))
			{
				//double val = clus->e;
				double val=1;
				h->Fill(i, (clus->z-off_z) * cos(i) + (clus->t-off_t) * sin(i),  val);  
				copy.Fill(i, (clus->z-off_z) * cos(i) + (clus->t-off_t) * sin(i), clus->e);
			}			
	  	}
	}
	

   TH2D ch;
   h->Copy(ch);
	

	h->Reset();
	
	//int nx=ch.GetNbinsX();
	//int ny=ch.GetNbinsY();   

	while(ch.GetMaximum()>1)
	{
		int x;
		int y;
		int z;
		int m = ch.GetMaximumBin();
		ch.GetBinXYZ(m,x,y,z);

		SaveHitGroup(&ch,(TH2D*)h, (double)ch.GetBinContent(x,y),(double)ch.GetBinContent(x,y), x,y);
    }    


	multimap<double, std::pair<double,double> > top_map;
	int tops=0;	


///////////////////
//make a list of all peaks....to speed up the TH1::GetMaximum(Bin) which is painfully slow....

	std::vector<int> peakbinx;
	std::vector<int> peakbiny;
	std::vector<double> peakx;
	std::vector<double> peaky;
	std::vector<double> peakvalue;
	std::vector<int> peakstillgood;
	
	for(int i=1;i<h->GetNbinsX()+1;i++)
	for(int j=1;j<h->GetNbinsY()+1;j++)
	{
		double thisc = h->GetBinContent(i,j);
		if(thisc<0.0001)continue;
		
		if(
			thisc >= h->GetBinContent(i+1,j)
			&&
			thisc >= h->GetBinContent(i,j+1)
			&&
			thisc >= h->GetBinContent(i-1,j)
			&&
			thisc >= h->GetBinContent(i,j-1)
			&&
			thisc >= h->GetBinContent(i+1,j+1)
			&&
			thisc >= h->GetBinContent(i+1,j-1)
			&&
			thisc >= h->GetBinContent(i-1,j+1)
			&&
			thisc >= h->GetBinContent(i-1,j-1)
		)
		{
			peakx.push_back(h->GetXaxis()->GetBinCenter(i));
			peaky.push_back(h->GetYaxis()->GetBinCenter(j));
			peakbinx.push_back(i);
			peakbiny.push_back(j);
			peakvalue.push_back(thisc);	
			peakstillgood.push_back(1);
		}
	
	}

/*
	int npeaks = peakstillgood.size();
	
	printf("found %d peaks",npeaks);
	
	while(npeaks>0)
	{
	
		double xv=0;
		double yv=0;
		double val=0;
		int cnt =0;
			
		GetPeakAreaAverage(xv, yv,val, peakbinx,peakbiny,peakx,peaky,peakvalue,peakstillgood);


		top_map.insert(make_pair((double)val, std::pair<double,double>(xv,yv)) );
		tops++;

		npeaks=0;
		for(unsigned int i=0;i<peakstillgood.size();i++)	
			npeaks+=peakstillgood[i]==1?1:0;
	}
	
	printf("%d peaks made\n",tops);
*/
///////////////

//there must be a better way....

	int npeaks=0;
	for(unsigned int i=0;i<peakstillgood.size();i++)	
		npeaks+=peakstillgood[i]==1?1:0;
	
	
	
	while(npeaks>0)//h->GetMaximum()>0)
	{
		double xv=0;
		double yv=0;
		double val=0;
		int cnt =0;
		
		int x;
		int y;
		//int z;
//		int m = h->GetMaximumBin();
//		h->GetBinXYZ(m,x,y,z);

		double max=0;
		int maxi=-1;
		for(int i=0;i<(int)peakvalue.size();i++)
		{
			if(peakvalue[i]>max && peakstillgood[i])
			{
				max = peakvalue[i];
				maxi=i;
			}
		}

		if(maxi<0)break;
		
		x=peakbinx[maxi];
		y=peakbiny[maxi];
		
		GetPeakAreaAverage(xv, yv,val, cnt, x, y, (TH2D*)h, peakvalue, peakstillgood, peakbinx, peakbiny);
		xv/=(double)cnt;
		yv/=(double)cnt;
		
		val=copy.GetBinContent(copy.FindBin(xv,yv));
		
		top_map.insert(make_pair((double)val, std::pair<double,double>(xv,yv)) );
		tops++;
		
		npeaks=0;
		for(unsigned int i=0;i<peakstillgood.size();i++)	
			npeaks+=peakstillgood[i]==1?1:0;
	
	}




	
	
//	printf("printing top %d\n",tops);
	
	multimap<double, std::pair<double,double> >::reverse_iterator it = top_map.rbegin();

	
	std::vector<HoughLine> *houghlines = view==2?&houghlinesU:&houghlinesV;
	
	for(int k=0;k<tops;k++)
	{
		double theta=it->second.first;
		double r = it->second.second;	
		HoughLine hl(theta, r, off_t, off_z);
		
	//	printf("peak at theta %f r %f off_t %f off_z %f\n",theta,r,off_t,off_z);
	
		LoadCloseHits(&hl,cluster_manager,view,0.0412);
			
		if(hl.ncluster>1)	
			houghlines->push_back(hl);

		it++;

		
	}
	
	
	sort(&(*houghlines)[0],&(*houghlines)[houghlines->size()],CompareLength);
	
	//should remove duplicates --- (caused by the same clusters being matched to different hough map peaks)
	std::vector<HoughLine> houghlinestemp;

        houghlinestemp = *houghlines;
	houghlines->clear();

	for(int i=0;i<(int)houghlinestemp.size();i++)
	{
		int dup=0;
		HoughLine * l = &houghlinestemp[i];
		for(int j=i+1;j<(int)houghlinestemp.size();j++)
		{
			HoughLine *r = &houghlinestemp[j];
			if(l && r && l->start_z == r->start_z && l->start_t == r->start_t && l ->end_z == r->end_z && l->end_t == r->end_t)
			{
				dup=1;
				break;
			}

		}

		if(!dup && houghlinestemp[i].ncluster>1)houghlines->push_back(houghlinestemp[i]);

	}


/* //this code is highly flawed...
	for(int i=0;i<houghlines->size();i++)
	{
		HoughLine * l = &(*houghlines)[i];
		HoughLine * r=0;
		if(i<houghlines->size()-1)r = &(*houghlines)[i+1];
		if((*houghlines)[i].ncluster>1)houghlinestemp.push_back((*houghlines)[i]);		
		if(l && r && l->start_z == r->start_z && l->start_t == r->start_t && l ->end_z == r->end_z && l->end_t == r->end_t)
		{

			while (i<houghlines->size())
			{
				i++;
				HoughLine * l = &(*houghlines)[i];
				HoughLine * r = &(*houghlines)[i+1];
				if(! (l->start_z == r->start_z && l->start_t == r->start_t && l ->end_z == r->end_z && l->end_t == r->end_t)) break;
			}
		}
	}

*/

	

	
	sort(&(*houghlines)[0],&(*houghlines)[houghlines->size()],CompareForwardAndClusters);
	
	//int max =-1;
	//if(houghlines->size()>10)max = houghlines->size()-10;
	for(int i=0;i<(int)houghlines->size();i++)
	{
		if((*houghlines)[i].ncluster<1)break;
//		printf("%d houghline with e %f clusters %d (t,z) (%f,%f) (%f,%f) chi2/ndf %f theta %f r %f phiz %f \n", i,(*houghlines)[i].sum_e, (*houghlines)[i].ncluster, (*houghlines)[i].start_t, (*houghlines)[i].start_z, (*houghlines)[i].end_t, (*houghlines)[i].end_z, (*houghlines)[i].chi2/(*houghlines)[i].ncluster,(*houghlines)[i].theta,(*houghlines)[i].r,(*houghlines)[i].phi);
		
									
	}
	
	

}


 
void PrimaryShowerFinder::SaveHitGroup(TH2D * his, TH2D * saver, double save_val, double with_val, int curx, int cury)
{
	if(curx<1 || cury < 1 || curx> his->GetNbinsX() || cury >his->GetNbinsY())return;
	

	
	double thisval = his->GetBinContent(curx,cury);

//	printf("saving %d %d thisval %f saveval %f withval %f\n",curx,cury,thisval,save_val, with_val);
	if(thisval==0)return;
	if(fabs(thisval-save_val)<0.001)
	{
		saver->SetBinContent(curx, cury,thisval);
		//thisval=0;
	}
	else if(thisval>=with_val)return;


	his->SetBinContent(curx,cury,0);
			
	SaveHitGroup(his,saver,save_val,thisval,curx+1,cury);
	SaveHitGroup(his,saver,save_val,thisval,curx,cury+1);
	SaveHitGroup(his,saver,save_val,thisval,curx-1,cury);
	SaveHitGroup(his,saver,save_val,thisval,curx,cury-1);
	SaveHitGroup(his,saver,save_val,thisval,curx+1,cury+1);
	SaveHitGroup(his,saver,save_val,thisval,curx-1,cury-1);
	SaveHitGroup(his,saver,save_val,thisval,curx+1,cury-1);
	SaveHitGroup(his,saver,save_val,thisval,curx-1,cury+1);
//	printf("save recursion done\n");
	

	//printf("saving points %d %d %f\n",curx,cury,thisval);

}

/*
void PrimaryShowerFinder::GetPeakAreaAverage(double &xv, double &yv, double &val, std::vector<int> &peakbinx, std::vector<int> &peakbiny, std::vector<double> &peakx, std::vector<double> &peaky, std::vector<double> &peakvalue, std::vector<int> &peakstillgood)
{
	xv=0;
	yv=0;
	val=0;
	
	//goal is to find all ajacent bins and group them into one bin....
	
	int peak=-1;
	
	//find the first available peak
	for(unsigned int i=0;i<peakstillgood.size();i++)
	{
		if(peakstillgood[i]==1)
		{
			peak=i;
			break;
		}
	}
	
	if(peak<0)return;
	
	val=peakvalue[peak];
	
	int binx=peakbinx[peak];
	int biny=peakbiny[peak];
	
	int lastaddcount=1;
	
	std::vector<int> addedx;
	std::vector<int> addedy;
	
	addedx.push_back(binx);
	addedy.push_back(biny);
	peakstillgood[peak]=0;
	
	while(lastaddcount>0)
	{
		lastaddcount=0;
		for(unsigned int i=0;i<peakstillgood.size();i++)
		{
			if(peakstillgood[i])
			{
				for(unsigned int j=0;j<addedx.size();j++)
				{
				
					if(peakbinx[i] == peakbinx[addedx[j]]+1  && peakbiny[i] == peakbiny[addedy[j]]
						||
						peakbinx[i] == peakbinx[addedx[j]]  && peakbiny[i] == peakbiny[addedy[j]]+1
						||
						peakbinx[i] == peakbinx[addedx[j]]-1  && peakbiny[i] == peakbiny[addedy[j]]
						||
						peakbinx[i] == peakbinx[addedx[j]]  && peakbiny[i] == peakbiny[addedy[j]]-1
						||
						peakbinx[i] == peakbinx[addedx[j]]+1  && peakbiny[i] == peakbiny[addedy[j]]+1
						||
						peakbinx[i] == peakbinx[addedx[j]]+1  && peakbiny[i] == peakbiny[addedy[j]]-1
						||
						peakbinx[i] == peakbinx[addedx[j]]-1  && peakbiny[i] == peakbiny[addedy[j]]+1
						||
						peakbinx[i] == peakbinx[addedx[j]]-1  && peakbiny[i] == peakbiny[addedy[j]]-1
					){
						peakstillgood[i]=0;
						addedx.push_back(peakbinx[i]);
						addedy.push_back(peakbiny[i]);
						lastaddcount++;
						break;
					}
					
					
					
			
				}
			}
		
		}	
	}
	
	
	//calculate the value
	
	double cnt=0;
	for(unsigned int i=0;i<addedx.size();i++)
	{
		xv+=peakx[addedx[i]];
		yv+=peaky[addedy[i]];
		cnt+=1;
	}
	
	xv/=cnt;
	yv/=cnt;
	

}*/


//assuming that the region is bounded by bins with 0 as their content
void PrimaryShowerFinder::GetPeakAreaAverage(double &x, double &y,double &val, int & cnt, int curx, int cury, TH2D * hist, std::vector<double> & peakvalue, std::vector<int> &peakstillgood, std::vector<int> &peakbinx,std::vector<int> & peakbiny)
{

	double thisval = hist->GetBinContent(curx,cury);

	if(thisval==0)return;
	
	val+=thisval;
	x+=hist->GetXaxis()->GetBinCenter(curx);
	y+=hist->GetYaxis()->GetBinCenter(cury);
	cnt++;
	hist->SetBinContent(curx,cury,0);

	/////
	
	for(int i=0;i<(int)peakbinx.size();i++)
	{
		if(peakbinx[i]==curx && peakbiny[i]==cury)
		{
			peakstillgood[i]=0;
			peakvalue[i]=0;
			break;
		}
	}
	
	
	/////


	
	GetPeakAreaAverage(x, y, val,cnt, curx+1, cury, hist, peakvalue, peakstillgood, peakbinx, peakbiny);
	GetPeakAreaAverage(x, y, val,cnt, curx, cury+1, hist, peakvalue, peakstillgood, peakbinx, peakbiny);
	GetPeakAreaAverage(x, y, val,cnt, curx, cury-1, hist, peakvalue, peakstillgood, peakbinx, peakbiny);
	GetPeakAreaAverage(x, y, val,cnt, curx-1, cury, hist, peakvalue, peakstillgood, peakbinx, peakbiny);
	GetPeakAreaAverage(x, y, val,cnt, curx+1, cury+1, hist, peakvalue, peakstillgood, peakbinx, peakbiny);
	GetPeakAreaAverage(x, y, val,cnt, curx-1, cury-1, hist, peakvalue, peakstillgood, peakbinx, peakbiny);
	GetPeakAreaAverage(x, y, val,cnt, curx+1, cury-1, hist, peakvalue, peakstillgood, peakbinx, peakbiny);
	GetPeakAreaAverage(x, y, val,cnt, curx-1, cury+1, hist, peakvalue, peakstillgood, peakbinx, peakbiny);

		
}

void PrimaryShowerFinder::FindHoughLineMatches()
{
	houghlineMatch.clear();

	//need to iterate over the hough lines in both view and look for matches

	for(unsigned int i=0;i<houghlinesU.size();i++)
	{
		//don't want to match to verticle lines
		if(fabs(houghlinesU[i].phi -3.141592/2)<0.1)continue;
		if(houghlinesU[i].ncluster<2)continue;
			HoughLine * l = &houghlinesU[i];

		for(unsigned int j=0;j<houghlinesV.size();j++)
		{
			//don't want to match to verticle lines
			if(fabs(houghlinesV[j].phi - 3.141592/2)<0.1)continue;
			if(houghlinesV[j].ncluster<2)continue;

			HoughLine * r = &houghlinesV[j];
			
			double overlap = l->end_z<r->end_z?l->end_z:r->end_z;
			overlap -= l->start_z > r->start_z ? l->start_z:r->start_z;
			
			//require overlap to exceed 90% of the length of the smaller line
			double l_z = l->end_z-l->start_z;
			double r_z = r->end_z-r->start_z;
			double smaller_z = l_z < r_z?l_z:r_z;
			if(smaller_z * 0.5 > overlap)
			{
			//	printf("failing %d %d with overlap %f > %f\n",i,j,smaller_z * 0.5,overlap);
				continue;
			}
			double smaller_e = l->sum_e < r->sum_e ? l->sum_e:r->sum_e;
			double larger_e = l->sum_e > r->sum_e ? l->sum_e:r->sum_e;

			//require smaller e to be more than 30% of larger e
			if(larger_e*0.3 > smaller_e)
			{
		//		printf("failing %d %d with energy %f > %f\n",i,j,larger_e*0.5 , smaller_e);
				continue;
			}


			//require a peak hit that is at least 20% of total e
		

			
			double max_hit_e=0;
			double sum_hit_e=0;
			for(unsigned int ka=0;ka<l->cluster_id.size();ka++)
			{
				Managed::ManagedCluster * clus = cluster_manager->GetCluster(l->cluster_id[ka]);
				if(!clus){
					//printf("missing cluster %d\n",l->cluster_id[ka]);
					continue;
				}
				sum_hit_e+=clus->e;
				
				if(clus->e>max_hit_e)max_hit_e=clus->e;
			
			}
			for(unsigned int ka=0;ka<r->cluster_id.size();ka++)
			{
				Managed::ManagedCluster * clus = cluster_manager->GetCluster(r->cluster_id[ka]);
				if(!clus){
					//printf("missing cluster %d\n",r->cluster_id[ka]);
					continue;
				}					
				sum_hit_e+=clus->e;
				
				if(clus->e>max_hit_e)max_hit_e=clus->e;
			
			}	
			if(sum_hit_e*0.2 > max_hit_e)
			{
			//	printf("continuing sum hit*.3 %f > maxhit %f\n",sum_hit_e*0.3, max_hit_e);
				continue;
			}
			if( max_hit_e<5)
			{
			//	printf("maxhit too small at %f\n", max_hit_e);
				continue;
			}
			
		//	printf("saving %d %d\n",i,j);
			houghlineMatch.push_back(make_pair(i,j));					
		}
	}

	std::map<double,int> orderer;
	for(unsigned int i=0;i<houghlineMatch.size();i++)
	{
		HoughLine * l = &houghlinesU[houghlineMatch[i].first];
		HoughLine * r = &houghlinesV[houghlineMatch[i].second];
		if(l->ncluster<2 || r->ncluster<2)continue;
		double weight = (l->sum_e+r->sum_e)*(cos(l->phi))*(cos(r->phi));
		weight *= l->ncluster / sqrt( (l->end_z-l->start_z)*(l->end_z-l->start_z)+(l->end_t-l->start_t)*(l->end_t-l->start_t));
		weight *= r->ncluster / sqrt( (r->end_z-r->start_z)*(r->end_z-r->start_z)+(r->end_t-r->start_t)*(r->end_t-r->start_t));

		orderer.insert(make_pair(weight,i));
	}		

	std::map<double,int>::iterator it_orderer;

	
	//printf("found %d preliminary matches\n",houghlineMatch.size());

	/*for(it_orderer = orderer.begin();it_orderer!=orderer.end();it_orderer++)
	//for(unsigned int i=0;i<houghlineMatch.size();i++)
	{
		int i= it_orderer->second;
		
		//printf("match %d %d\n",houghlineMatch[i].first,houghlineMatch[i].second);
		//HoughLine * l = &houghlinesU[houghlineMatch[i].first];
		//HoughLine * r = &houghlinesV[houghlineMatch[i].second];
	
		//printf("\t U s %f %f e %f %f sume %f clust %d phi %f\n",l->start_z,l->start_t, l->end_z,l->end_t,l->sum_e,l->ncluster,l->phi);
		//printf("\t V s %f %f e %f %f sume %f clust %d phi %f\n",r->start_z,r->start_t, r->end_z,r->end_t,r->sum_e,r->ncluster,r->phi);
	}
*/
}

//make a chain using the primary shower
Chain *  PrimaryShowerFinder::MakeShowerChain(int view)
{
	if(houghlineMatch.size()<1)return 0;
	
	Chain * c = new Chain();

	//take the first match as the best for now...
	HoughLine * l = 0;
	/*if(view ==2)l=&houghlinesU[houghlineMatch[houghlineMatch.size()-1].first];
	else if(view ==3)l=&houghlinesV[houghlineMatch[houghlineMatch.size()-1].second];
	else return 0;
	*/

/*	
	if(view ==2)l=&houghlinesU[houghlinesU.size()-1];
	else if(view ==3)l=&houghlinesV[houghlinesV.size()-1];
	else return 0;
*/	

	if(view ==2)l=&houghlinesU[houghlineMatch[houghlineMatch.size()-1].first];
	else if(view ==3)l=&houghlinesV[houghlineMatch[houghlineMatch.size()-1].second];
	else return 0;
	
	
	l->primary=1;
	
	//printf("making chain for view %d\n",view);
	//printf("look at houghline %d\n",view==2?houghlineMatch[houghlineMatch.size()-1].first:houghlineMatch[houghlineMatch.size()-1].second);
	for(int i=0;i<l->ncluster;i++)
	{
	//	printf("getting cluster %d\n",l->cluster_id[i]);
		Managed::ManagedCluster *clus = cluster_manager->GetCluster(l->cluster_id[i]);
		if(!clus)continue;
		c->insert(clus->t,clus->z,clus->e,clus->id);
	//	printf("inserting %f %f %f %d\n",clus->t,clus->z,clus->e,clus->id);	
	}


	return c;
}


void PrimaryShowerFinder::DumpHoughLines(int view)
{
	
	if(!MsgService::Instance()->IsActive("PrimaryShowerFinder",Msg::kDebug))return;
	
	printf("\nView %d\n",view);
	std::vector<HoughLine>* houghlines=0;
	if(view==2)houghlines=&houghlinesU;
	if(view==3)houghlines=&houghlinesV;
	if(!houghlines)return;

	for(int i=(*houghlines).size()-1;i>-1;i--)
	{
		if((*houghlines)[i].ncluster<1)continue;
		printf("houghline with e %f clusters %d (t,z) (%f,%f) (%f,%f) chi2/ndf %f  offz %f expt %f phi %f\n", (*houghlines)[i].sum_e, (*houghlines)[i].ncluster, (*houghlines)[i].start_t, (*houghlines)[i].start_z, (*houghlines)[i].end_t, (*houghlines)[i].end_z, (*houghlines)[i].chi2/(*houghlines)[i].ncluster,(*houghlines)[i].offset_z,(*houghlines)[i].GetExpectedT((*houghlines)[i].offset_z),(*houghlines)[i].phi);
	}

	printf("\n");
}
