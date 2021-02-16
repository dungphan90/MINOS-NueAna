#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterManager.h"
#include "MessageService/MsgService.h"

#include <sstream>

#include <math.h>
using namespace Managed;

ClassImp(ClusterManager)

CVSID("$Id: ClusterManager.cxx,v 1.1 2009/06/19 14:32:38 scavan Exp $");

ClusterManager::ClusterManager()
{

		clusters_are_made=0;
		
		loc_map.clear();
		cluster_map.clear();
		cluster_map_u.clear();
		cluster_map_v.clear();
		hits.clear();
		clusters.clear();	
		
		clusters_to_delete.clear();
		
		minu=100000;
		maxu=-100000;
		minv=100000;
		maxv=-100000;
		minz=100000;
		maxz=-100000;
		clustersaver=0;
		
		hitmanager=0;
		
		needMapRebuild=0;

}


//return a vector of cluster ids 
//for clusters with a z position within 1/2 plane of z
std::vector<int> ClusterManager::FindClustersInZ(double z,int view=0)
{
	double planedist = 0.04;
	
	
	std::map<double, std::map<double, int>  >*cmap = 0;
	if(view==2)cmap=&cluster_map_u;
	else if(view==3)cmap=&cluster_map_v;
	else cmap=&cluster_map;
	
	std::map<double, std::map<double, int>  >::iterator it = cmap->lower_bound(z-planedist);
	
	std::vector<int> rec;
	for(;it!=cmap->end();it++)
	{
		if(it->first > z+planedist)break;
		std::map<double, int>::iterator s_iter;
		for(s_iter=it->second.begin();s_iter!=it->second.end(); s_iter++)
		{
			rec.push_back(s_iter->second);
		}	
	
	}

	return rec;
}


std::map<double, std::map<double, int>  > * ClusterManager::GetClusterMap(int view)
{

	//need rebuild?
	
	if(needMapRebuild)RebuildClusterMaps();



	if(view==2)return &cluster_map_u;
	if(view==3)return &cluster_map_v;
	return &cluster_map;
}



void ClusterManager::RebuildClusterMaps()
{


	if(clusters_to_delete.size()>0)
	{
		std::vector<Managed::ManagedCluster> cluster_temp;
		for(unsigned int i=0;i<clusters.size();i++)
		{
			int keep=1;
			for(unsigned int j=0;j<clusters_to_delete.size();j++)
			{
				if(clusters_to_delete[j]==clusters[i].id)
				{
					keep=0;
					break;
				}
			}
			if(!keep)continue;
			
			cluster_temp.push_back(clusters[i]);
		}
	
		clusters=cluster_temp;
		cluster_temp.clear();
		clusters_to_delete.clear();
	}

	cluster_map.clear();
	cluster_map_u.clear();
	cluster_map_v.clear();

	for(unsigned int i=0;i<clusters.size();i++)
	{
		cluster_map[clusters[i].z][clusters[i].t]=clusters[i].id;
		if(clusters[i].view==2)cluster_map_u[clusters[i].z][clusters[i].t]=clusters[i].id;
		if(clusters[i].view==3)cluster_map_v[clusters[i].z][clusters[i].t]=clusters[i].id;	   					   		
	}

	needMapRebuild=0;
}

//save a cluster to the cluster saver if it has been defined
//return the id of the new cluster so that objects can be updated....
int ClusterManager::SaveCluster(int cluster_id,double energy_to_save,int status)
{
	if(cluster_id<0)
	{
//		printf("!!!! attempt to save an already saved cluster!\n");
		return cluster_id;//its already a saved cluster!
	}
	//tell them we need a new map because we may be removing clusters....
	needMapRebuild=1;

	ManagedCluster *cluster = GetCluster(cluster_id);
	
	if(!cluster)return 0;
	
	/////////////////////////
	//allow the ability to fake a save cluster for tempory work
	if(!clustersaver)
	{
//		printf("SAVING WITHOUT SAVER!\n");
		//can only allocation the full energy in the cluster...
		if(energy_to_save-cluster->e >= 0)
		{
//			printf("setting e to 0, was %f\n",cluster->e);
			SplitEnergy(cluster,cluster->e,0);
			if(hitmanager)AdjustCluster(cluster);
			return 0;
		
		}	
		
			SplitEnergy(cluster,energy_to_save,0);
			if(hitmanager)AdjustCluster(cluster);
	
		return 0;
	}
	////////////////////////
	
	
	

	//can only allocation the full energy in the cluster...
	if(energy_to_save-cluster->e >= 0 || energy_to_save==0)
	{
		energy_to_save=cluster->e;
		int newid= clustersaver->SaveCluster(cluster);
		ManagedCluster *newcluster=GetCluster(newid);
		if(newcluster)
		{
			newcluster->SetStatus(status);
			SplitEnergy(cluster,newcluster->e,0);
		}else{
			SplitEnergy(cluster,0,0);//throw away a cluster since the new cluster can't be made... probably because the old cluster had 0 e anyways.
		}
		AdjustCluster(cluster);
		
	//	if(newcluster)printf("saving full cluster with e %f  old has %f\n",newcluster->e,cluster->e);
		return newid;
	}
			
			
	//if we are here we want to save only part of the cluster!
	ManagedCluster newcluster = *cluster;
			
	SplitEnergy(&newcluster,energy_to_save,1);
	SplitEnergy(cluster,energy_to_save,0);
			
	newcluster.SetStatus(status);
	int newid = clustersaver->SaveCluster(&newcluster);
			
	AdjustCluster(cluster);

//	printf("saving split save has %f old has %f\n",GetCluster(newid)->e,cluster->e);
			
	return newid;


}

//once we have adjusted the energy of a cluster (because we saved part of it into another cluster
//we need to adjust the cluster that is in the cluster manager... 
//   as well as modifying the hits in the hitsmanager
//pass in the cluster taken from this cluster manager, but with the hits already modified
void ClusterManager::AdjustCluster(Managed::ManagedCluster *cluster)
{
	if(!cluster)return;
	ManagedCluster * current_cluster = GetCluster(cluster->id);
	if(!current_cluster)return;
	
	MSG("ClusterManager",Msg::kDebug)<<"cluster e "<<cluster->e<<" curcluster "<<current_cluster->e<<"\n";
	
	//is the passed object a part of the current cluster vector, or was it copied?
	//if its copied, make sure it gets copied back into the version in the vector
//	if(fabs(cluster->e - current_cluster->e)<1e-9)
//		return;
	
	//now we need to update the hits!
	int modified_cluster=0;
	for(unsigned int i=0;i<cluster->hit_id.size();i++)
	{
		ManagedHit * hit = hitmanager->FindHit(cluster->hit_id[i]);
		
		double setenergy = 0;
		if(hit)setenergy = hit->SetEnergy(cluster->hite[i]);
		
		MSG("ClusterManager",Msg::kDebug)<<"setting energy to "<<cluster->hite[i]<<" and it set "<<setenergy<<"\n";
		if(fabs(setenergy-current_cluster->hite[i])>1e-9)
		{
			MSG("ClusterManager",Msg::kDebug)<<"modifying cluster\n";
			current_cluster->hite[i]=setenergy;
			modified_cluster=1;
		}
	}
		
	if(modified_cluster)current_cluster->Finalize();
	
	if(current_cluster->e<0.001)clusters_to_delete.push_back(current_cluster->id);

}

//modify the passed cluster to make it have the specified amount of energy if matchenergy=1
//if matched energy=0, take away the given energy 
void ClusterManager::SplitEnergy(Managed::ManagedCluster * cluster,double energy, int matchenergy)
{

	MSG("ClusterManager",Msg::kDebug)<<"splitting cluster with "<<cluster->e<<" to have "<<energy <<" match "<<matchenergy<<"\n";

	double totale=cluster->e;
	if(totale<1e-9)return;
	for(unsigned int i=0;i<cluster->hite.size();i++)
	{
	
		double newe =0;
		if(matchenergy)newe=cluster->hite[i]/totale*energy;
		else newe=cluster->hite[i]/totale*(totale-energy);
		
	/*	printf("hit with %f should have %f\n",cluster->hite[i],newe);
		
		if(matchenergy)
		{
			//needs to also adjust the hit stored in the hit manager here.....
			Managed::ManagedHit *hit = hitmanager->FindHit(cluster->hit_id[i]);		
			double sete = hit->SetEnergy(newe);
			cluster->hite[i]=sete;
			printf("setting to %f\n",sete);
		}else{

			//add save the rest as a new hit
			if(newe>0)
			{
				printf("newhit energy %f\n",newe);
				Managed::ManagedHit *hit = hitmanager->FindHit(cluster->hit_id[i]);	
				int newid = hitmanager->InsertHit(hit->GetView(), hit->GetPlane(), hit->GetStrip(), hit->GetZ(), hit->GetT(), newe);
				cluster->hite[i]=newe;
				cluster->hit_id[i]=newid;
				printf("setting to %f\n",newe);
			}else{
				cluster->hite[i]=0;
				cluster->hit_id[i]=-1;
				printf("setting to 0");
			}
		}
		
	*/
	cluster->hite[i]=newe;	

//		if(matchenergy)hite[i]=hite[i]/totale*energy;
//		else hite[i]=hite[i]/totale*(totale-energy);

	}
	cluster->Finalize();

}


//get a vector of indicies in the clusters vector corresponding to a given view
std::vector<int>  ClusterManager::GetViewIndex(int view)
{
	std::vector<int> ret;
	
	for(unsigned int i=0;i<clusters.size();i++)
		if(clusters[i].view==view)ret.push_back(i);
		
	return ret;

}

ClusterManager::~ClusterManager()
{}


void ClusterManager::DumpClusters()
{
	ostringstream os;

	std::map<double, std::map<double, int> >::iterator p_iter;
	std::map<double, int>::iterator s_iter;

	os<<"dumping clusters (all)...\n";
           
	for(p_iter=cluster_map.begin();p_iter!=cluster_map.end(); p_iter++)
	{
		std::map<double, int>::iterator s_iter;
		for(s_iter=p_iter->second.begin();s_iter!=p_iter->second.end(); s_iter++)
		{
			ManagedCluster *mc = GetCluster(s_iter->second);
			if(!mc){
				os<<"missing cluster id "<<s_iter->second<<"\n";
				continue;
			}
			os<<"id "<<s_iter->second<<"   z "<<mc->z<<" t "<<mc->t<<" e "<<mc->e<<" dz "<<mc->dz<<" dt "<<mc->dt<<" view "<<mc->view<<"\n";
		
	
		}
	}


	os<<"dumping clusters (u)...\n";
           
	for(p_iter=cluster_map_u.begin();p_iter!=cluster_map_u.end(); p_iter++)
	{
		std::map<double, int>::iterator s_iter;
		for(s_iter=p_iter->second.begin();s_iter!=p_iter->second.end(); s_iter++)
		{
			ManagedCluster *mc = GetCluster(s_iter->second);
			if(!mc){
				os<<"missing cluster id "<<s_iter->second<<"\n";
				continue;
			}
			os<<"z "<<mc->z<<" t "<<mc->t<<" e "<<mc->e<<" dz "<<mc->dz<<" dt "<<mc->dt<<" view "<<mc->view<<"\n";
			os<<"id "<<s_iter->second<<"\n";
	
		}
	}
	
	os<<"dumping clusters (v)...\n";
           
	for(p_iter=cluster_map_v.begin();p_iter!=cluster_map_v.end(); p_iter++)
	{
		std::map<double, int>::iterator s_iter;
		for(s_iter=p_iter->second.begin();s_iter!=p_iter->second.end(); s_iter++)
		{
			ManagedCluster *mc = GetCluster(s_iter->second);
			if(!mc){
				os<<"missing cluster id "<<s_iter->second<<"\n";
				continue;
			}
			os<<"z "<<mc->z<<" t "<<mc->t<<" e "<<mc->e<<" dz "<<mc->dz<<" dt "<<mc->dt<<" view "<<mc->view<<"\n";
			os<<"id "<<s_iter->second<<"\n";
	
		}
	}
	
	
	os<<"done....\n\n";
	
	
	MSG("ClusterManager",Msg::kInfo)<<os;
}




void ClusterManager::MakeClusters(double min_cluster_e, double min_strip_e, double max_t_skip)
{

	//see if we have hits
	if(hits.size()==0 && hitmanager==0)return;
	
	LoadHits();

	
	if(hits.size()==0)return;

	////////////////////////////////////////////////////////////
	///// Make initial clustering for vertex finding and shower location guess

//	double thresh=0.4; ///for a cluster
//	double sthresh=0.4;//for an individual strip   0.2 mips ~ 1 pe   
						// roughly gaussian 7.361 w/ rms 1.252 pe/mip - lowest bin is 5 pe/mip
	//int skip=1;	
//	double skipt=0.05; //strip width is 0.025?

	
	double thresh=min_cluster_e;
	double sthresh=min_strip_e;
	double skipt=max_t_skip;



	std::vector<ManagedCluster> theclusters;


	std::map<int, std::map<int, int> >::iterator p_iter;
	std::map<int, int>::iterator s_iter;
	double loc_z=-1;

	//int last_strip=-1;
	int numstrip=0;
	double sum_e_t=0;
	double sum_e=0;	
	double lastt=-100000.0;
	int thisview=0;
           
	ManagedCluster c;  
	c.ResetIDCounter();
	c.AdvanceID();
         
         
	// cout <<"idcounter 1 " <<c.idcounter<< " "<<c.id<<endl;
           
           
	for(p_iter=loc_map.begin();p_iter!=loc_map.end(); p_iter++)
	{
		std::map<int, int>::iterator s_iter;
		for(s_iter=p_iter->second.begin();s_iter!=p_iter->second.end(); s_iter++)
		{
			if(hits[s_iter->second].GetERemaining()<sthresh)continue;
			if(numstrip==0) lastt=hits[s_iter->second].GetT();

		 	if((fabs(hits[s_iter->second].GetT()-lastt)>skipt|| loc_z!=hits[s_iter->second].GetZ()) && numstrip>0)
			{  //we have finished with this cluster, so record it and reset info
				sum_e_t/=sum_e;  //weighted average t
			   	if(thresh<sum_e)
			   	{

			   		cluster_map[loc_z][sum_e_t]=c.id;
			   		if(thisview==2)cluster_map_u[loc_z][sum_e_t]=c.id;
			   		if(thisview==3)cluster_map_v[loc_z][sum_e_t]=c.id;
			   					   		
			   		c.Finalize();
					theclusters.push_back(c);

				}



				//printf("----cluster at %f %f %f %d   id %d\n\n", loc_z, sum_e_t,sum_e, thisview,c.id);

			   	sum_e_t=0;
			   	sum_e=0;
				//   loc_z=-1;
			   	numstrip=0;	
			   	thisview=0;
				c.Reset();
				c.AdvanceID();
				
				   // cout <<"idcounter 2 " <<c.idcounter<< " "<<c.id<<endl;
					
				if(loc_z!=hits[s_iter->second].GetZ())
				{
					lastt=-100000.0;
				
				}
				//last_strip=-1;
			}



			
            loc_z=hits[s_iter->second].GetZ();
			lastt=hits[s_iter->second].GetT();
			//last_strip=s_iter->first;
			numstrip++;
			sum_e+=hits[s_iter->second].GetERemaining();
			sum_e_t+=hits[s_iter->second].GetERemaining()*hits[s_iter->second].GetT();
			thisview=hits[s_iter->second].GetView();
			c.view=hits[s_iter->second].GetView();
			c.Insert(hits[s_iter->second].GetZ(),hits[s_iter->second].GetT(),hits[s_iter->second].GetERemaining(),hits[s_iter->second].GetPlane(),hits[s_iter->second].GetStrip(),hits[s_iter->second].GetID());



			//printf("%d - %d %d %f %f %f %d\n",numstrip,p_iter->first, s_iter->first, hits[s_iter->second].GetERemaining(),hits[s_iter->second].GetZ(),hits[s_iter->second].GetT(),hits[s_iter->second].GetView());
		 }
	}
	//get the last cluster///
	
	if(numstrip>0)
	{
		sum_e_t/=sum_e;  //weighted average t
        if(thresh<sum_e)
        {
        	cluster_map[loc_z][sum_e_t]=c.id;
			if(thisview==2)cluster_map_u[loc_z][sum_e_t]=c.id;
			if(thisview==3)cluster_map_v[loc_z][sum_e_t]=c.id;
		
        	c.Finalize();
			theclusters.push_back(c);
		}
	
	}
	MSG("ClusterManager",Msg::kDebug)<<"----cluster at "<< loc_z<<" "<<sum_e_t<<" "<<sum_e<<" "<<thisview<<" id "<<c.id<<"\n";
	
	/*
	//uncomment this to bypass the sharing

	for(unsigned int k=0;k<theclusters.size();k++)	
	{
		theclusters[k].Finalize();
		clusters.push_back(theclusters[k]);
	}

	return;
	*/
	////
	
	c.Reset();
	c.AdvanceID();
					
	//cout <<"idcounter 3 " <<c.idcounter<< " "<<c.id<<endl;

	    
	std::vector<ManagedCluster> tmponepeak;
	std::vector<ManagedCluster> tmpmultpeak;
	//go through clusters, looking for multiple peaks....
	for(unsigned int k=0;k<theclusters.size();k++)	
	{
		std::map<double,int>::iterator it;
		it=theclusters[k].tsortmap.begin();
	
		int npeak=0;

		

		std::vector<int> maxidx;
		std::vector<double> maxe;
		
		std::vector<int> sortedidx;
		std::vector<double> sortede;
		
		MSG("ClusterManager",Msg::kDebug)<<"clu idx "<< k<<" id "<<theclusters[k].id <<" with "<<(int)theclusters[k].hite.size() <<" hits \n";

		if(theclusters[k].hite.size()==1)
		{
			theclusters[k].Finalize();
	
			if(theclusters[k].e>0.001)
				tmponepeak.push_back(theclusters[k]);
			continue;

		
		}



		
		///////fixed code
		int plastidx=-1;
		double plaste=0;
		int lastidx=it->second;
		double laste = theclusters[k].hite[it->second];
					
		MSG("ClusterManager",Msg::kDebug)<<"\t z "<<theclusters[k].hitz[it->second]<<" t "<<theclusters[k].hitt[it->second]<<" e "<<theclusters[k].hite[it->second]<<"\n";

		it++;		
		for(;it!=theclusters[k].tsortmap.end();it++)
		{
		
			double e=theclusters[k].hite[it->second];
			
			MSG("ClusterManager",Msg::kDebug)<<"\t z "<<theclusters[k].hitz[it->second]<<" t "<<theclusters[k].hitt[it->second]<<" e "<<theclusters[k].hite[it->second]<<"\n";
			
			if(  e < laste&&(laste>plaste || plastidx<0) )
			{	
				maxidx.push_back(lastidx);
				maxe.push_back(laste);
				MSG("ClusterManager",Msg::kDebug)<<"\t\tpeak e "<<laste<<"\n";
				npeak++;			
			}else{
				//these hits needs to be shared!
				sortedidx.push_back(lastidx);
				sortede.push_back(laste);
			}
			plaste=laste;
			plastidx=lastidx;
			laste=e;
			lastidx=it->second;
					
		}		
			
		//last hit
		if(plastidx>-1 && laste>plaste)
		{
			maxidx.push_back(lastidx);
			maxe.push_back(laste);
			MSG("ClusterManager",Msg::kDebug)<<"\t\tpeak e "<<laste<<"\n";
			npeak++;
		}else{
			//these hits needs to be shared!
			sortedidx.push_back(lastidx);
			sortede.push_back(laste);		
		}
		
		
		//////////
		
		
		
		
		
		
		
		MSG("ClusterManager",Msg::kDebug)<<"cluster "<<theclusters[k].id<< " npeak "<<npeak<<"\n";
		
		if(npeak<2) //1 hit cluster will have 0 peaks!
		{
			theclusters[k].Finalize();
	
			if(theclusters[k].e>0.001)
				tmponepeak.push_back(theclusters[k]);
			continue;
		}
		
		
		
		//split up the cluster!
		for(unsigned int i=0;i<maxidx.size();i++)
		{
			c.Reset();
			c.AdvanceID();
		
		    //cout <<"idcounter 4 " <<c.idcounter<< " "<<c.id<<endl;
		
			double t=theclusters[k].hitt[maxidx[i]];
			double z=theclusters[k].hitz[maxidx[i]];
			double maxe=theclusters[k].hite[maxidx[i]];
			int plane=theclusters[k].hitplane[maxidx[i]];
			int strip=theclusters[k].hitstrip[maxidx[i]];
			int hitid=theclusters[k].hit_id[maxidx[i]];

			MSG("ClusterManager",Msg::kDebug)<<"filling m "<<z<<" "<<t<<" "<<maxe<<" "<<plane<<" "<<strip<<"\n";

			c.Insert(z, t, maxe, plane, strip, hitid);
		
			for(unsigned int j=0;j<sortedidx.size();j++)
			{
		
				//calculate normalizing denominator....
				double denom=0;
				for(unsigned int m=0;m<maxidx.size();m++)
				{
					double maxt=theclusters[k].hitt[maxidx[m]];
					double sortedt=theclusters[k].hitt[sortedidx[j]];
		
					denom+=theclusters[k].hite[maxidx[m]]/ ( (maxt-sortedt) * (maxt-sortedt) );
				}		
		
		
		
				double st = theclusters[k].hitt[sortedidx[j]];
				double e= sortede[j];// / ( (st - t ) * (st -t) );
				e*=maxe/((t-theclusters[k].hitt[sortedidx[j]]) * (t-theclusters[k].hitt[sortedidx[j]]))/denom;
				
					double t=theclusters[k].hitt[sortedidx[j]];
					double z=theclusters[k].hitz[sortedidx[j]];
					int plane=theclusters[k].hitplane[sortedidx[j]];
					int strip=theclusters[k].hitstrip[sortedidx[j]];
					int hitid=theclusters[k].hit_id[sortedidx[j]];
									
					MSG("ClusterManager",Msg::kDebug)<<"filling s "<<z<<" "<<t<<" "<<plane<<" "<<strip<<" --- "<<st<<" --- "<<maxe/denom<<" tot "<<sortede[j]<<" used "<<e<<"\n";
				
				if(e>0.001)
				{
					c.Insert(z, t, e, plane, strip,hitid);
				}
			}
			
			c.view = theclusters[k].view;
			c.Finalize();
			
			if(c.hitt.size()>0 && c.e>0.001)
			{
				tmponepeak.push_back(c);
				MSG("ClusterManager",Msg::kDebug)<<"new cluster stored with size "<<c.hitt.size()<<" e "<<c.e<<" in view "<<c.view<<"\n";
			}
		}
	
	}
	

	
	
	cluster_map.clear();
	cluster_map_u.clear();
	cluster_map_v.clear();
	
	for(unsigned int k=0;k<tmponepeak.size();k++)
	{
		cluster_map[tmponepeak[k].z][tmponepeak[k].t]=tmponepeak[k].id;
		tmponepeak[k].Finalize();
		clusters.push_back(tmponepeak[k]);
		if(tmponepeak[k].view==2)cluster_map_u[tmponepeak[k].z][tmponepeak[k].t]=tmponepeak[k].id;
		if(tmponepeak[k].view==3)cluster_map_v[tmponepeak[k].z][tmponepeak[k].t]=tmponepeak[k].id;
		MSG("ClusterManager",Msg::kDebug)<<"saving cluster id "<<tmponepeak[k].id<<"\n";
		
		
	//	minz=minz < tmponepeak[k].GetZ() ? minz : tmponepeak[k].GetZ();
	//	maxz=maxz > tmponepeak[k].GetZ() ? maxz : tmponepeak[k].GetZ();
	//	mint=mint < tmponepeak[k].GetT() ? mint : tmponepeak[k].GetT();
	//	maxt=maxt > tmponepeak[k].GetT() ? maxt : tmponepeak[k].GetT();
		
		
		
	}
	
	//printf("ClusterManager dimensions z,t %f %f z,t %f %f\n",minz,mint, maxz,maxt);
	
}




void ClusterManager::LoadHits()
{
	if(!hitmanager)return;
	Reset();
	
	std::vector<ManagedHit> avail_hits = hitmanager->GetAvailableHits();
	
	for(unsigned int i=0;i<avail_hits.size();i++)
	{
			
			hits.push_back(avail_hits[i]);
			loc_map[avail_hits[i].GetPlane()][avail_hits[i].GetStrip()]=hits.size()-1;
	
			double sz = avail_hits[i].GetZ();
			double st = avail_hits[i].GetT();
			
			minz=minz < sz ? minz : sz;
			maxz=maxz > sz ? maxz : sz;
			
			int view=avail_hits[i].GetView();
			if(view==2)
			{
				minu=minu < st ? minu : st;
				maxu=maxu > st ? maxu : st;
			}else if(view==3)
			{
				minv=minv < st ? minv : st;
				maxv=maxv > st ? maxv : st;			
			}
	}

}





void ClusterManager::AddStrip(int myplane, int mystrip, double myenergy, double st, double sz, int view)
{
	ManagedHit h(view,myplane,mystrip,sz,st,myenergy);
	
	hits.push_back(h);
	loc_map[myplane][mystrip]=hits.size()-1;
	
	
	minz=minz < sz ? minz : sz;
	maxz=maxz > sz ? maxz : sz;



			if(view==2)
			{
				minu=minu < st ? minu : st;
				maxu=maxu > st ? maxu : st;
			}else if(view==3)
			{
				minv=minv < st ? minv : st;
				maxv=maxv > st ? maxv : st;			
			}

	
}



void ClusterManager::Reset()
{
		clusters_are_made=0;
		
		loc_map.clear();
		cluster_map.clear();
		cluster_map_u.clear();
		cluster_map_v.clear();
		hits.clear();
		clusters.clear();
		
		minu=100000;
		maxu=-100000;
		minv=100000;
		maxv=-100000;
		minz=100000;
		maxz=-100000;
		
		needMapRebuild=0;
		clusters_to_delete.clear();
		
		inuse.clear();
}



Managed::ManagedCluster* ClusterManager::GetCluster(int cid)
{
	if(cid==0)return 0;
	if(cid<0) //means it is a saved cluster
	{
		if(clustersaver)return clustersaver->GetCluster(cid);
		else return 0;
	}

	ManagedCluster * c =0;
	for(unsigned int i=0;i<clusters.size();i++)
	{
		if(clusters[i].id==cid)c=&clusters[i];
	}

	for(unsigned int i=0;i<inuse.size();i++)
	{
		if(inuse[i]==cid)return 0;
	}

	return c;
}

Managed::ManagedCluster* ClusterManager::GetSavedCluster(int cid)
{
	return clustersaver->GetCluster(cid);
}

//take 2 clusters, merge them, and make a new cluster from them (returning the new cluster id)
//the original 2 clusters must not be used (and are removed from the ClusterManager)
int ClusterManager::MergeClusters(int cid1, int cid2)
{
	ManagedCluster *c1 = GetCluster(cid1);
	ManagedCluster *c2 = GetCluster(cid2);
	
	if(!c1 && c2)return cid2;
	if(!c2 && c1)return cid1;
	if(!c1 && !c2)return -1;
	
	
	ManagedCluster newcluster;
	MSG("ClusterManager",Msg::kDebug)<<"merging "<<cid1<<" "<<cid2<< " making "<<newcluster.id<<"\n";
	
	//copy over the old hits...
	for(unsigned int i=0;i<c1->hit_id.size();i++)
		newcluster.Insert(c1->hitz[i],c1->hitt[i],c1->hite[i],c1->hitplane[i],c1->hitstrip[i],c1->hit_id[i]);
	for(unsigned int i=0;i<c2->hit_id.size();i++)
		newcluster.Insert(c2->hitz[i],c2->hitt[i],c2->hite[i],c2->hitplane[i],c2->hitstrip[i],c2->hit_id[i]);

//	printf("%d %d\n",c1->id,c2->id);


	//the pointers c1 and c2 seem to be disrupted by the assignment of the newcluster to the maps... why?
	//because the pointers are into a vector... which has its elements shifted!
	
	EraseCluster(cid1);
	EraseCluster(cid2);
		
	cluster_map[newcluster.z][newcluster.t]=newcluster.id;
	newcluster.Finalize();
	clusters.push_back(newcluster);
	if(newcluster.view==2)cluster_map_u[newcluster.z][newcluster.t]=newcluster.id;
	if(newcluster.view==3)cluster_map_v[newcluster.z][newcluster.t]=newcluster.id;


	


	return newcluster.id;
}
		

//for internal use only...
//hits owned by this cluster should first be assigned elsewhere before erasing the cluster!		
void ClusterManager::EraseCluster(int cid)
{

	MSG("ClusterManager",Msg::kDebug)<<"erasing "<<cid<<"\n";

    std::map<double, std::map<double, int> >::iterator p_iterr;
    std::map<double, int >::iterator s_iterr;
	for(p_iterr=cluster_map.begin();p_iterr!=cluster_map.end(); p_iterr++)
    {	
    	int done =0;
    	for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
        {
        	if(s_iterr->second==cid)
        	{
        		p_iterr->second.erase(s_iterr);
        		done=1;
        		break;
			}
		}
		if(done)break;
	}	

	for(p_iterr=cluster_map_u.begin();p_iterr!=cluster_map_u.end(); p_iterr++)
    {	
    	int done =0;
    	for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
        {
        	if(s_iterr->second==cid)
        	{
        		p_iterr->second.erase(s_iterr);
        		done=1;
        		break;
			}
		}
		if(done)break;
	}	
	
	for(p_iterr=cluster_map_v.begin();p_iterr!=cluster_map_v.end(); p_iterr++)
    {	
    	int done =0;
    	for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
        {
        	if(s_iterr->second==cid)
        	{
        		p_iterr->second.erase(s_iterr);
        		done=1;
        		break;
			}
		}
		if(done)break;
	}	

	
	std::vector<Managed::ManagedCluster>::iterator it;
	for(it=clusters.begin();it!=clusters.end();it++)
	{
		if(it->id==cid)
		{
			clusters.erase(it);
			break;
		}
	}

}




void ClusterManager::FillClusterMap(std::map<double, std::map<double, std::pair<double, int> > > * cluster_map	)
{
	for(unsigned int i=0;i<clusters.size();i++)
	{
		(*cluster_map)[clusters[i].z][clusters[i].t]=std::pair<double,int>(clusters[i].e, clusters[i].view);
	}

}





