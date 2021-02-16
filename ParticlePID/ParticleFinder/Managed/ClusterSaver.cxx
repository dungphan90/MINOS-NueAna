#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterSaver.h"

#include <math.h>
using namespace Managed;
ClassImp(ClusterSaver)

ClusterSaver::ClusterSaver()
{
		Reset();
		
}


ClusterSaver::~ClusterSaver()
{

		cluster_map.clear();
		clusters.clear();	
		cluster_map_u.clear();
		cluster_map_v.clear();	
		
		clusters_to_delete.clear();
}


int ClusterSaver::SaveCluster(Managed::ManagedCluster *cluster)
{

	if(cluster->e<0.0001)return 0;//don't save an empty cluster

	Managed::ManagedCluster c = *cluster;
	c.id = --save_id;
	cluster_map[c.z][c.t]=c.id;
	c.Finalize();
	clusters.push_back(c);
	
	cluster_map[c.z][c.t]=c.id;
	
	if(c.view==2)cluster_map_u[c.z][c.t]=c.id;
	if(c.view==3)cluster_map_v[c.z][c.t]=c.id;


	if(c.z<minz)minz=c.z;
	if(c.z>maxz)maxz=c.z;
	if(c.view==2)
	{
		if(c.t<minu)minu=c.t;
		if(c.t>maxu)maxu=c.t;
	}
	if(c.view==3)
	{
		if(c.t<minv)minv=c.t;
		if(c.t>maxv)maxv=c.t;
	}

	if(c.t<mint)mint=c.t;
	if(c.t>maxt)maxt=c.t;

	return c.id;
}


std::map<double, std::map<double, int>  > * ClusterSaver::GetClusterMap(int view)
{

	//need rebuild?
	
	if(needMapRebuild)RebuildClusterMaps();



	if(view==2)return &cluster_map_u;
	if(view==3)return &cluster_map_v;
	return &cluster_map;
}



void ClusterSaver::RebuildClusterMaps()
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




void ClusterSaver::DumpClusters()
{

	RebuildClusterMaps();
/*	
	printf("array dump\n");
	for(int i=0;i<clusters.size();i++)
	{		ManagedCluster *mc = &clusters[i];
			if(!mc){printf("missing cluster id %d\n",i);continue;}
			printf("%d  z %f t %f e %f dz %f dt %f view %d \n",i,mc->z,mc->t,mc->e,mc->dz,mc->dt,mc->view);	
	}
*/
	std::map<double, std::map<double, int> >::iterator p_iter;
	std::map<double, int>::iterator s_iter;

	printf("dumping clusters... %d found in array\n",(int)clusters.size());
           
	for(p_iter=cluster_map.begin();p_iter!=cluster_map.end(); p_iter++)
	{
		std::map<double, int>::iterator s_iter;
		for(s_iter=p_iter->second.begin();s_iter!=p_iter->second.end(); s_iter++)
		{
			ManagedCluster *mc = GetCluster(s_iter->second);
			if(!mc){printf("missing cluster id %d\n",s_iter->second);continue;}
			printf("z %f t %f e %f dz %f dt %f view %d \n",mc->z,mc->t,mc->e,mc->dz,mc->dt,mc->view);
			printf("id %d\n",s_iter->second);
	
		}
	}

	printf("done....\n\n");
}



void ClusterSaver::Reset()
{

		
		cluster_map.clear();
		clusters.clear();
		cluster_map_u.clear();
		cluster_map_v.clear();		
		mint=100000;
		maxt=-100000;
		minz=100000;
		maxz=-100000;
		minu=100000;
                maxu=-100000;
                minv=100000;
                maxv=-100000;
		nClusters=0;

		save_id=-1;
		needMapRebuild=0;
		clusters_to_delete.clear();		
}



Managed::ManagedCluster* ClusterSaver::GetCluster(int cid)
{
	ManagedCluster * c =0;
	for(unsigned int i=0;i<clusters.size();i++)
	{
		if(clusters[i].id==cid)c=&clusters[i];
	}

	return c;
}



void ClusterSaver::FillClusterMap(std::map<double, std::map<double, std::pair<double, int> > > * cluster_map	)
{
	for(unsigned int i=0;i<clusters.size();i++)
	{
		(*cluster_map)[clusters[i].z][clusters[i].t]=std::pair<double,int>(clusters[i].e, clusters[i].view);
	}

}


std::map<std::pair<int,int>, double> ClusterSaver::GetStripEnergy()
{
	std::map<std::pair<int,int>, double> ret;

	for(unsigned int i=0;i<clusters.size();i++)
	{
		for(unsigned int j=0;j<clusters[i].hitplane.size();j++)
		{
			if(clusters[i].GetStatus()<=-10)continue;

			//printf("id %d p %d s %d e %f\n",clusters[i].id,clusters[i].hitplane[j],clusters[i].hitstrip[j],clusters[i].hite[j]);
	
			ret[std::pair<int,int>(clusters[i].hitplane[j],clusters[i].hitstrip[j])]+=clusters[i].hite[j];

		}
	}
	return ret;
}


void ClusterSaver::recomputeBounds()
{

                mint=100000;
                maxt=-100000;
                minz=100000;
                maxz=-100000;
                minu=100000;
                maxu=-100000;
                minv=100000;
                maxv=-100000;
		nClusters=0;

        for(unsigned int i=0;i<clusters.size();i++)
        {
                for(unsigned int j=0;j<clusters[i].hitplane.size();j++)
                {
                        if(clusters[i].GetStatus()<=-10)continue;

			ManagedCluster c = clusters[i];
			if(c.e<1e-6)continue;

			nClusters++;

		    	if(c.z<minz)minz=c.z;
		        if(c.z>maxz)maxz=c.z;
		        if(c.view==2)
        		{               
                		if(c.t<minu)minu=c.t;
                		if(c.t>maxu)maxu=c.t;
        		}
        		if(c.view==3)
        		{
       		        	if(c.t<minv)minv=c.t;
                		if(c.t>maxv)maxv=c.t;
        		}

        		if(c.t<mint)mint=c.t;
        		if(c.t>maxt)maxt=c.t;

                    

                }
        }

}



