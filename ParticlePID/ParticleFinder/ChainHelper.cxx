#include "MessageService/MsgService.h"
#include "NueAna/ParticlePID/ParticleFinder/ChainHelper.h"
#include <math.h>
#include <map>
#include <algorithm>

#include <sstream>



CVSID("$Id: ChainHelper.cxx,v 1.4 2010/07/30 03:46:50 rhatcher Exp $");


ClassImp(ChainHelper)



ChainHelper::ChainHelper() : parents(0),muonlikechains(0),emlikechains(0),totalchains(0), max_z_gap(0.3), max_slope(0.0),
								found_max_path(0)
{
	Reset();

	working.clear();
	finished.clear();
	pending_t.clear();
	pending_z.clear();
	pending_e.clear();
	pending_cluster_id.clear();
	maxpath.clear();
	
	
	Chain::ResetCounter();
	
	idhelper.clear();
	
	detector=Detector::kUnknown;
	
}

int ChainHelper::NewChain() // for manually adding a chain to finished...
{
	Chain c;
	AddFinishedChain(c);
	return c.myId;

}


void ChainHelper::AddFinishedChain(Chain c)
{
	//printf("ChainHelper adding chain %d\n",c.myId);
	finished.push_back(c);
	idhelper[c.myId]= finished.size()-1 ;	
}


ChainHelper::~ChainHelper()
{
	working.clear();
	finished.clear();
	pending_t.clear();
	pending_z.clear();
	pending_e.clear();
	pending_cluster_id.clear();
	idhelper.clear();
	maxpath.clear();
}


Chain* ChainHelper::GetChain(int id)
{
	for(unsigned int i=0;i<finished.size();i++)
	if(id == finished[i].myId)return &finished[i];


/*	std::map<int,int>::iterator it;
	it = idhelper.find(id);
	
	if(it!=idhelper.end())
		return &finished[it->second];   ///does not check to see if it is valid!!!!!
*/	
	return 0;
}


void ChainHelper::DeleteChain(int id)
{
	
	//cout << "deleting "<<id<<endl;
	
	Chain* todel = GetChain(id);
	if(!todel)return;
	
	if(todel->parentChain>-1)
	{
		Chain*p = GetChain(todel->parentChain);
		
		for(unsigned int i=0;i<p->children.size();i++)
		{
			if(p->children[i]==id)
			{
				p->children.erase(p->children.begin()+i);
			
				//cout<<"\t removing child reference from parent "<<todel->parentChain<<endl;
			
				break;	
			}
		}
	}
	
	
	for(unsigned int i=0;i<todel->children.size();i++)
	{
		Chain * child = GetChain(todel->children[i]);
		if(!child)continue;
			child->parentChain=-1;
	}
	
	
	for(unsigned int i=0;i<finished.size();i++)
	{
		if(finished[i].myId==id)
		{
			finished.erase(finished.begin()+i);	
			break;
		}
	}
	
	
	//rebuild helper!!!!
	
	idhelper.clear();
	for(unsigned int i=0;i<finished.size();i++)
		idhelper[finished[i].myId]=i;
		
	found_max_path=0;  //in case we deleted a chain in the max path
	
}



void ChainHelper::Reset()
{
	working.clear();
	finished.clear();
	pending_t.clear();
	pending_z.clear();
	pending_e.clear();
	pending_cluster_id.clear();
	
		Chain::ResetCounter();
		
		found_max_path=0;
		parents=0;
		muonlikechains=0;
		emlikechains=0;

	vtx_z=0;
	vtx_t=0;

}


void ChainHelper::add_to_plane(double it, double iz, double ie, int cluster_id)
{
	pending_t.push_back(it);
	pending_z.push_back(iz);
	pending_e.push_back(ie);
	pending_cluster_id.push_back(cluster_id);
	
	//reset found items, since we have added data
	found_max_path=0;

}

void ChainHelper::matchChains()
{

	found_max_path=0;

	std::vector<int> longs;
	
	for(unsigned int i=0;i<finished.size();i++)
	{
		if(finished[i].entries>3) longs.push_back(i);
	}
	
	for(unsigned int i=0;i<longs.size();i++)
	for(unsigned int j=i+1;j<longs.size();j++)
	{
		if(finished[longs[j]].start_z-finished[longs[i]].end_z > 1.0)continue;
		if(finished[longs[j]].end_z-finished[longs[j]].start_z <0.0)continue;
		if(finished[longs[i]].end_z-finished[longs[i]].start_z <0.0)continue;	
		if(finished[longs[j]].parentChain>-1)continue;	
		
		//else see if there is a match....
		
		double dz=finished[longs[j]].start_z-finished[longs[i]].end_z;
		
//		double t1 = finished[longs[i]].end_t+finished[longs[i]].back_slope*dz;
//		double t2 = finished[longs[j]].start_t-finished[longs[j]].front_slope*dz;

		double t1 = finished[longs[i]].interpolate(dz);
		double t2 = finished[longs[j]].interpolate(dz);

		
		
		double dt1 = fabs(finished[longs[j]].start_t - t1);
		double dt2 = fabs(finished[longs[i]].end_t - t2);
		
		
		MSG("ChainHelper",Msg::kDebug) <<"matching ... " << dt1 <<" "<<dt2<<endl;
				
		if(dz<1.0 && (dt1<0.1 || dt2<0.2)  )
		{
			finished[longs[i]].children.push_back(finished[longs[j]].myId);
			finished[longs[j]].parentChain=finished[longs[i]].myId;
		
		}	
	
	}




}


void ChainHelper::process_plane()
{

	int wsize=working.size();


	std::vector<int> todel(wsize);
	for(int i=0;i<wsize;i++)todel[i]=0;
		
	
	std::map<int, std::vector<int> > toadd;
	//std::vector<int> toadd(wsize);
	//for(int i=0;i<wsize;i++) toadd[i]=0;
	
	std::vector<Chain> tmpchain;
	


	//printf("%d clusters for this plane\n",pending_t.size());
	
	for(unsigned int j=0;j<pending_t.size();j++)
	{
	
		double iz=pending_z[j];
		double it=pending_t[j];
		
			
		//if we are in the near... and in the calorimeter, we should increaze the z gap size
		double max_z_gap_temp=max_z_gap;
		if(detector==Detector::kNear && iz > 6)
		{
			max_z_gap_temp*=5;
		}
	
		double closest_d=100000;
		int closest_index=-1;
		int match=0;
		
		double closest_t_unmatched=100000;
		double closest_t=100000;
		
		int singlehitmatch=0;
		double singlehite=0;
		
		for(unsigned int i=0;i<working.size();i++)
		{
			if(fabs(iz-working[i].end_z)<0.0001) //same plane
			{
			//	printf("same plane %f\n",iz);
				continue;
			}
		
		
			if(fabs(iz-working[i].end_z)>max_z_gap_temp)
			{
				//printf("removing chain %d from working list - it's too far away\n",working[i].myId);
				todel[i]=1;
				continue;
			}
			
			double tslope = (working[i].end_t-it)/(working[i].end_z-iz);
		
			///max dist
				double ddist = sqrt((working[i].end_t-it)*(working[i].end_t-it) + (working[i].end_z-iz)*(working[i].end_z-iz));
				
				if(ddist>max_z_gap_temp)continue;
			
			
			///end max dist
		
		
		
			closest_t_unmatched=closest_t_unmatched<fabs(working[i].end_t-it) ?closest_t_unmatched:fabs(working[i].end_t-it);
		
		//	if(working[i].children.size()>0)continue; //should be finalized...
		
			if(fabs(tslope-working[i].last_slope)<max_slope/(fabs(working[i].end_z-iz)/0.0708) || !working[i].good_slope() || working[i].entries<2 )// || fabs(working[i].last_slope) > 100000)
			{
				match=1;
				
				if(working[i].good_slope() && fabs(working[i].last_slope)<10000)match=2;
				
				
				double dist = closest_t;
				
		//		if(working[i].good_slope())
		//		 dist = sqrt((working[i].end_t-it)*(working[i].end_t-it)+(tslope-working[i].last_slope)*(tslope-working[i].last_slope) );
		//		 else
				 dist = fabs(working[i].end_t-it);
				
				//closest_t=closest_t<fabs(working[i].end_t-it) ?closest_t:fabs(working[i].end_t-it);
				
				
				closest_t=closest_t<dist?closest_t:dist;
				
				
				
				double ldt=(working[i].end_t-it)*(working[i].end_t-it)+(working[i].end_z-iz)*(working[i].end_z-iz);
				
				double kinky=0.0;
				if(working[i].entries>1)
				{
					int a = working[i].entries-2;
					
					
					if (! fabs(working[i].z[a+1]-working[i].z[a]) < 0.00001)
					{
						double slope = (working[i].t[a+1]-working[i].t[a]) /  (working[i].z[a+1]-working[i].z[a]);
						double off = working[i].t[a];
						double z = working[i].z[a];
						kinky =  fabs ( ( (iz - z) * slope + off ) - it);
						
						kinky = fabs ( working[i].interpolate(iz)-it);
						
					}
		
				}
				
			
				if(ldt<closest_d ||singlehitmatch==1)  //|| match==2)
				if(kinky < 0.1)
				{
					//if this is a match to a single hit...
					//make sure that there is not a single hit match with greater energy!				
					//is single hit match?
					if(working[i].entries==1)
					{
						if(singlehitmatch==0)
						{
							singlehitmatch=1;
							singlehite=working[i].e[0];
							closest_d=ldt;
							closest_index=i;
						}else{
							if(singlehite <working[i].e[0])
							{
								singlehite = working[i].e[0];
								closest_d=ldt;
								closest_index=i;
							}
						}
					
					}else{
						closest_d=ldt;
						closest_index=i;
						singlehitmatch=0;
					}
				}
			}
			
			//printf("chain id %d  cluster idx %d  match %d\n",working[i].myId,j,match);
		}

		
		int notyetdone=1;
		//see if best match is more than 1 plane away....
		//if it is, see if there is a working chain that had its tail at 1 plane away....
		//but a hit in that chain at the same distance as this match is a better candidate....
		if(closest_index>-1)
		{
			double zdist = fabs (working[closest_index].end_z-iz);
			double zplanedist = zdist / 0.035;
			if(zplanedist>1.4)
			{
				std::vector<int> matchidx;
				std::vector<int> matchentry;
				std::vector<double> matchdist;
				//iterate over working chains, look for hit in same plane as current best guess
				for(unsigned int i=0;i<working.size();i++)
				{
					if(closest_index==(int)i)continue;//we want to consider other chains...
				
					//consider chains with 2+ hits
					if(working[i].entries<2)continue;
					if(fabs(iz-working[i].end_z) * 1.2 <zdist)
					{
						//see if chain has hit in same plane as current candidate
						for(int j=0;j<working[i].entries;j++)
						{
							if(fabs(working[i].z[j] - working[closest_index].end_z)< 0.0175) //half a plane difference....
							{
								matchidx.push_back(i);
								matchentry.push_back(j);
								double dist2 = (working[i].z[j]-iz)*(working[i].z[j]-iz) + (working[i].t[j]-it)*(working[i].t[j]-it) ;
								matchdist.push_back(dist2);
							}
						
						}
					
					}			
				}
				
				int tmpi=-1;
				for(unsigned int i =0 ;i<matchdist.size();i++)
				{
					if(matchdist[i]< closest_d)
					{
						tmpi=i;
						closest_d=matchdist[i];
					}
				}
				
				if(tmpi>-1)
				{
					//make a new chain with this hit.. and attach it to the chain
					Chain c;
					MSG("ChainHelper",Msg::kDebug) <<"starting new chain, attaching to old chain "<< c.myId << ", " << working[matchidx[tmpi]].myId <<endl;
					c.add_to_back(pending_t[j],pending_z[j],pending_e[j],pending_cluster_id[j]);
					//printf("!!! %d\n",pending_cluster_id[j]);
					
					
					tmpchain.push_back(AttachAt(&working[matchidx[tmpi]], &c,working[matchidx[tmpi]].z[matchentry[tmpi]]));
					tmpchain.push_back(c);
				
					notyetdone=0;
				}
				
				
				
			} 
		}


	//this should consider direct pointing over closest match....
	//	if(closest_index<0)
	//	{
		//make sure that there are no chains clearly pointing to this hit
		double pdistr=1000;
		int tmp_closest_index=-1;
		for(unsigned int i=0;i<working.size();i++)
		{
			if(closest_index==(int)i)continue;//we want to consider other chains...
				
			//consider chains with 2+ hits
			if(working[i].entries<2)continue;
		
			//double tpos = working[i].back_offset + working[i].back_slope*iz; 
			double tpos =working[i].interpolate(iz);
		
		
			if(fabs(it-tpos)<0.025*1.5)
			{
				double rdist = (working[i].end_z-iz)*(working[i].end_z-iz)+(working[i].end_t-it)*(working[i].end_t-it);
				
				if(rdist < pdistr)
				{
					tmp_closest_index=i;
					pdistr=rdist;
				}
			}
			
		}

		
		if(pdistr<max_z_gap_temp && tmp_closest_index>-1)closest_index=tmp_closest_index;
	//	}


		if(notyetdone)
		{
		
			if(closest_index>-1)
			{
				toadd[closest_index].push_back(j);
			}else{

				//are we close to the vertex?
				//printf("pending z t %f %f vtx %f %f\n",pending_z[j],pending_t[j],vtx_z,vtx_t);
				if(vtx_z==0 || (fabs(pending_t[j]-vtx_t)<0.2 && fabs(pending_z[j]-vtx_z)<0.2)){
			
				Chain c;
				MSG("ChainHelper",Msg::kDebug) <<"no match, starting new chain "<< c.myId <<endl;
				c.add_to_back(pending_t[j],pending_z[j],pending_e[j],pending_cluster_id[j]);
				//printf("!!! %d\n",pending_cluster_id[j]);
				//working.push_back(c);
				tmpchain.push_back(c);
				}
			}

		}

	}
	
	
	for(unsigned int i=0;i<tmpchain.size();i++)	
	working.push_back(tmpchain[i]);




//	for(int i=0;i<wsize;i++)
//	{

    std::map<int, std::vector<int> >::iterator iter;
	for(iter=toadd.begin();iter!=toadd.end();iter++)
	{

		//int s = toadd[i].size()>0;
		int s = iter->second.size();
		
		
		if(s==1)
		{
			int item=iter->second[0];
			
			if(working[iter->first].children.size()>0)
			{
				Chain *t=split(&working[iter->first]);
				t->add_to_back(pending_t[item],pending_z[item],pending_e[item],pending_cluster_id[item]);
				working.push_back(*t);
				delete t;t=0;
				MSG("ChainHelper",Msg::kDebug) <<"---splitting chain "<< working[iter->first].myId <<" to make "<< t->myId <<endl;
			
			}else{
				working[iter->first].add_to_back(pending_t[item],pending_z[item],pending_e[item],pending_cluster_id[item]);
			}
			
			
			MSG("ChainHelper",Msg::kDebug)<< "adding to chain "<<working[iter->first].myId<<endl;
		}
		
		if(s>1)
		{
			//need to copy the chains
			for(int k=0;k<s;k++) //s-1
			{
				Chain *t=split(&working[iter->first]);
				int item=iter->second[k];
				t->add_to_back(pending_t[item],pending_z[item],pending_e[item],pending_cluster_id[item]);
				working.push_back(*t);
				delete t;t=0;
				MSG("ChainHelper",Msg::kDebug)<< "---splitting chain " << working[iter->first].myId << " to make " <<  t->myId <<endl;
			}
		//	int item=iter->second[s-1];
	//		working[iter->first].add_to_back(pending_t[item],pending_z[item],pending_e[item],pending_cluster_id[item]);
//			MSG("ChainHelper",Msg::kDebug)<<"adding to chain"<<endl;
		}
	}


	
	
	pending_t.clear();
	pending_z.clear();
	pending_e.clear();
	pending_cluster_id.clear();
	found_max_path=0;

}

Chain * ChainHelper::split(Chain * d){

	Chain * r = new Chain();
	r->parentChain=d->myId;
	r->start_t=d->end_t;
	r->start_z=d->end_z;
	r->level=d->level+1;

	int dsize=d->t.size()-1;
	r->add_to_back(d->t[dsize],d->z[dsize],d->e[dsize],d->cluster_id[dsize]);

	d->children.push_back(r->myId);


//cout <<"????????? new id "<<r->myId << " split from "<<d->myId<<endl;

	return r;
}

Chain  ChainHelper::AttachAt(Chain *c, Chain * daughter, double splitpointz)
{
//	if(c->end_z < splitpointz || c->start_z > splitpointz || c->entries<1)
//	{
	if(c->entries<1 || splitpointz < c->start_z || splitpointz > daughter->end_z){  //its ok if the split point is past the parent chain... then we just attach to the end!

		cout<<"!!!!!!!!!THIS SHOULDN'T HAPPEN --- CALL TO ATTACHAT IN CHAINHELPER WITH SPLITPOINT OUTSIDE OF PARENT CHAIN!!!!"<<endl;
		cout<<"requested to connect "<<daughter->myId <<" to "<<c->myId <<" at "<<splitpointz<<endl;
		cout<<"daughter\n";
		daughter->PrintChain();
		cout<<"parent\n";
		c->PrintChain();
		cout<<"\n\n";
		return Chain(); // can't split up the chain if we are not in the chain...
	}
	std::vector<double> t;
		for(unsigned int i=0;i<c->t.size();i++)t.push_back(c->t[i]);
	std::vector<double> e;
		for(unsigned int i=0;i<c->e.size();i++)e.push_back(c->e[i]);
	std::vector<double> z;
		for(unsigned int i=0;i<c->z.size();i++)z.push_back(c->z[i]);
		
	std::vector<int> cid;
		for(unsigned int i=0;i<c->cluster_id.size();i++)cid.push_back(c->cluster_id[i]);	
		
	//int pid = c->parentChain;
	std::vector<int> children;
		for(unsigned int i=0;i<c->children.size();i++)children.push_back(c->children[i]);
	
	
	
	
	double close = fabs(z[0]-splitpointz);
	int cidx=0;
	//find the split point
	for(unsigned int i=1;i<z.size();i++)
	{
		if(fabs(z[i]-splitpointz)<close)
		{
			close = fabs(z[i]-splitpointz);
			cidx=i;
		}
	}

	Chain newc;
	
	c->ClearHits();
	for(int i=0;i<=cidx;i++)
	{
		c->add_to_back(t[i],z[i],e[i],cid[i]);
	}
	
	c->children.clear();
	c->children.push_back(newc.myId);
	
	for(unsigned int i=cidx;i<z.size();i++)
	{
		newc.add_to_back(t[i],z[i],e[i],cid[i]);
	}
	
	newc.parentChain = c->myId;
	newc.children = children;
	
	newc.level = c->level+1;
	
	AttachAsChild(c, daughter);
	
	
	
	return newc;


}


Chain * ChainHelper::SplitAt(Chain * c, double splitpointz)
{


	//printf("Request to split chain %d\n",c->myId);

	if(c->end_z < splitpointz || c->start_z > splitpointz || c->entries<1) 
	{
		printf("CANT SPLIT CHAIN - split point out of chain range\n");
		return c; // can't split up the chain if we are not in the chain...
	}
	
	std::vector<double> t;
		for(unsigned int i=0;i<c->t.size();i++)t.push_back(c->t[i]);
	std::vector<double> e;
		for(unsigned int i=0;i<c->e.size();i++)e.push_back(c->e[i]);
	std::vector<double> z;
		for(unsigned int i=0;i<c->z.size();i++)z.push_back(c->z[i]);
		
	std::vector<int> cid;
		for(unsigned int i=0;i<c->cluster_id.size();i++)cid.push_back(c->cluster_id[i]);		
		
	//int pid = c->parentChain;
	std::vector<int> children;
		for(unsigned int i=0;i<c->children.size();i++)children.push_back(c->children[i]);
	


	
	
	double close = fabs(z[0]-splitpointz);
	int cidx=0;
	//find the split point
	for(unsigned int i=0;i<z.size();i++)
	{
		if(fabs(z[i]-splitpointz)<close)
		{
			close = fabs(z[i]-splitpointz);
			cidx=i;
		}
	}
	
	
	if (cidx==0 || cidx == (int)z.size()-1)return c; //cant split a chain at the first or last hit
	
	//printf("splitting at %d\n",cidx);




	Chain newc;
	
	c->ClearHits();
	for(int i=0;i<=cidx;i++)
	{
		c->add_to_back(t[i],z[i],e[i],cid[i]);
	}






	c->children.clear();
	c->children.push_back(newc.myId);

	
	for(unsigned int i=cidx;i<z.size();i++) //duplicate the head cluster
	{
		newc.add_to_back(t[i],z[i],e[i],cid[i]);
	}
	
	newc.parentChain = c->myId;
	for(unsigned int child = 0;child < children.size();child++)
		newc.children.push_back(children[child]);
	
	newc.level=c->level+1;
	
	//printf("storing....\n");


	int  retid = c->myId;
	AddFinishedChain(newc);

	
	
	//printf("done\n");
	return GetChain(retid);  //reget the pointer...... adding to the finished vector can change the addresses...

}

		
void ChainHelper::finish()
{
	for(unsigned int i=0;i<working.size();i++)
	{
		AddFinishedChain(working[i]);
		
	}
	
	
	for(unsigned int i=0;i<finished.size();i++)
	{
		if(finished[i].parentChain<0)parents++;
		if(finished[i].muonfrac>0.8)muonlikechains++;
		if(finished[i].emlike>0.7)emlikechains++;
	}
	working.clear();
	
	std::sort(finished.begin(), finished.end(), LTId() ); 
	
	//matchChains();
	
	totalchains=finished.size();
	
	
}




void ChainHelper::print_finished()
{

	if(!MsgService::Instance()->IsActive("ChainHelper",Msg::kDebug))return;
	
	cout<<"chain helper print_finished()\n";
	
	
	for(unsigned int i=0;i<finished.size();i++)
	{
		
		std::ostringstream os;
	
	
		os <<"\n chain " <<finished[i].myId <<", parent " <<finished[i].parentChain <<", level " <<  finished[i].level <<", children ";

		
		for (unsigned int j=0;j<finished[i].children.size();j++)
			os<< finished[i].children[j] << " ";
			os << "\n";
		
			os << " muonfrac " << finished[i].muonfrac <<", emlike " << finished[i].emlike <<", numpeaks " <<finished[i].num_peaks <<"  strick inc " << finished[i].strict_increasing<<" dec " << finished[i].strict_decreasing<<"  entries " << finished[i].entries<<endl;
		
		
		
	
			os  << "start (t,z)  (" << finished[i].start_t <<", "<< finished[i].start_z <<")   end (t,z)  (" << finished[i].end_t <<", "<<  finished[i].end_z <<")  avg slope " << finished[i].avg_slope<< " avg offset " << finished[i].avg_offset<<" last slope " << finished[i].last_slope<<"  sum_e " <<finished[i].sum_e <<" avg_e " << finished[i].sum_e/finished[i].entries << "\n\t"; 
			
			os << "front slope " << finished[i].front_slope << " front offset " << finished[i].front_offset << " back slope " << finished[i].back_slope <<  " back offset " << finished[i].back_offset<<"\n\t";
			
		
			for(unsigned int k=0;k<finished[i].t.size();k++)
			{
				os << "(  " << finished[i].t[k] <<", "<< finished[i].z[k] << ", "<< finished[i].e[k] <<" - "<< finished[i].cluster_id[k]<<") ";
			}			
		
	
		
		cout<< os.str() << endl<<endl;

	}


}



//find the continguous set of chains with the largest energy
std::vector<int> ChainHelper::FindMaxPath()
{
	if(found_max_path)return maxpath; //dont rerun it if it is already done

	//make a list of all "parent" chains (not a daughter of another chain)
	//for each in the list, 
	
	double maxe=0;
	std::pair< std::vector<int>, double > maxp;
	for(unsigned int i=0;i<finished.size();i++)
	{
		if(finished[i].parentChain<0)
		{
		
			std::pair< std::vector<int>, double > f = FindMaxPath(GetChain(finished[i].myId));
			if(f.second>maxe)
			{
				maxe=f.second;
				maxp=f;
			}		
		}
	}
	
	
	ostringstream os;
	
	os << "\nmaxpath found with energy "<< maxe <<endl;

	os << "\t Chains : ";
	for(int i=maxp.first.size()-1;i>-1;i--)
	{
		int idx = maxp.first[i];
		os << idx <<" ";
	}
	os<<endl;	

	os << "\t Chain path : ";
	
	for(int i=maxp.first.size()-1;i>-1;i--)
	{
		int idx = maxp.first[i];
		
		os << " from (" << GetChain(idx)->start_z << " " <<  GetChain(idx)->start_t <<") to (" <<GetChain(idx)->end_z<< " " <<GetChain(idx)->end_t<< ") - "; 
	}
	
	os << "\n\t Chain energies : ";
	
	if(maxp.first.size()>0)
		os <<" "<< GetChain(maxp.first[maxp.first.size()-1])->e[0] << " - ";
		
	for(int i=maxp.first.size()-1;i>-1;i--)
	{
		int idx = maxp.first[i];
		
		
			for(unsigned int a=1;a<GetChain(idx)->e.size();a++) //start with 2nd in each chain to avoid double printing of first hit!
			os << " " << GetChain(idx)->e[a] << " - "; 
			
			os << " --- ";
					
					
	}

	MSG("ChainHelper",Msg::kDebug) << os.str() << endl;
		
	
	//since this vector was made recursively, the entries are backwards (from the end of the path to the beginning...) so reverse it
	
	std::vector<int> rev;
	for(int i=maxp.first.size()-1;i>-1;i--)
	rev.push_back(maxp.first[i]);
	
	
	found_max_path=1;  
	
	maxpath = rev;
	
	return rev;
}


//for recursively 
std::pair< std::vector<int>, double> ChainHelper::FindMaxPath(Chain *c)
{


//	printf("find max path chain id %d\n",c->myId);	
		std::vector<int> v;
		double e;

	if (c->children.size()<1)
	{
		v.push_back(c->myId);
		e=c->sum_e;
		return std::make_pair(v,e);
	}else{
	
	
		double maxe=0;
		std::pair< std::vector<int>, double > maxp;
		for(unsigned int i=0;i<c->children.size();i++)
		{
			Chain *nextchain = GetChain(c->children[i]);
			if(!nextchain)
			{
				MSG("ChainHelper",Msg::kWarning)<<"Call to child chain that does not exists!...."<<endl;
				continue;
			}
			std::pair< std::vector<int>, double > f = FindMaxPath(nextchain);
			//printf("maxpath on %d has e %f\n",nextchain->myId,f.second);
			if(f.second>maxe)
			{
				maxe=f.second;
				maxp=f;
			}		
		}

	
		maxp.first.push_back(c->myId);
		maxp.second+=c->sum_e;
		
		
		
	
		return maxp;
	}

}


void ChainHelper::insert(double it, double iz, double ie, int my_cluster_id)
{
	//find best chain match, and attach it
	int match=0;
	
	
	double closest_dt=100000;
	double closest_index=-1;
	
	
	double max_z_gap_temp=max_z_gap;
	
	if(detector==Detector::kNear && iz > 6)
	{
		max_z_gap_temp*=5;
	}
	
	std::vector<int> todel;
	for(unsigned int i=0;i<working.size();i++)
	{
		if(fabs(iz-working[i].end_z)>max_z_gap_temp)
		{
			//too far away in z, so mark this chain as finished
			AddFinishedChain(working[i]);
			todel.push_back(i);
			continue;
		}
		
		
		double tslope = (working[i].end_t-it)/(working[i].end_z-iz);
		
		if(fabs(tslope-working[i].last_slope)<max_slope)
		{
			match=1;
			double ldt=fabs(working[i].end_t-it);
			if(ldt<closest_dt)
			{
				closest_dt=ldt;
				closest_index=i;
			}
		}
	}
	
	
	/////here the problem is
	/// if there are two in the same z that should be added, the first will be added, 
	/// the second will be comparing to the one just added, not the one in the last plae
	
	
	for(unsigned int i=0;i<todel.size();i++)
	{
		working.erase(working.begin()+todel[i]-i); //need the extra i because items past deletion point will be shifted by 1 to the left for each previous deletion
	}
	
	
	
	
	
	//if no match found, start new chain
	if(!match)
	{
		Chain c;
		c.add_to_back(it,iz,ie,my_cluster_id);
		working.push_back(c);	
	}

}



void ChainHelper::AttachAsChild(Chain * parent, Chain * child)
{
	child->parentChain = parent->myId;
    parent->children.push_back(child->myId);	
	MSG("ChainHelper",Msg::kDebug)<<"attaching "<<parent->myId<<" to "<< child->myId  <<"setting child level from "<<child->level<<" to "<<1+parent->level <<endl;
	AdjustLevel(child, 1 + parent->level); //increase the level of children 

	found_max_path=0;
}


void ChainHelper::AdjustLevel(Chain * c, int level)
{
	if(!c)return;
	c->level=level;
	for(unsigned int i=0;i<c->children.size();i++)
	{
		Chain *d=GetChain(c->children[i]);
		if(d)
			AdjustLevel(d,level+1);
		else
			MSG("ChainHelper",Msg::kError)<<"bogus children list - parent id " << c->myId << " request for child with id "<< c->children[i]<<endl;
	}
} 




std::vector<int> ChainHelper::GetAllChildren(Chain *c){
	
	std::vector<int>ret;
	for(unsigned int i=0;i<c->children.size();i++)
	{
		Chain *me = GetChain(c->children[i]);
		if(!me)continue;
		std::vector<int> t= GetAllChildren(me);
		for(unsigned int j=0;j<t.size();j++)
		{
			ret.push_back(t[j]);
		}
	}

	ret.push_back(c->myId);
	return ret;

}


void ChainHelper::ChangeHitId(int oldid, int newid)
{
	for(unsigned int i=0;i<finished.size();i++)
	{
		Chain *c = &(finished[i]);
		for(unsigned int j=0;j<c->cluster_id.size();j++)
		{
			if(c->cluster_id[j]==oldid)c->cluster_id[j]=newid;
		}
	}
}



std::vector<foundpath> ChainHelper::FindPaths(Chain*c)
{
	std::vector<foundpath> a;
	if(!c)
	{
	//	printf("ChainHelper::FindPaths called without a chain!\n");
		return a;
	}
	//printf("finding paths for chain %d has %d children\n",c->myId,c->children.size());

	

	if(c==0)return a;

	if(c->children.size()<1)
	{
		foundpath b;
		b.end_t=c->end_t;
		b.end_z=c->end_z;
		b.energy=c->sum_e;
		b.muonlike=c->muonfrac;
		b.start_t=c->start_t;
		b.start_z=c->start_z;
		b.path.push_back(c->myId);
		b.entries=c->entries;
		a.push_back(b);
		return a;	
	}

	int foundsingle=0;

	for(unsigned int i=0;i<c->children.size();i++) 
	{

		std::vector<foundpath> d = FindPaths(GetChain(c->children[i]));
		//see 2+ of these paths are childless
		//if so, see if they have the same number of hits
		//if so, remove one of these for now....
		
		int issingle=0;
		
		if(d.size()==1)
			if(d[0].path.size()==1)
				if(GetChain(d[0].path[0])->entries<3)
					issingle=1;
		
		
		if(issingle && ! foundsingle) //only store one!
			a.push_back(d[0]);
		else if(!issingle)//store them all
			for(unsigned int j=0;j<d.size();j++)
				a.push_back(d[j]);
		if(issingle)foundsingle=1;
	}


	for(unsigned int i=0;i<a.size();i++)
	{
		a[i].start_t = c->start_t;
		a[i].start_z = c->start_z;
		a[i].energy+=c->sum_e;
		a[i].muonlike = (a[i].muonlike*a[i].entries + c->muonfrac*c->entries) / (a[i].entries+c->entries); //this may give unfair (2x) weight to duplicated head entries in children!
		a[i].entries+=c->entries -1; //head entry is duplicated in child!
		a[i].path.insert(a[i].path.begin(),c->myId);  //perhaps should replace with a queue!
	}

	return a;
}



