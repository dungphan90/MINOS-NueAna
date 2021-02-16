#include "MessageService/MsgService.h"

#include "NueAna/ParticlePID/ParticleFinder/Finder.h"
//#include "NueAna/ParticlePID/ParticleFinder/EMFitter.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle3D.h"

#include "NueAna/ParticlePID/ParticleFinder/ParticleType.h"
#include "NueAna/ParticlePID/ParticleFinder/LongMuonFinder.h"
#include "NueAna/ParticlePID/ParticleFinder/PrimaryShowerFinder.h"
#include "Conventions/Detector.h"

#include "TPolyLine.h"
#include "TMarker.h"
#include "TMath.h"
#include <algorithm>

#include <sstream>
#include <math.h>
#include "ShwFit.h"

#include <limits>





CVSID("$Id: Finder.cxx,v 1.10 2014/02/18 03:37:06 rhatcher Exp $");







Finder::Finder():  vtx_u(0.0), vtx_v(0.0), vtx_z(0.0)
{

        sorter_map.clear();
        DoMRCC=0;





        
        //cluster_map.clear();
        meupergev=0;
}

Finder::~Finder()
{

        sorter_map.clear();
        //cluster_map.clear();

}

void Finder::Reset()
{
        sorter_map.clear(); 


        loc_map.clear(); 
        plane.clear(); 
        strip.clear(); 
        energy.clear(); 
        particles.clear(); 
        t.clear();
        z.clear(); 
        view.clear();  
        //cluster_map.clear();

        true_vtx_u=0.0;
        true_vtx_v=0.0;
        true_vtx_z=0.0;


        vtx_u=0.0;
        vtx_v=0.0;
        vtx_z=0.0;
        shareholder.Reset();
}


void Finder::AddStrip(int myplane, int mystrip, double myenergy, double st, double sz, int myview)
{
        mypoh->event.nstrips++;


        mypoh->strips.AddStrip(myplane,mystrip,myenergy,st,sz,myview);

        plane.push_back(myplane); 
        strip.push_back(mystrip); 
        energy.push_back(myenergy); 
        t.push_back(st); 
        z.push_back(sz);
        view.push_back(myview);

        mypoh->event.large_minz = (mypoh->event.large_minz < sz) ? mypoh->event.large_minz : sz;
        mypoh->event.large_maxz = (mypoh->event.large_maxz > sz) ? mypoh->event.large_maxz : sz;
        if(myview==2)
        {
                mypoh->event.large_minu = (mypoh->event.large_minu < st) ? mypoh->event.large_minu : st;
                mypoh->event.large_maxu = (mypoh->event.large_maxu > st) ? mypoh->event.large_maxu : st;
        }else if(myview==3){
                mypoh->event.large_minv = (mypoh->event.large_minv < st) ? mypoh->event.large_minv : st;
                mypoh->event.large_maxv = (mypoh->event.large_maxv > st) ? mypoh->event.large_maxv : st;        
        }


        //printf("adding %f %f %d\n",st,sz,myview);
    sorter_map.insert(std::pair <double,int>(myenergy, (int)energy.size()-1))
;
        
    loc_map[myplane][mystrip]=plane.size()-1;




}


void Finder::ClearXTalk()
{
        printf("DONT USE Finder::ClearXTalk()!\n");
        exit(1);

        double thresh = 3; //number of mips below which we don't care about looking for xtalk....

        for(unsigned int i=0;i<plane.size();i++)
        {
                sorter_map.insert(std::pair <double,int>(energy[i],i));
                loc_map[plane[i]][strip[i]]=i;
        }
        

        
        double reassignedxtalke=0.0;
        
        int foundxtalk=0;
        std::map<int, std::map<int, int> >::iterator p_iter;
        for(p_iter=loc_map.begin();p_iter!=loc_map.end(); p_iter++)
    {
        std::map<int, int>::iterator s_iter;
        for(s_iter=p_iter->second.begin();s_iter!=p_iter->second.end(); s_iter++)
        {
                        if (energy[s_iter->second]<thresh)continue;
                        
                        std::map<int, int>::iterator s_iter2;
                        for(s_iter2=p_iter->second.begin();s_iter2!=p_iter->second.end(); s_iter2++)
                {
                                if(s_iter==s_iter2)continue;
                                int sp=abs(s_iter->first-s_iter2->first);
        //                      printf("sep %d view %d plane %d high %f low %f\n",sp,view[s_iter->second],plane[s_iter->second],energy[s_iter->second],energy[s_iter2->second]);
                                
                                if( sp>9 && sp<14)//sp == 13)
                                {
                                        if(energy[s_iter2->second] < 0.3 * energy[s_iter->second]) //suspected crosstalk hit is < 30% of the high hit
                                        {
                                                reassignedxtalke+=energy[s_iter2->second];
                                                
                                                energy[s_iter->second]+=energy[s_iter2->second];
                                                energy[s_iter2->second]=0.0;
                                                foundxtalk=1;
                                        }
                                }
                        
                        }
                }
        }


//      if(reassignedxtalke>0)printf("XTALK FILTER --- %f reassigned!\n",reassignedxtalke);


        if(foundxtalk)
        {
                std::vector<int>tplane;
                std::vector<int>tstrip;
                std::vector<int>tview;
                std::vector<double>tt;
                std::vector<double>tz;
                std::vector<double>tenergy;
        
        
                for(unsigned int i=0;i<energy.size();i++)
                {
                        if(energy[i]>0.001)
                        {
                                tplane.push_back(plane[i]); 
                                tstrip.push_back(strip[i]); 
                                tenergy.push_back(energy[i]); 
                                tt.push_back(t[i]); 
                                tz.push_back(z[i]);
                                tview.push_back(view[i]);
                        }
                }
                
                plane=tplane;
                strip=tstrip;
                view=tview;
                t=tt;
                z=tz;
                energy=tenergy;
                
        }


        sorter_map.clear();
        loc_map.clear();

}



void Finder::SetStatus(Chain *c, int status)
{

        for(unsigned int i=0;i<c->cluster_id.size();i++)
        {
                mypoh->clusterer.GetCluster(c->cluster_id[i])->SetStatus(status);
        //      printf("setting status cluster id %d to %d\n",c->cluster_id[i],status);
        }
}


void Finder::Process(ParticleObjectHolder & p)
{
        if(!meupergev)
        {
                cout <<"At Finder::Process... meupergev has not been set!\n";
                exit(1);
        }

        
        mypoh=&p;

//don't clear the xtalk here... use the XTalkFilter module!     
//      if(p.GetHeader().GetVldContext().GetDetector()==Detector::kFar)ClearXTalk(); //only far has xtalk
        

/*      
        for(unsigned int i=0;i<plane.size();i++)
        {
                sorter_map.insert(std::pair <double,int>(energy[i],i));
                loc_map[plane[i]][strip[i]]=i;
        
                Managed::ClusterManager * clusterhelper = (view[i]==2) ? &mypoh->clusterer : (view[i]==3) ? &mypoh->clusterer : 0;
                if(clusterhelper)clusterhelper->AddStrip(plane[i], strip[i], energy[i], t[i], z[i], view[i]);

        }
        
*/

        
//mypoh->clusterer=clustermanager_u;
mypoh->clusterer=clustermanager_v;
        
//      mypoh->clusterer.MakeClusters();
//      mypoh->clusterer.MakeClusters();

        mypoh->event.nclusters =  mypoh->clusterer.clusters.size();
//      mypoh->event.nclusters += mypoh->clusterer.clusters.size();
        

        std::vector<Managed::ManagedCluster> clusters = mypoh->clusterer.clusters;
/*      for(int i=0;i<clusters.size();i++)
        {
                printf("clu id %d view %d z %f t %f e %f\n",clusters[i].id,clusters[i].view,clusters[i].z,clusters[i].t,clusters[i].e);
        
        
        }
*/
        //currently required for the display... should clean it up later!
        mypoh->clusterer.FillClusterMap(&mypoh->cluster_map);
        mypoh->clusterer.FillClusterMap(&mypoh->cluster_map);
        

        //look for long muons
        mypoh->clusterer.MakeClusters(0.02,0.05,0.01);

        LongMuonFinder lmf(mypoh->GetHeader().GetVldContext().GetDetector());
        int foundlongmuon = lmf.FindLongMuon(&mypoh->clusterer);

        int didMRCC=0;

        if(foundlongmuon )
        {
                        mypoh->eventquality.foundlongmuon=1;
                                Particle3D *lp=lmf.GetParticle3D();
                                
                                if(lp)
                                {
                                        std::vector<Particle3D*> pout1;
                                        pout1.push_back(lp);    
                
                                        FinalizeParticles3D(pout1); //final computation of values, etc...
                        
                                        if(!DoMRCC)
                                        {
                                                for(unsigned int i=0;i<pout1.size();i++)
                                                {
                                                        mypoh->AddParticle3D(*(pout1[i]));
                                                }
                                        }else{
                                                MRCCInfo * mrccinfo = mypoh->GetMRCCInfo();
                                                mrccinfo->removedmuon = *(pout1[0]); 
                                                mrccinfo->stage=1;
                                                didMRCC=1;
                                                
                                        }               
                                }
                                
                                //we need to store the chains to prevent further use...even if doing MRCC!
                                
                                int cu = mypoh->chu.NewChain();
                                Chain * tempcu = mypoh->chu.GetChain(cu);
                                *tempcu = *lmf.GetChain(2);

                                int cv = mypoh->chv.NewChain();
                                Chain * tempcv = mypoh->chv.GetChain(cv);
                                *tempcv = *lmf.GetChain(3);
                                
                                if(DoMRCC)
                                {
                                        SetStatus(tempcu,-10);
                                        SetStatus(tempcv,-10);
                                }
                        
        }


        //cout <<"*********Currently have "<<mypoh->particles3d.size()<<" particles in POH\n";
        if(DoMRCC)foundlongmuon=0;//we don't want to use the muon vertex in the reco!

        if(lmf.FoundSingleViewLongMuon())
        {
                mypoh->eventquality.single_view_long_muon=1;
                //printf("!!!! Found a single good muon in only one view!\n");
        
        }


        //if we don't have a long muon, then this is probably either a NC or nue event
        //so lets look for a primary shower to get a fix on the event vertex
        mypoh->clusterer.MakeClusters(0.2,0.05,0.05);


        mypoh->clusterer.DumpClusters();

        int foundprimaryshower =0;

        PrimaryShowerFinder*psf = &mypoh->psf;  //mypoh has the psf for the display...
        psf->chu=&mypoh->chu;
        psf->chv=&mypoh->chv;

        //if foundlongmuon
        if(foundlongmuon)
        {
                Particle3D *lp=lmf.GetParticle3D();
                //printf("FOUND LONG MUON setting psf vertex %f %f %f\n",lp->start_u, lp->start_v, lp->start_z);


                psf->SetVertex(lp->start_u, lp->start_v, lp->start_z);
        }

        foundprimaryshower = psf->FindPrimaryShower(&mypoh->clusterer);


        if(foundprimaryshower)
        {
                mypoh->eventquality.foundprimaryshower=1;
        
                                Particle3D *lp=psf->GetParticle3D();
                                
                                if(lp)
                                {
                                        std::vector<Particle3D*> pout1;
                                        pout1.push_back(lp);    
                
                                        FinalizeParticles3D(pout1); //final computation of values, etc...

                
                                        if(!DoMRCC || (pout1.size()>0 && pout1[0]->particletype != Particle3D::muon) || didMRCC)
                                        {
                                                for(unsigned int i=0;i<pout1.size();i++)
                                                {
                                                        mypoh->AddParticle3D(*(pout1[i]));
                                                }
                                        }else{
                                                MRCCInfo * mrccinfo = mypoh->GetMRCCInfo();
                                                mrccinfo->removedmuon = *(pout1[0]); 
                                                mrccinfo->stage=2;
                                                
                                        }               

                                }
                
                //chains are being added in psf now...          
                /*              
                                int cu = mypoh->chu.NewChain();
                                Chain * tempcu = mypoh->chu.GetChain(cu);
                                *tempcu = *psf->GetChain(2);

                                int cv = mypoh->chv.NewChain();
                                Chain * tempcv = mypoh->chv.GetChain(cv);
                                *tempcv = *psf->GetChain(3);
                */              
                        
        }


        if(psf->FoundSingleViewPrimaryShower())
        {
                mypoh->eventquality.single_view_primary_shower=1;
                //printf("!!!! Found a single primary shower in only one view!\n");
        
        }




        //if foundprimaryshower
        if(foundprimaryshower)
        {
                Particle3D *lp=psf->GetParticle3D();
                vtx_u = lp->start_u;
                vtx_v = lp->start_v;
                vtx_z = lp->start_z;
        }       
        
        //if foundlongmuon
        if(foundlongmuon)
        {
                Particle3D *lp=lmf.GetParticle3D();
                
                if(foundprimaryshower)
                {
                        Particle3D *lpe=psf->GetParticle3D();
                        if(lpe->start_z > lp->start_z)
                        {
                                vtx_u = lp->start_u;
                                vtx_v = lp->start_v;
                                vtx_z = lp->start_z;
                        }
                }else{
                                vtx_u = lp->start_u;
                                vtx_v = lp->start_v;
                                vtx_z = lp->start_z;            
                
                }
        }       
        

        //cout <<"*********Currently have "<<mypoh->particles3d.size()<<" particles in POH\n";

//mypoh->clusterer.GetClusterSaver()->DumpClusters();

        //now look in clusters for tracks
        //remember to recluster!
        //!!!! do not remake clusters unless we are starting the chains from scratch!!!!
//      mypoh->clusterer.MakeClusters(0.2,0.05,0.05);


        MSG("Finder",Msg::kInfo) <<"----U----"<<endl;   


        
        ChainHelper * chu= &mypoh->chu;
        //Managed::ClusterManager  * clu = &mypoh->clusterer;
/*      chu->max_z_gap=0.3;
        chu->max_slope=3;
        chu->SetDetector(mypoh->GetHeader().GetVldContext().GetDetector());
        MakeChains(clu,chu,2);
*/      
        MSG("Finder",Msg::kInfo) <<"----V----"<<endl;   
        
        ChainHelper * chv= &mypoh->chv;
        //Managed::ClusterManager  * clv = &mypoh->clusterer;
/*      chv->max_z_gap=0.3;
        chv->max_slope=3;
        chv->SetDetector(mypoh->GetHeader().GetVldContext().GetDetector());
        MakeChains(clv,chv,3);
*/

 
        //first guess of vertex 
//      FindVertex(chu,chv);
//      MSG("Finder",Msg::kInfo) << "first guess - vertex at (u,v,z)  ("<< vtx_u <<", "<< vtx_v <<", "<< vtx_z <<")"<<endl;
//      MSG("Finder",Msg::kInfo) << "first guess - true vertex at (u,v,z)  ("<<true_vtx_u <<", "<< true_vtx_v <<", "<< true_vtx_z <<")"<<endl;

        chu->print_finished();
        chv->print_finished();



        MergeChains(chu);
        MergeChains(chv);

        


        FindVertex(chu,chv);
        
        
        //if foundprimaryshower
        if(foundprimaryshower)
        {
                Particle3D *lp=psf->GetParticle3D();
                vtx_u = lp->start_u;
                vtx_v = lp->start_v;
                vtx_z = lp->start_z;
                
        }       
        
        //if foundlongmuon
        if(foundlongmuon)
        {
                Particle3D *lp=lmf.GetParticle3D();
                
                if(foundprimaryshower)
                {
                        Particle3D *lpe=psf->GetParticle3D();
                        if(lpe->start_z > lp->start_z)
                        {
                                vtx_u = lp->start_u;
                                vtx_v = lp->start_v;
                                vtx_z = lp->start_z;
                        }
                }else{
                
                        vtx_u = lp->start_u;
                        vtx_v = lp->start_v;
                        vtx_z = lp->start_z;            
                }
        }       
                
        
        MSG("Finder",Msg::kInfo) << "vertex at (u,v,z)  ("<< vtx_u <<", "<< vtx_v <<", "<< vtx_z <<")"<<endl;
        MSG("Finder",Msg::kInfo) << "true vertex at (u,v,z)  ("<<true_vtx_u <<", "<< true_vtx_v <<", "<< true_vtx_z <<")"<<endl;




//      RemoveNonVertexPointingChains(chu,2);
//      RemoveNonVertexPointingChains(chv,3);   



        chu->print_finished();
        chv->print_finished();

//      mypoh->clusterer.GetClusterSaver()->DumpClusters();

        //cout <<"*********Currently have "<<mypoh->particles3d.size()<<" particles in POH\n";
        Weave(chu,chv);
        //cout <<"*********Currently have "<<mypoh->particles3d.size()<<" particles in POH\n";


        std::vector<Particle3D * > p3d;
/*      for(int i=0;i<mypoh->particles3d1->GetEntries();i++)
                p3d.push_back( (Particle3D * )mypoh->particles3d1->At(i) );
*/

        p3d.clear();
        for(unsigned int i=0;i<mypoh->particles3d.size();i++)
                p3d.push_back(&mypoh->particles3d[i]);
        
        DumpParticle3D(p3d);
        //cout <<"*********Currently have "<<mypoh->particles3d.size()<<" particles in POH\n";
                                
/*      for (unsigned int i=0;i<p3d.size();i++)
        {
                //ostringstream s;
        
                Particle3D * p = p3d[i];
                if (p==0)continue;
        
//              printf("saving with parb %f\n",p->emfit_b);
        }
*/


        //and store event params now...
        mypoh->cluster_saver.recomputeBounds();
        mypoh->event.minu=mypoh->cluster_saver.minu;
        mypoh->event.maxu=mypoh->cluster_saver.maxu;
        mypoh->event.minv=mypoh->cluster_saver.minv;
        mypoh->event.maxv=mypoh->cluster_saver.maxv;

        mypoh->event.minz=mypoh->cluster_saver.minz;
        mypoh->event.maxz=mypoh->cluster_saver.maxz;

        //recalculate the number of strips and event energy

        std::map<std::pair<int,int>, double> stripmap =mypoh->cluster_saver.GetStripEnergy();
        std::map<std::pair<int,int>, double>::iterator it;
        
        int stripcount=0;
        double  stripe=0;

        for(it = stripmap.begin();it!=stripmap.end();it++)
        {
                if(it->second<1e-6)continue;
                stripcount++;
                stripe+=it->second;
//              printf("%d %d %f\n",it->first.first,it->first.second,it->second);
        }

        double olde = mypoh->event.visenergy;
        mypoh->event.nstrips=stripcount;
        mypoh->event.visenergy=stripe;
        mypoh->event.nclusters =  mypoh->cluster_saver.nClusters;


//      printf("old e %f new e %f strips %d\n",olde,stripe,stripcount);


        if(stripe - olde > 1)
        {
                MSG("Finder",Msg::kError) << "Event Energy MISMATCH! Strip SumE = "<<olde<<" new ClusterSaver Energy "<<stripe<<"\n";   

        }
}


//this will attempt to merge chains which point to eachother
void Finder::MergeChains(ChainHelper *ch)
{

        std::vector<std::pair<int,int> >match;

        //iter over head chains
        for(int i=0;i<(int)ch->finished.size();i++)
        {
                if(ch->finished[i].parentChain>-1)continue; //make it a head chain
                if(ch->finished[i].entries<2)continue;
                
        
                double bestdist=10000;
                //double maxtdist=0.1;
                int bestid=-1;
                double bestzdist=10000;
                
                //iter over ends of chains
                for(int j=0;j<(int)ch->finished.size();j++)
                {
                        if(i==j)continue;
                        if(ch->finished[j].children.size()>0)continue; //make it and ending chain
                        if(ch->finished[j].entries<2)continue;
                        
                        Chain *e  = &ch->finished[j]; //end of one chain
                        Chain *b = &ch->finished[i];  //beginning of another chain
                        
                        if(b->start_z-e->end_z < 0.035*2)continue; // require at least 2 planes away...
                //      if(e->start_z - b->start_z< 0.035*2)continue; // require at least 2 planes away...
                        
                        double d = e->end_t - (e->end_z * b->front_slope + b->front_offset);
                
                //      double d = maxtdist;
                        
                        //need to make sure that the chain to be attached is vertex pointing...
                        

/*
//////////////////////                  
                        double canslope=b->front_slope;
                        double canoffset=b->front_offset;
                        double vpslope=e->back_slope;
                        double vpoffset=e->back_offset;
                        
                        
                        double t = ( vpoffset * canslope - vpslope *canoffset ) / (canslope - vpslope);
                        double z = (canoffset - vpoffset) / (vpslope - canslope);                       
                        
                        
                        printf("vp %d a %f b %f  can %d a %f b %f ------ t %f z %f\n",e->myId,e->a, e->b,b->myId,b->a, b->b,t,z);
                        printf("D t %f z%f \n",t,z);
                                
                        if( z > e->end_z && z < b->start_z) // its a possible match! --- see if its the closest match...
                        if( t > e->end_t && t < b->start_t ||  t < e->end_t && t > b->start_t) //make sure its not kinky
                        {
                        
                                printf("found match at z %f   t %f\n",z,t);
                                
                                d = e->end_t - (e->end_z * b->front_slope + b->front_offset);
                        
                        }
//////////////////////
*/              
                        
                        
//                      if(d<maxtdist)
//                      {
                                double totdist = sqrt( d*d + (e->end_z-b->start_z)*(e->end_z-b->start_z));
                        
                                double zdist  = fabs(e->end_z-b->start_z);
                                if(totdist<bestdist || zdist < bestzdist)
                                {
                                        int points=0;
                                        double canslope=b->front_slope;
                                        double canoffset=b->front_offset;
                                        double vpslope=e->back_slope;
                                        double vpoffset=e->back_offset;
                                        double t = 0;
                                        double z = 0;

                                        if(canslope-vpslope>1e-10)
                                        {
                                                t = ( vpoffset * canslope - vpslope *canoffset ) / (canslope - vpslope);
                                                z = (canoffset - vpoffset) / (vpslope - canslope);
                                        }else{
                                                continue;
                                        }

                                //      if( z > e->end_z && z < b->start_z) points=1;
                                        
                                        if(fabs(e->back_slope - b->front_slope)<0.5)points=1;
                                
                        //              cout <<"back slope "<< e->back_slope<<" front slope "<<b->front_slope<<endl;    
                                        
                                        //between supermodules in far?
                                        double maxtotdist=0.5;
                                        if(mypoh->GetHeader().GetVldContext().GetDetector() == Detector::kFar && fabs(e->end_z - 14.6)<0.2 && fabs(b->start_z-16)<0.2)maxtotdist+=16-14.6;
                                
                                        //in back of near (partially instrumented every fourth plane)?
                                        if(mypoh->GetHeader().GetVldContext().GetDetector() == Detector::kNear && e->end_z >6 && b->start_z>6)maxtotdist*=4.;
                                
                                
                                        //make sure that the connection would not be kinky
                                        if(points)
                                        {
                                                if((e->end_z-b->start_z)==0)continue;// that would be kinky anyways...
                                                double connecting_slope = (e->end_t-b->start_t) / (e->end_z-b->start_z);
                                                
                                                if(fabs(connecting_slope - b->front_slope)>=0.5)points=0;
                                                if(fabs(e->back_slope - connecting_slope)>=0.5)points=0;                                        
                                        }
                                
                                        
                                        
                                        if(!points /*&& totdist>maxtotdist*/)
                                        {
                                                MSG("Finder",Msg::kDebug) << "failing on points "<<e->myId <<" to "<< b->myId<<" totdist "<< totdist<<" zdist "<<zdist <<" tdist "<<d <<" intersection at z t"<<z<<" "<<t<<"\n";
                                                continue;
                                        }
                                        
                                        
                                        
                                        bestdist=totdist;
                                        bestid=e->myId;
                                        bestzdist=zdist;
                                        
                                        MSG("Finder",Msg::kDebug) << "found better match "<<e->myId <<" to "<< b->myId<<" totdist "<< totdist<<" zdist "<<zdist <<" tdist "<<d <<"\n";
                                        
                                }
                        }       
                        
//              }
                
                
                if(bestid>-1)
                {
                        match.push_back(std::pair<int,int>(bestid,ch->finished[i].myId));
                
                        MSG("Finder",Msg::kDebug) <<"!!!matched "<<bestid<<" "<<ch->finished[i].myId<<"\n";
                
                }
        }
        
        
        for(int i=0;i<(int)match.size();i++)
        {
                ch->AttachAsChild(ch->GetChain(match[i].first),ch->GetChain(match[i].second));
        
        }



}





//after finding a vertex, run this function to remove chains that either do not point to the vertex  
// or do not point to another chain that points to the vertex...
// if the chain points to another chain that points to the vertex, then connect them!
//
//  this will use the chains slope and offset... but if there are multiple chains, with the first few only having a few hits
//  then a new slope and offset should be constructed....   (not for now)
//
void Finder::RemoveNonVertexPointingChains(ChainHelper *ch, int view)
{

        double max_dt = 0.1; // max distance of chain pointing to vertex z from vertex in t

        double vt =0.0;
        if(view==2) vt=vtx_u;
        if(view==3) vt=vtx_v;


        std::vector<int> marked;
        std::vector<int> vtxpointing;
        std::vector<int> head;



        for(unsigned int i=0;i<ch->finished.size();i++)
        {       
                ostringstream os;
                
                //here is where we can reverse chains if they point the wrong way
                //chain begin should be closer to the vertex than the end....
                
                
                
                
                ///     
/*              
                os << "Chain id " << ch->finished[i].myId << " entries: " << ch->finished[i].entries<< " extrap posvtx: " <<  ch->finished[i].avg_slope * vtx_z + ch->finished[i].avg_offset << " vtxt " << vt << "  diff " <<fabs(vt - (ch->finished[i].avg_slope * vtx_z + ch->finished[i].avg_offset) ) ;
        
                if(ch->finished[i].entries>1)
                if(ch->finished[i].parentChain<0 || (ch->finished[i].level==1 && ch->GetChain(ch->finished[i].parentChain)->entries==1) ) //only look at head chains...         
                if(max_dt < fabs(vt - (ch->finished[i].avg_slope * vtx_z + ch->finished[i].avg_offset) ) )
                {
                        //its too far away...
                        //throw it out
                        marked.push_back(ch->finished[i].myId);
                        
                        os << "not pointing to vertex ";
                        
                }else{
                        vtxpointing.push_back(ch->finished[i].myId);
                }
                MSG("Finder",Msg::kInfo) << os.str() << endl; 
*/


                int close_enough=0;
                double diff = fabs(vt - (ch->finished[i].front_slope * vtx_z + ch->finished[i].front_offset) );
                close_enough = max_dt > diff;
                
os << "Chain id " << ch->finished[i].myId << " entries: " << ch->finished[i].entries<< " extrap posvtx: " <<  ch->finished[i].front_slope * vtx_z + ch->finished[i].front_offset << " vtxt " << vt<<" diff "<<diff<<endl;
                                
                //using interpolation is probably just better for finding hits to attach...
                diff = fabs(vt - ch->finished[i].interpolate(vtx_z) );
                if(!close_enough)close_enough = max_dt > diff;

                os << "Chain id " << ch->finished[i].myId << " entries: " << ch->finished[i].entries<< " extrap posvtx: " <<  ch->finished[i].interpolate(vtx_z)  << " vtxt " << vt<<" diff "<<diff<<endl;              

                //using interpolation is probably just better for finding hits to attach...
                diff = fabs(vt - (ch->finished[i].avg_slope * vtx_z + ch->finished[i].avg_offset) );
                if(!close_enough)close_enough = max_dt > diff;
                
                


                os << "Chain id " << ch->finished[i].myId << " entries: " << ch->finished[i].entries<< " extrap posvtx: " <<  (ch->finished[i].avg_slope * vtx_z + ch->finished[i].avg_offset)  << " vtxt " << vt <<" diff "<<diff;

                /*        
                if(ch->finished[i].entries>1)
                if(ch->finished[i].parentChain<0 || (ch->finished[i].level==1 && ch->GetChain(ch->finished[i].parentChain)->entries==1) ) //only look at head chains...         
                if(!close_enough)
                {
                        //its too far away...
                        //throw it out
                        marked.push_back(ch->finished[i].myId);
                        
                        os << " not pointing to vertex ";
                        
                }else{
                        vtxpointing.push_back(ch->finished[i].myId);
                }
                */
                // RWH 2014-02-17 Above is what used to be here ...
                // but the compiler is completely confused about dangling else 
                // conditions and can't decide which to attach to which if
                // this is current best guess of what (might have been)/was intended
                if (ch->finished[i].entries>1) {
                  if (ch->finished[i].parentChain<0 || 
                      (ch->finished[i].level==1 && ch->GetChain(ch->finished[i].parentChain)->entries==1) ) {
                    //only look at head chains...         
                    if (!close_enough) {
                      //its too far away...
                      //throw it out
                      marked.push_back(ch->finished[i].myId);
                      
                      os << " not pointing to vertex ";
                      
                    } else {
                      vtxpointing.push_back(ch->finished[i].myId);
                    }
                  }
                }
                // RWH end of problematic area

                MSG("Finder",Msg::kInfo) << os.str() << endl;


        }
        
        
        std::vector<int> todel;
        
        for (unsigned int i=0;i<marked.size();i++)
        {
                //see if the chain is pointing to a vertexpointing chain....  if so, attach it...
                //see if the projections cross somewhere between the start of the vertexpointing chain and the start of the candidate chain
                
                Chain * can = ch->GetChain(marked[i]);
                
                if(can->entries<2)continue;
                
                
                double bestdist = 100000;
                int bestpid = -1;
                double attachpoint=0;
                
                for(unsigned int j=0;j<vtxpointing.size();j++)
                {
                        //printf("looking to see if it points to a good chain for attachment....\n");
                
                        Chain *vp = ch->GetChain(vtxpointing[j]);
                        //printf("A\n");        
                        if(vp->entries<2)continue;
                
                        //printf("B\n");        
                        if(fabs(vp->a - can->a) < 0.000001)continue;  //dont want a fpe
                        //printf("C\n");        
                        if(can->start_z < vp->start_z)continue;
                        
        //              double t = ( vp->b * can->a - vp->a * can->b ) / (can->a - vp->a);
        //              double z = (can->b - vp->b) / (vp->a - can->a);
        
                
//                      double t = ( vp->a * can->b - vp->b * can->a ) / (can->b - vp->b);
//                      double z = (can->a - vp->a) / (vp->b - can->b);
                

                        double canslope=can->front_slope;
                        double canoffset=can->front_offset;
                        double vpslope=vp->front_slope;
                        double vpoffset=vp->front_offset;
                        
                        
                        double t = ( vpoffset * canslope - vpslope *canoffset ) / (canslope - vpslope);
                        double z = (canoffset - vpoffset) / (vpslope - canslope);                       
                        
                        
                        //printf("vp a %f b %f  can a %f b %f ------ t %f z %f\n",vp->a, vp->b,can->a, can->b,t,z);
                        
                        
                        
                //      z+=vp->start_z;
                        //printf("D t %f z%f \n",t,z);  
                        if( z > vp->start_z && z < can->start_z) // its a possible match! --- see if its the closest match...
                        {
                        
                                //printf("found match at z %f   t %f\n",z,t);
                                
                                double bd = bestdist+1  ;
                        
                                attachpoint =z;
                        
                                if( z < vp->end_z) // attach to middle of chain...
                                {
                                        //printf("!!!! need to attach to middle of chain!!!\n");
                                        bd = sqrt ( (can->start_z-z)*(can->start_z-z) + (can->start_t - t)*(can->start_t - t));
                                        bestpid = vp->myId;
                                        
                                }else{
                                        //printf("attach to end of chain\n");
                                        bd = sqrt ( (vp->end_z-can->start_z)*(vp->end_z-can->start_z) + (vp->end_t - can->start_t)*(vp->end_t - can->start_t));
                                }
                                if(bd<bestdist)
                                {
                                        bestdist=bd;
                                        bestpid = vp->myId;
                                }
                        }
                }
                
                
        
        
                //if not, delete it
                if(bestpid>-1)
                {
                        if(attachpoint < ch->GetChain(bestpid)->end_z) //need to break apart vertex pointing chain
                        {
                                Chain * np = ch->SplitAt(ch->GetChain(bestpid),attachpoint);  //returns the chain that now ends at attachpoint
                
                                ch->AttachAsChild(np,can);
                //              can->parentChain = np->myId;
                //              np->children.push_back(can->myId);
                        
                        }else{   //just attach to the end
                        
                                ch->AttachAsChild(ch->GetChain(bestpid),can);
                //              can->parentChain = bestpid;
                //              ch->GetChain(bestpid)->children.push_back(can->myId);
                        }
                }else{
                        todel.push_back(marked[i]);
                }
        }
        


        
        
        for(unsigned int i=0;i<todel.size();i++)
                ch->DeleteChain(todel[i]);
                
        
        
        
        ostringstream os;
        os << "remaing chains: ";
        for(unsigned int i=0;i<ch->finished.size();i++)
        {
                os << ch->finished[i].myId<< " ";
        }
        
        MSG("Finder",Msg::kInfo) << os.str() << endl; 
}





void Finder::FindMuons(ChainHelper *ch)
{


//travel over chain - if we have at least two muon like hits, make a new muon chain going out to last muon hit
        //this new chain starts at first hit in the chain, and is not part of the primary chain
        //all muon like hits have all energy moved to new chain
        //larger hits have some average muon energy removed!
        
        

        std::vector<int> path = ch->FindMaxPath();
        if(path.size()<1)return;


        std::vector< std::pair<int,int> > muonlike;
        std::vector<double> muonlikee;
        
        int lastconsec=0;
        
        double avgmuon=0;
        int amuon=0;
        
        for(int i=0;i<(int)path.size();i++)
        {
        
                MSG("Finder",Msg::kDebug) <<"path chain "<<path[i]<<endl;
                int head=(i>0)?1:0;
                for(int j=head;j<ch->GetChain(path[i])->entries;j++)
                {
                        double e = ch->GetChain(path[i])->e[j];
                        if(e<2 && e>0.5){
                                if((i!=(int)path.size()-1 || j!=ch->GetChain(path[i])->entries-1) || amuon>0)  
                                //don't add if the only muonlike hit is at the end of a long chain
                                //it could be a low energy shower hit!
                                {
                                        muonlike.push_back(std::make_pair(path[i],j));
                                        muonlikee.push_back(e);
                                        avgmuon+=e;
                                        amuon++;
                                        lastconsec++;
                                }
                        
                        }else{
                        
                          /*
                                if(i==(int)path.size()-1 && j>0)
                                if(e<6.0*ch->GetChain(path[i])->e[j-1])  //try not to take overlapping large hits
                                {
                                        if(lastconsec>0)
                                        {
                                                muonlike.push_back(std::make_pair(path[i],j));
                                                muonlikee.push_back(e);
                                        }
                                }else{
                                        lastconsec=0;
                                }       
                          */
                          // RWH 2014-02-17 Above is what used to be here ...
                          // but the compiler is completely confused about dangling else 
                          // conditions and can't decide which to attach to which if
                          // this is current best guess of what (might have been)/was intended
                          if (i==(int)path.size()-1 && j>0) {
                            if (e<6.0*ch->GetChain(path[i])->e[j-1]) {
                              //try not to take overlapping large hits
                              
                              if (lastconsec>0) {
                                muonlike.push_back(std::make_pair(path[i],j));
                                muonlikee.push_back(e);
                              }
                            } else {
                              lastconsec=0;
                            }       
                          }
                          // RWH end of problematic area

                        }
                }
        }
        
        if(amuon==0)return; //no muon hits found!
        
        
        avgmuon/=amuon;
        
//      double avgmuon;
//      for(int i=0;i<muonlikee.size();i++)
//              avgmuon+=muonlikee[i];
//      avgmuon/=muonlikee.size();
        
        //subtract out found muon hits
//      for(int i=0;i<muonlike.size();i++)
//              ch->GetChain(muonlike[i].first)->e[muonlike[i].second]=0.0;
        
        
        //subtract out avgmuon from all hits greater than avgmuon
        //make new chain, with found levels, and avglevels for rest
        int cid = ch->NewChain();
        Chain * c = ch->GetChain(cid);



//      Chain c;

        for(unsigned int i=0;i<path.size();i++)
        {
                int head=i>0?1:0;
                
                Chain * d = ch->GetChain(path[i]);
                
                if(d==0)
                {
                        MSG("Finder",Msg::kError) <<"Bad path chain used - chain id "<<path[i]<<endl;
                        continue;
                
                }
                
                for(int j=head;j<d->entries;j++)
                {
        //              double e = 0;
/*              
                        if(d->e[j]>avgmuon)
                        {       d->e[j]-=avgmuon;
                                e=avgmuon;
                        }else{
                                e=d->e[j];
                                
                        
                        
                        
                //      double e = d->e[j];
                        for(int k=0;k<muonlike.size();k++)
                        {
                                if(muonlike[k].first==path[i] && muonlike[k].second==j)
                                {
                                        e= e > muonlikee[k] ? muonlikee[k] : e;
                                }
                        }
                                d->e[j]-=e;
                        }
*/

                        
                        //if(d->e[j]>avgmuon)d->e[j]-=avgmuon;
                        //double e = d->e[j];
                        
                        int found =0;
                        double e=0;
                        for(unsigned int k=0;k<muonlike.size();k++)
                        {
                                if(muonlike[k].first==path[i] && muonlike[k].second==j)
                                {
                                        e= muonlikee[k];
                                        found=1;
                                        break;
                                }
                        }
                        if(found) //we already taged it as muon like - so move the hit
                        {
                                //e=d->e[j];
                                d->e[j]=0;
                        }else{
                                e=avgmuon;
                                d->e[j]-=avgmuon;
                        }


                
                        if(e>0)
                        {       c->add_to_back(d->t[j],d->z[j],e,-1);
                                MSG("Finder",Msg::kDebug) <<"muon adding "<< path[i] <<", "<< j<<" energy "<<e <<endl;
                        }
                }
                
                
                d->Recalc(); //remove the hits that we set e to 0 from the chain 
        }
        
        
//print out chain ids....

//      int cid = ch->NewChain();
  //  Chain * e =ch->GetChain(cid);
  //  *e = c;
        
        
        //attach new chain to current head....
// //   c.parentChain=ch->GetChain(path[0]).myId;

        //ch->AttachAsChild(ch->GetChain(path[0]),c);
        
        c->particletype=Particle3D::muon; //say that we think it is a muon
        
        //add new chain to finished
//      Chain d = *c;
        ch->particles.push_back(*c);    
        delete c;       

}



void Finder::MakeChains(Managed::ClusterManager  *cl, ChainHelper * ch, int view)
{

        if(mypoh->event.vtx_z>0)
        {
                ch->vtx_z=mypoh->event.vtx_z;
                ch->vtx_t=view==2?mypoh->event.vtx_u:mypoh->event.vtx_v;
                
        }


     std::map<double, std::map<double, int> >::iterator p_iterr;
     std::map<double, int >::iterator s_iterr;

     double lastz=-1;
     for(p_iterr=cl->GetClusterMap(view)->begin();p_iterr!=cl->GetClusterMap(view)->end(); p_iterr++)
     {
        //for the first one...
                if(lastz==-1)lastz=p_iterr->first;
                
                
                if(fabs(p_iterr->first - lastz) > 0.01)
                {
                        ch->process_plane();
                        lastz=p_iterr->first;
                        //printf("process plane\n");
                }
      

         for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
         {

      
                Managed::ManagedCluster * c = cl->GetCluster(s_iterr->second);
                if(!c)
                {
                        //printf("!!! bad cluster map! finder.cxx around 897\n");
                        continue;
                }
                //printf("t %f z %f e %f id %d\n",s_iterr->first, p_iterr->first, c->e,c->id);
                ch->add_to_plane(s_iterr->first, p_iterr->first, c->e,c->id);   
             }
// ch->process_plane();
        
     }
     
     ch->process_plane();
 
        ch->finish();
        ch->print_finished();
        ch->FindMaxPath();

}


void Finder::FindVertex(ChainHelper * cu, ChainHelper * cv)
{
        std::vector<int> mpu = cu->FindMaxPath();
        std::vector<int> mpv = cv->FindMaxPath();
        
        if(mpu.size() <1 || mpv.size() <1)return;
        
        
        double cu_t = cu->GetChain(mpu[0])->start_t;
        double cu_z = cu->GetChain(mpu[0])->start_z;

        double cv_t = cv->GetChain(mpv[0])->start_t;
        double cv_z = cv->GetChain(mpv[0])->start_z;    
        
        
        ///////////////
        //try iterating through, measuring dt between adjacent hits in chain.
        //go through until dt(0,1)-dt(1,2)<thresh
        
        
/*      int forcedadvancedu=0;
        int forcedadvancedv=0;
*/
        

MSG("Finder",Msg::kDebug)<< "cuz " << cu_z <<" cvz "<< cv_z <<endl;
        
        
        ////try extrapolating vertex of all......
/*      
        if((forcedadvancedu && forcedadvancedv) || (!forcedadvancedu && !forcedadvancedv))
        {
        
                if(cv_z<cu_z)
                        cu_z=cv_z;
                else
                        cv_z=cu_z;
        }else{
                if(forcedadvancedu)
                        cv_z=cu_z;
                if(forcedadvancedv)
                        cu_z=cv_z;
        }
*/              



                if(cv_z<cu_z)
                        cu_z=cv_z;
                else
                        cv_z=cu_z;
                        

                int chuid=0;
                
        //      if(cu->GetChain(mpu[0]).start_z<cu_z)
                for(unsigned int i=0;i<mpu.size();i++)
                {
                        if ( ( cu->GetChain(mpu[i])->start_z <= cu_z && cu->GetChain(mpu[i])->end_z >= cu_z )
                                 || (cu->GetChain(mpu[i])->start_z >= cu_z && cu->GetChain(mpu[i])->end_z <= cu_z))
                        {
                                chuid=i;
                                
                        }
                }
                
                if(cu->GetChain(mpu[0])->entries<2 && mpu.size()>1)chuid=1;
                
                
                int chvid=0;
                
        //      if(cv->GetChain(mpv[0]).start_z<cv_z)
                for(unsigned int i=0;i<mpv.size();i++)
                {
                        if ( ( cv->GetChain(mpv[i])->start_z <= cv_z && cv->GetChain(mpv[i])->end_z >= cv_z )
                                 || (cv->GetChain(mpv[i])->start_z >= cv_z && cv->GetChain(mpv[i])->end_z <= cv_z))
                        {
                                chvid=i;
                                
                        }
                }
        
                if(cv->GetChain(mpv[0])->entries<2 && mpv.size()>1)chvid=1;
                //printf("chain has %d entries, path has %d entries\n",cv->GetChain(mpv[0])->entries,mpv.size());

                /*
                printf("uchain\n");
                for(unsigned int i=0;i<mpu.size();i++)printf("%d\n",mpu[i]);
                printf("vchain\n");
                for(unsigned int i=0;i<mpv.size();i++)printf("%d\n",mpv[i]);
        */      
                
                        Chain * uuu = cu->GetChain(mpu[chuid]);
                        Chain * vvv =  cv->GetChain(mpv[chvid]);
                
                
                
                        cu_z = uuu->start_z;
                        cv_z = vvv->start_z;
                
                        if(cv_z<cu_z)
                        {
                                cu_z=cv_z;
                        }
                        else
                        {
                                cv_z=cu_z;

                        }
                
                
                        double tmpuz=0;
                        double tmpvz=0;
                        
                //      double slopeu = cu->GetChain(mpu[chuid])->avg_slope;
                //      double bzu = cu->GetChain(mpu[chuid])->start_z;
                //      double btu = cu->GetChain(mpu[chuid])->start_t;
                        //double bt = cu->GetChain(mpu[chuid])->weighted_t;
                        
                                        //see if there is a big energy difference between the first 2 hits... if so, we probably want the 2nd hit
                                        if(uuu->entries>1 && uuu->e[0] < uuu->e[1]*0.2 )
                                        {
                                                //bz = uuu->z[1];
                                                //bt = uuu->t[1];
                                                tmpuz=uuu->z[1];
                                                cu_z=uuu->z[1];
                                        }
                                        else if(uuu->entries>2 && uuu->e[1] < uuu->e[2]*0.2 )
                                        {
                                                //bz = uuu->z[2];
                                                //bt = uuu->t[2];
                                                tmpuz=uuu->z[2];
                                                cu_z=uuu->z[1];
                                        }
                        
                        


        
        
                        
//                      double slopev = cv->GetChain(mpv[chvid])->avg_slope;
//                      double bzv = cv->GetChain(mpv[chvid])->start_z;
//                      double btv = cv->GetChain(mpv[chvid])->start_t;

                        //double slopev = cv->GetChain(mpv[chvid])->front_slope;
                        //double offsetv = cv->GetChain(mpv[chvid])->front_offset;
                        //double slopeu = cu->GetChain(mpu[chuid])->front_slope;
                        //double offsetu = cu->GetChain(mpu[chuid])->front_offset;

                                                
                                        if(vvv->entries>1 && vvv->e[0] < vvv->e[1]*0.2 )
                                        {
                                                //bz = vvv->z[1];
                                                //bt = vvv->t[1];
                                                tmpvz=vvv->z[1];
                                                cv_z=vvv->z[1];
                                        }
                                        else if(vvv->entries>2 && vvv->e[1] < vvv->e[2]*0.2 )
                                        {
                                        //      bz = vvv->z[2];
                                        //      bt = vvv->t[2];
                                                tmpvz=vvv->z[2];
                                                cv_z=vvv->z[2];
                                        }

                        //bt = cv->GetChain(mpv[chvid])->weighted_t;
                        
                
        //////

                
        vtx_z=cu_z<cv_z?cu_z:cv_z;
        
//      if(tmpuz>0 && tmpvz == 0){ vtx_z=tmpuz; cu_z = vtx_z; }
//      if(tmpvz>0 && tmpuz == 0){ vtx_z=tmpvz; cv_z = vtx_z;}
//      if(tmpuz>0 && tmpvz > 0){ vtx_z=tmpuz<tmpvz ? tmpuz:tmpvz;  cu_z=vtx_z; cv_z=vtx_z; }           
                
        //vtx_z-= 0.04;         
        
        
                
                        //if we have a hit at the right z spot, use that t, otherwise extrapolate
                        if(fabs(vtx_z-uuu->z[0])<0.01)
                        {
                                cu_t = uuu->t[0];
                        }else{
                                
                                double oldcut=cu_t;                     
                        //      cu_t = slopeu*(-bzu+cu_z) + btu;
        ////                    cu_t = slopeu * vtx_z + offsetu;
                
                        cu_t = uuu->interpolate(vtx_z);
                
                        
                                MSG("Finder",Msg::kDebug) << "using chain "<< mpu[chuid]<<" "<<cu->GetChain(mpu[chuid])->myId <<" for slope, moving ut from "<< oldcut <<" to "<< cu_t <<endl;
                        }
        

                        if(fabs(vtx_z-vvv->z[0])<0.01)
                        {
                                cv_t = vvv->t[0];
                        }else{
                                                
                                double oldcut=cv_t;
                                //cv_t = slopev*(-bzv+cv_z) + btv;              

                                //cv_t = slopev * vtx_z + offsetv;

                        cv_t = vvv->interpolate(vtx_z);


                                MSG("Finder",Msg::kDebug) << "using chain "<< mpv[chvid]<<" "<<cv->GetChain(mpv[chvid])->myId <<" for slope, moving ut from "<< oldcut <<" to "<< cv_t <<endl;
                        }
        
        

        
        //we should see if there is a chain path that passes close to this vertex
        //that is, if the found vertex is between the beginning and end of this chain, the chain should be split at the vertex, have the front part reversed
        //and the vertex should be adjusted to fall on this chain
        
        for(int q=0;q<2;q++)
        {
        
        ChainHelper * ch = q==0? cu:cv;
        
        for(unsigned int i=0;i<ch->finished.size();i++)
        {
                if(ch->finished[i].start_z > vtx_z) continue;
                std::vector<foundpath> paths = ch->FindPaths(ch->GetChain(ch->finished[i].myId));
                
                //find the chain where we would need to split....
                for(unsigned int p =0;p<paths.size();p++)
                {
        
                        //make sure the path is sufficiently far from the split point
                        if(! ( paths[p].start_z < vtx_z - 0.04 && paths[p].end_z > vtx_z + 0.04) )break;
        
                        int getp=-1;
                        for(unsigned int j=0;j<paths[p].path.size();j++)
                        {
                                getp = paths[p].path[j];
                                if(ch->GetChain(getp)->start_z < vtx_z -0.03 && ch->GetChain(getp)->end_z >= vtx_z+0.03)break;
                                getp=-1;
                        }
                        if(getp<0)continue;

                
        //              cout <<"WANT TO SPLIT "<<getp << " AT " <<vtx_z<<endl;
                
                        Chain * m = ch->GetChain(getp);
                        //if its close to the vertex, ajust vertex position
                        
                        int pt_bef=-1;
                        int pt_aft=-1;
                        for(int j=1;j<m->entries;j++)
                        {
                                if(m->z[j-1]<=vtx_z && m->z[j]>=vtx_z)
                                {
                                        pt_bef=j-1;
                                        pt_aft=j;
                                        break;
                                }
                        }
                
//                      cout<<"recalculating t from idx "<<pt_bef<<" "<<pt_aft<<endl;
                
                        if(pt_bef>-1 && pt_aft>-1)
                        {
                                double dt = m->t[pt_aft]-m->t[pt_bef];
                                double dz = m->z[pt_aft]-m->z[pt_bef];
                        
                //              cout << "dt "<<dt <<" dz "<<dz<<endl;
                        
                                double vt=0;
                                if(dt<0.0001)vt=m->t[pt_aft];
                                
                                else vt = (vtx_z-m->z[pt_bef]) * dt/dz + m->t[pt_bef];
                                
                                double oldt = q==0?cu_t:cv_t;
                                
                //              cout <<"old t "<<oldt << " new t "<<vt<<endl;
                                
                                if(fabs(oldt-vt)<0.05)
                                {
                                        if(q==0)cu_t=vt;
                                        if(q==1)cv_t=vt;
                                }
                        }
                
                        
                        
                        Chain * d = ch->SplitAt(m,vtx_z);
                        
                        if(d->children.size()<1)continue; //unable to split the chain...
                        
                        //d->parentChain=-1;
                        d->level=0;
                        for(unsigned int chi =0;chi<d->children.size();chi++)
                        {
                                Chain * cld = ch->GetChain(d->children[chi]);
                                if(cld){
                                        cld->parentChain=-1;
                                }
                        }
                        
                        Chain*cp = ch->GetChain(d->parentChain);
                        if(cp){
                                        for(unsigned int i=0;i<cp->children.size();i++)
                                        {
                                        if(cp->children[i]==d->myId)
                                        {
                                                cp->children.erase(cp->children.begin()+i);
                        
                        //                      cout<<"\t reversal causing removing of child reference from parent "<<d->parentChain<<endl;
                        
                                                break;  
                                        }
                                        }       
                        }
                        d->parentChain=-1;
                        
                        d->children.clear();
                        d->Reverse();
                        
                }
        
        }
        
        
        }
        

        
//      vtx_z-= 0.04;//0.025;  //event is usually before first hit.... make it better later, based on if this looks like a proton (indicating event forward of strip), or if this is a low hit strip(indicating hit well inside steel)
        
        vtx_u=cu_t;
        vtx_v=cv_t;

}




void Finder::FindIsolatedHits()
{

        printf("dont use Finder::FindIsolatedHits()\n");
        exit(1);

        double minthreshold = 0.5;  //peak must be at least this big
        double peakfraction = 0.85;  //fraction of energy in peak in this single strip

        //iterate until all isolated peaks are exhaused....


      

        for(std::multimap<double,int>::reverse_iterator iter = sorter_map.rbegin(); iter!=sorter_map.rend();iter++)
        {

        //for now, just look for max

        int maxi=-1;
        double maxe=0;

        int p=0;
        int c=0;

        double st=0;
        double sz=0;



        double locale=iter->first;

        maxi = iter->second;
        int i = maxi;
                        maxe=energy[i];
                        p=plane[i];
                        c=strip[i];
                        st=t[i];
                        sz=z[i];
        MSG("Finder",Msg::kDebug) <<"isolated hit peak at plane, cell " << p << " "<< c << " with energy " << maxe <<endl;



        if(maxe<minthreshold)break; //since iteration is in e order, all further entries will also be below threshold




        for(unsigned int i=0;i<plane.size();i++)
        {
//              if(plane[i]==p && abs(strip[i]-c)==1)locale+=energy[i];
//              if(strip[i]==c && abs(plane[i]-p)==1)locale+=energy[i];

                if(plane[i]==p && abs(strip[i]-c)==1)locale+=energy[i];


        }



        MSG("Finder",Msg::kDebug) << "peak "<< maxe/locale << " "<< maxe << " "<< p << " "<< sz << " " << sz/(1.0*p) <<endl;
        if(maxe/locale>peakfraction) //add particle
        {
                Particle ptemp;
                ptemp.begPlane=p;
                ptemp.endPlane=p;
                ptemp.begStrip=c;
                ptemp.endStrip=c;
                ptemp.stripEnergy=maxe - (locale-maxe) /4;  //assign the energy above the local average
                ptemp.type=2112; //2112 neut 2212 prot
                ptemp.endz=sz;
                ptemp.begz=sz;
                ptemp.endt=st;
                ptemp.begt=st;

                ptemp.plane.push_back(p);
                ptemp.strip.push_back(c);
                ptemp.energy.push_back(ptemp.stripEnergy);
                


                MSG("Finder",Msg::kDebug) <<"p at "<< p <<" "<< c <<" "<< st <<" "<< sz<<endl;

                particles.push_back(ptemp);

        
                //adjust the energy in the peak map, removing the value assigned to this hit, and replacing it with remaining value
//              sorter_map.erase((--iter).base());
//              sorter_map.insert(std::pair <double,int>((locale-maxe) /4, (int)energy.size()));        





        }else{continue;}

                
        }


}



std::vector<std::pair<int,int> > Finder::matchViews(std::vector<foundpath> pu, std::vector<foundpath> pv)
{
        std::vector<std::pair<int,int> >matches;
        
//      cout <<"views sizes "<<pu.size()<<" "<<pv.size()<<endl;
        for(unsigned int i=0;i<pu.size();i++)
        {

                std::vector<std::pair<int,int> >zmatches;

MSG("Finder",Msg::kDebug) <<"--- "<<i<<endl;

                //if(pu[i].entries<2)continue;
        
                std::vector<double> zoverlap;   
        
                for(unsigned int j=0;j<pv.size();j++)
                {

MSG("Finder",Msg::kDebug) <<"----- "<<j<<endl;

                //      if(pv[j].entries<2)continue;
                
                        if((pu[i].entries==1 && pv[j].entries>2 )|| (pu[i].entries>2 && pv[j].entries==1 ))
                        {
                        
                                zoverlap.push_back(0.0);
                                continue; //don't match single hit with path longer than 2
                        }
                        
                        
        MSG("Finder",Msg::kDebug)<< "looking at "<<i << "  "<<j<<endl;
        MSG("Finder",Msg::kDebug) << " pu ent " << pu[i].entries << " pv ent " << pv[j].entries<<endl;
                
                
                        double o=0.0;
                        if(pu[i].start_z>pv[j].end_z || pu[i].end_z < pv[j].start_z) //no overlap!
                        {
                                zoverlap.push_back(o);
                                
                        MSG("Finder",Msg::kDebug) << "no overlap " << pu[i].start_z << " to " << pu[i].end_z << " vs " << pv[j].start_z << " to "<< pv[j].end_z <<endl;
                                
                        //      continue;
                        }else{
                        
                        double start = pu[i].start_z < pv[j].start_z ? pv[j].start_z : pu[i].start_z;
                        double end = pu[i].end_z< pv[j].end_z ? pu[i].end_z : pv[j].end_z;
                        
                        o = fabs(end-start);
                        
                        double sm1= fabs(pu[i].end_z - pu[i].start_z);
                        double sm2= fabs(pv[j].end_z - pv[j].start_z);
                        double sm=0;
                        sm = sm1 > 0 ? sm1 : sm;
                        sm = sm2 > 0 && (sm2 < sm || sm ==0)? sm2 : sm;
                        
                        sm = sm1>sm2?sm1:sm2;
                        
                MSG("Finder",Msg::kDebug)<< "sm1 "<<sm1<<" sm2 "<<sm2<<" pl " << o << " sm " << sm <<endl;
                        
                        if(sm<0.00001)continue;
                        
                        if(o==0) //special case --- we already made sure that there was overlap ... so this happens when one of the views has 1 hit, so its pathlength is 0
                        o=sm=1;
                        
                        zoverlap.push_back(o / sm);  // recording fraction of overlap over smaller path
                        
        MSG("Finder",Msg::kDebug)<< "overlap : " << o/sm<< " smallest path length "<<sm <<endl;
                        }
                }


                //find the best
                double bestz=0.0;
                for(unsigned int j=0;j<zoverlap.size();j++)
                        bestz = bestz > zoverlap[j]  ? bestz : zoverlap[j] > 0.999 ? bestz : zoverlap[j]; //if its a complete overlap, we don't want to compare against it...
                
                        //require minimum overlap of 50%
                        if(bestz < 0.4 && bestz!=0.0) continue;
                
                
                //find all matches within 30% of best
                for(unsigned int j=0;j<zoverlap.size();j++)
                {
                        if((bestz - zoverlap[j]) < 0.3 * bestz || zoverlap[j] > 0.999)
                        {
                                zmatches.push_back(std::pair<int,int>(i,j));
                                MSG("Finder",Msg::kDebug)<< "\tkeeping "<<i<<" "<<j<<endl;
                        }
                }
                
                if(zmatches.size()==0)continue;
                
/*              if(zmatches.size()>1) //we need to pick a single match, so now look at energy and pick the closest!
                {
                        double beste = 100000.0; // distance in energy from ith pu
                        int besteidx=-1;
                        for(unsigned int j=0;j<zmatches.size();j++)
                        {
                                double de =  fabs(pu[zmatches[j].first].energy - pv[zmatches[j].second].energy);
                                if(beste > de)
                                {
                                        beste = de;
                                        besteidx=j;
                                }
                        
                        }
                        
                        if(besteidx>-1)
                        matches.push_back(zmatches[besteidx]);
                        
                }else{
                        matches.push_back(zmatches[0]);
                }               
*/



                if(zmatches.size()>1) //we need to pick a single match, so now look at energy / # entries and pick the closest!
                {
                        int closest=100000;
                        int cidx=-1;
                        for(unsigned int j=0;j<zmatches.size();j++)
                        {
                        //      int d = fabs(pu[zmatches[j].first].energy / pu[zmatches[j].first].entries - pv[zmatches[j].first].energy / pv[zmatches[j].second].entries);
                                
                        //      d=d/fabs(pu[zmatches[j].first].entries - pv[zmatches[j].second].entries);
                                
                                
                                int d = (int)(TMath::Abs(pu[zmatches[j].first].entries - pv[zmatches[j].second].entries));

                                if(d<closest)
                                {
                                        closest=d;
                                        cidx=j;
                                }
                        
                        }
                        
                        if(cidx>-1)
                        {
                                matches.push_back(zmatches[cidx]);
                                continue;                               
                        }


                }



        



                if (zmatches.size()>1) {
                        MSG("Finder",Msg::kWarning) << "Degeneracy not resolved in making 3d particle tracks!"<<endl;
                }
                for(unsigned int j=0;j<zmatches.size();j++)
                        matches.push_back(zmatches[j]);

        }


        return matches;
}


void Finder::Weave(ChainHelper * chu, ChainHelper * chv)
{
        std::vector<foundpath> pu = GetPaths(chu);
        std::vector<foundpath> pv = GetPaths(chv);

        if(pu.size()<1 || pv.size()<1)return;

        
        DumpPaths(pu);
        DumpPaths(pv);

        //////////////////////////////
        //we should clean each found path vector to remove paths which differ only by the last chain
        //as happens when a muon brems
        
        for(int i=0;i<2;i++){
        
        std::vector<foundpath>::iterator it1;
        std::vector<foundpath>::iterator it2;   
        std::vector<foundpath> paths=i==0?pu:pv;
        if(paths.size()<2)continue;
        
        double minzlength=1.0; //distance in m before we can delete...
        
        it1 =paths.begin();
        while(it1!=paths.end())
        //for(it1 =paths.begin();it1!=paths.end();it1++)
        {
                if(it1->path.size()<2)
                {
                        it1++;
                        continue;
                }
                int diddel=0;
                
                it2=it1;
                it2++;
                while(it2!=paths.end())
                //for(it2 =it1;it2!=paths.end();it2++)
                {
                        diddel=0;
                        if(it2->path.size()<2)
                        {
                                it2++;
                                continue;
                        }

                        
                        //MSG("Finder",Msg::kDebug)<<"paths start / end /it1 /it2 "<<&(paths.begin())<<" "<<&(paths.end())<<" "<<&it1<<" "<<&it2<<endl;
                        
        
                        
                        //MSG("Finder",Msg::kDebug)<<"Comparing with sizes "<<it1->path.size()<<" "<<it2->path.size()<<endl;
                        
                        int close=1;
                        
                //      int larger=it1->path.size()>it2->path.size()?0:1;
                        int larger=it1->end_z-it1->start_z > it2->end_z-it2->start_z ? 0:1;
                
                
                        std::vector<int> maxpi = larger==0?it1->path:it2->path;
                        std::vector<int> maxpj = larger==1?it1->path:it2->path; 
                        
                        for(unsigned int k=0;k<maxpj.size()-1;k++)
                        {
                        MSG("Finder",Msg::kDebug)<<"chains "<<maxpi[k]<<" "<<maxpj[k]<<endl;
                                if(maxpi[k]!=maxpj[k])
                                {
                                        close=0;
                                        break;
                                }
                        }
                        
                        //delete the smaller one?
                        std::vector<foundpath>::iterator itsmall=larger==0?it2:it1;
                        double dist = itsmall->end_z-itsmall->start_z;
                        if(close && dist > minzlength)
                        {
                        /*MSG("Finder",Msg::kError)*/
                        /*cout<<"removing close chain (";
                        for(unsigned int q=0;q<it1->path.size();q++)cout <<it1->path[q]<<" ";
                        cout<<") (";
                        for(unsigned int q=0;q<it2->path.size();q++)cout <<it2->path[q]<<" ";
                        cout<<") \n";                   
                        */

                                int both=it1==it2?1:0;
                                if(larger==1)
                                {
                                        it1=paths.erase(it1);
                                        if(both)it2=it1;
                                }
                                if(larger==0)
                                {
                                        it2=paths.erase(it2);
                                        if(both)it2=it1;
                                }
                                
                                //something doesnt work right....!
                                //do this sucky hack
                                diddel=1;
                                
                                                //MSG("Finder",Msg::kError)<<"current pointers start / end /it1 /it2 "<<&(paths.begin())<<" "<<&(paths.end())<<" "<<&it1<<" "<<&it2<<endl;
                        }
                        
                        if(diddel)
                        {
                                it1=paths.begin();
                                it2=it1;
                        }
                        it2++;
                        
                }       
                it1++;
                
        }
                if(i==0)pu=paths;
                if(i==1)pv=paths;
        }//loop over views
        ///////////////////////////////

        if(pu.size()<1 || pv.size()<1)
        {
                MSG("Finder",Msg::kError)<<"Removal of similar paths results in empty view!\n"; 
                return;
        }


        //cout <<"*********Currently have "<<mypoh->particles3d.size()<<" particles in POH\n";

        DumpPaths(pu);
        DumpPaths(pv);
        
        
        //now try to find best matches!
        std::vector<std::pair<int,int> >matches = matchViews(pu,pv);    
        for(unsigned int i=0;i<matches.size();i++)
        {
                MSG("Finder",Msg::kDebug) << "found pair " << matches[i].first << " " <<matches[i].second<<endl;
        }

        //and check the other view first as well!
        std::vector<std::pair<int,int> >matcheso = matchViews(pv,pu);
        for(unsigned int i=0;i<matcheso.size();i++)
        {
                MSG("Finder",Msg::kDebug) << "found pair " << matcheso[i].second << " " <<matcheso[i].first<<endl;
        }
        for(unsigned int i=0;i<matcheso.size();i++)
        {
                int found =0;
                std::pair<int,int> tmp = make_pair(matcheso[i].second, matcheso[i].first); //need to flip the pair
                for(unsigned int j=0;j<matches.size();j++)
                {
                        
                        if(matches[j]==tmp)
                        {
                                found=1;
                                break;
                        }
                }
                if(!found)matches.push_back(tmp);
        
        }
        

        
//      printf("!!1 matches size %d\n",matches.size());
        
        
        if(matches.size()<1)return;
        
        ostringstream os;
        os << "Candidate matches on z overlap : ";
        for(unsigned int i=0;i<matches.size();i++)
        {
                os << matches[i].first <<"-"<<matches[i].second<<" ";
        }
        MSG("Finder",Msg::kDebug) << os.str()<<endl;    



        //when getting rid of a match... make sure its the best match!  (as we can share a chain in one view!)
        //need to see if all matches are "fair"
        std::vector<std::pair<int,int> >matchestmp;
        std::vector<int>multu;
        std::vector<int>multv;
        for(unsigned int i=0;i<matches.size();i++)
        {
                int mult0=0;
                int mult1=0;
                for(unsigned int j=0;j<matches.size();j++)
                {
                        if(i==j)continue;
                        if(matches[i].first==matches[j].first)mult0++;
                        if(matches[i].second==matches[j].second)mult1++;
                }
                if(! ( mult0>0 && mult1>0 )) // no degeneracy
                {
                        matchestmp.push_back(matches[i]);
                        multu.push_back(mult0);
                        multv.push_back(mult1);
                }else if( mult0>0 && mult1>0 ) //degeneracy, but see if we really want it
                {
                        //see if the reason why both of these are degenerate is because they are really long
                        //for example, this can happen to a long muon in a DIS CC
                        //we want each of these to be the longest out of all matches using each of these paths
                        int pass=1;
                        
                        for(unsigned int k=0;k<matches.size();k++)
                        {
                                if(matches[i].second!=matches[k].second)continue;
                                if(pu[matches[i].first].end_z-pu[matches[i].first].start_z < pu[matches[k].first].end_z-pu[matches[k].first].start_z)
                                {
                                        pass=0;
                                        break;
                                }
                        }
                        for(unsigned int k=0;k<matches.size();k++)
                        {
                                if(matches[i].first!=matches[k].first)continue;
                                if(pv[matches[i].second].end_z-pv[matches[i].second].start_z < pv[matches[k].second].end_z-pv[matches[k].second].start_z)
                                {
                                        pass=0;
                                        break;
                                }
                        }
                
                        if(pass)
                        {
                                matchestmp.push_back(matches[i]);
                                multu.push_back(mult0);
                                multv.push_back(mult1);
                        }
                }
                
                MSG("Finder",Msg::kDebug) << "pair " << matches[i].first <<" " << matches[i].second << " matches " <<mult0 <<" "<< mult1<<"\n";
        }
        matches=matchestmp;

        //cout <<"*********Currently have "<<mypoh->particles3d.size()<<" particles in POH\n";
        
        //we now have matches...
        //for each, construct a 3d path view
        
        
        std::vector<Particle3D> p3v;
        for(unsigned int i=0;i<matches.size();i++)
        {
                MSG("Finder",Msg::kDebug) << "Making 3d particle with paths "<<matches[i].first<<" "<<matches[i].second<<"\n";
        
                Particle3D p3 = Make3DParticle(pu[matches[i].first].path, pv[matches[i].second].path, chu, chv,multu[i],multv[i]) ;
                
                //we also need to save the clusters....
                
                for(unsigned int j=0;j<p3.view.size();j++)
                {
                        Managed::ClusterManager * cluster_manager = &mypoh->clusterer;//p3.view[j]==2?&clustermanager_u:&clustermanager_v;
                        if(p3.chainhit[j]>-1)
                        {
                                ChainHelper *ch = p3.view[j]==2?chu:chv;
                                int oldid = ch->GetChain(p3.chain[j])->cluster_id[p3.chainhit[j]];
                                int newid = cluster_manager->SaveCluster(oldid,p3.e[j],3);
                                
                                Managed::ManagedCluster * cluster = cluster_manager->GetCluster(newid);
                                if(cluster)
                                {
                                        //printf("SSSS saved %d  with e %f new id %d with e %f \n",oldid,p3.e[j],newid,cluster->e);
                                
                                        p3.e[j]=cluster->e;
                                
                                        p3.chainhit[j]=newid;
                                        
                                        ch->ChangeHitId(oldid,newid);
                        
                                }else{
                                        //printf("NNNN cluster %d doesn't exist! removing e %f from particle \n",oldid,p3.e[j]);
                                        p3.e[j]=0;
                                        p3.chainhit[j]=-1;
                                
                                        
                                }
                        }
                        
                                
                
                }
                
                if(p3.entries && p3.sum_e>1e-5)
                {
                        //printf("adding p3d with e %f\n",p3.sum_e);
                        p3v.push_back(p3);
                }
        }

        
//printf("!!2 p3 size %d \n",p3v.size());

        
        
        
        
        
        

        std::vector<Particle3D* > p3d;  
        for(unsigned int i=0;i<p3v.size();i++)
                p3d.push_back(&p3v[i]); 
        
//      p3d=SetShared(p3d);     

//      shareholder.Reset();
//      shareholder.BulkInsert(p3d);

        /*
        DumpParticle3D(p3d);


        std::vector<Particle3D*> pout;
        std::vector<Particle3D*> toshare;
        for(unsigned int i=0;i<p3d.size();i++)
        {
                //cout<<"on particle3d "<<i<<endl;
        
                std::vector<Particle3D* > po = ProcessParticle3D1(p3d[i]);
                for(unsigned int j=0;j<po.size();j++)
                {
                        //if(po[j]->hasShared()<1)
                                pout.push_back(po[j]);  
                        //else toshare.push_back(po[j]);
                }
        
        }
        */
        
/*      
        cout <<"PROCESSING COMPLETE\n";
        cout <<"not shared\n";
                DumpParticle3D(pout);
        cout <<"has shared\n";

                DumpParticle3D(toshare);
                
        cout << "UNSHARING ENERGY\n";
*/      
        

/*

        int its=5; //max iterations to try to remove shared hits!
        while(toshare.size()>0 && its>0)
        {
                std::vector<Particle3D*> tmp;
                
                //some code to process them, returning a vector....
                
                for(unsigned int i=0;i<toshare.size();i++)
                        tmp.push_back(toshare[i]);
                toshare.clear();
                
                RemoveSharedHits(tmp);
                
                for(unsigned int i=0;i<tmp.size();i++)
                {
                        if(tmp[i]->hasShared())
                                toshare.push_back(tmp[i]);
                        else pout.push_back(tmp[i]);
                
                }
                
                its--;
        }
        
        if(toshare.size()>0)
        MSG("Finder",Msg::kWarning)<<"particles lost .... could not uniquely assign shared energy !"<<endl;
        
        */


/*
        std::vector<Particle3D*> pout;
        for(unsigned int i=0;i<p3d.size();i++)
        {
                std::vector<Particle3D* > po = ProcessParticle3D(p3d[i]);       //returned vector may not have shared hits....
                
                std::vector<Particle3D*> forsharing;
                
                for(unsigned int j=0;j<po.size();j++)
                        pout.push_back(po[j]);
                
                for(unsigned int j=0;j<pout.size();j++)
                        forsharing.push_back(pout[j]);
                
                for(unsigned int j=i+1;j<p3v.size();j++)
                        forsharing.push_back(p3d[j]);

                        
                forsharing=SetShared(forsharing); //to reset if a hit is shared....
                
                
        }
        
        
        
*/      
/*
        cout <<"CLEANING\n";
        for(unsigned int i=0;i<pout.size();i++)pout[i]->Clean(); //remove 0 energy hits...
        DumpParticle3D(pout);
*/      



        //we need to recalculate values after sharing hits
        //but don't try to find more unique muons
        //we just want to make sure that the type hasn't changed!
        std::vector<Particle3D*> pout1;
        for(unsigned int i=0;i<p3d.size();i++)
        {
                std::vector<Particle3D* > po = ProcessParticle3D1(p3d[i],0);
                for(unsigned int j=0;j<po.size();j++)
                {
                        pout1.push_back(po[j]); 
                }
        
        }
        
        //printf("weave found %d\n",pout1.size());

//      pout1=p3d;
        //look for straglers (large hit, not associated with any chain.... condider these to be neutrons...)
        FindNeutrons(pout1);    
        //printf("weave found %d\n",pout1.size());
        

        FinalizeParticles3D(pout1); //final computation of values, etc...
        //printf("weave found %d\n",pout1.size());


        for(unsigned int i=0;i<pout1.size();i++)
        {
                if(pout1[i]->sum_e>0.01)
                        mypoh->AddParticle3D(*(pout1[i]));
        }
        
        //printf("weave found %d\n",pout1.size());

                
        DumpParticle3D(pout1);  

                
        //RecordLostHits(pout1);        
                
//              printf("!!3 pout size %d\n",pout.size());
                
                
//              printf("poh entries %d\n", mypoh->particles3d1->GetEntries());
                                

//      for(unsigned int i=0;i<p3d.size();i++)
//              mypoh->AddParticle3D(*(p3d[i]));        



        //do a check to see if the total particle energy is less than the total event energy...
        double parte=0;
        for(unsigned int i=0;i<p3d.size();i++) parte+=p3d[i]->sum_e;
        
        if(parte-mypoh->event.visenergy>.001)printf("!!!!!!!!!!!!!!!!!!\n********************\nenergy mismatch!!!! particle sum %f   energy total %f\n!!!!!!!!!!!!!!!!!!\n********************\n",parte,mypoh->event.visenergy);
        
        

        
}




void Finder::FindNeutrons(std::vector<Particle3D*> & pout)
{
        std::vector<Particle3D*> pout1=pout;
        //mark all clusters as unused
        Managed::ClusterManager  *clhu = &mypoh->clusterer;
        Managed::ClusterManager  *clhv = &mypoh->clusterer;

/*
        ChainHelper *cht = &mypoh->chu;
        for(int i=0;i<cht->finished.size();i++)
        {
                printf("chain %d \n",cht->finished[i].myId);
                for(int j=0;j<cht->finished[i].cluster_id.size();j++)
                        printf("%d ",cht->finished[i].cluster_id[j]);
                printf("\n");
        
        }


        cht =& mypoh->chv;
        for(int i=0;i<cht->finished.size();i++)
        {
                printf("chain %d \n",cht->finished[i].myId);
                for(int j=0;j<cht->finished[i].cluster_id.size();j++)
                        printf("%d ",cht->finished[i].cluster_id[j]);
                printf("\n");
        
        }
*/              


//cant do this check here.. because other particles might already be in poh and are not in pout 
        for(unsigned int i=0;i<clhu->clusters.size();i++)clhu->clusters[i].inuse=0;
        for(unsigned int i=0;i<clhv->clusters.size();i++)clhv->clusters[i].inuse=0;
        
        //go through existing partcle3d's marking used clusters as used

        for(unsigned int i=0;i<mypoh->particles3d.size();i++)pout1.push_back(&(mypoh->particles3d[i]));

        for(unsigned int i=0;i<pout1.size();i++)
        {
                for(int j=0;j<pout1[i]->entries;j++)
                {
                        ChainHelper *ch = pout1[i]->view[j]==2? &mypoh->chu : pout1[i]->view[j]==3?& mypoh->chv:0;
                        if(!ch)continue;
                        
                        //printf("c %d ch %d\n",pout1[i]->chain[j],pout1[i]->chainhit[j]);
                        
                        Chain * c = ch->GetChain(pout1[i]->chain[j]);
                        if(!c)
                        {
                                //cout<<"chain information lost\n";
                                continue;
                        }
                        
                        if(pout1[i]->chainhit[j]<0)continue;  //not available
                        
                        int id = c->cluster_id[pout1[i]->chainhit[j]];
                        
                        //cout <<"chain/hit "<<pout1[i]->chain[j] << " " <<pout1[i]->chainhit[j]<<endl;
                        
                        Managed::ClusterManager  * clh = pout1[i]->view[j]==2? clhu : pout1[i]->view[j]==3? clhv:0;
                        if(!clh)continue;
                        
                        Managed::ManagedCluster *cluster=clh->GetCluster(id);
                        if(cluster)cluster->inuse=1;
                        //cout<<"in use\n";
                
                }
        
        }

/*

        for(unsigned int i=0;i<mypoh->particles3d.size();i++)pout1.push_back(&(mypoh->particles3d[i]));

        for(unsigned int i=0;i<pout1.size();i++)
        {
                for(unsigned int j=0;j<pout1[i]->entries;j++)
                {
                        ChainHelper *ch = pout1[i]->view[j]==2? &mypoh->chu : pout1[i]->view[j]==3?& mypoh->chv:0;
                        if(!ch)continue;
                        
                        printf("c %d ch %d\n",pout1[i]->chain[j],pout1[i]->chainhit[j]);
                        
                        Chain * c = ch->GetChain(pout1[i]->chain[j]);
                        if(!c)
                        {
                                //cout<<"chain information lost\n";
                                continue;
                        }
                        
                        if(pout1[i]->chainhit[j]<0)continue;  //not available
                        
                        int id = c->cluster_id[pout1[i]->chainhit[j]];
                        
                        cout <<"chain/hit "<<pout1[i]->chain[j] << " " <<pout1[i]->chainhit[j]<<endl;
                        
                        Managed::ClusterManager  * clh = pout1[i]->view[j]==2? clhu : pout1[i]->view[j]==3? clhv:0;
                        if(!clh)continue;
                        
                        Managed::ManagedCluster *cluster=clh->GetCluster(id);
                        if(cluster)cluster->inuse=1;
                        cout<<"in use\n";
                
                }
        
        }


*/

        //find unused clusters with energy > some thresh
        //make these into particles...
/*      
        cout << "//////////////////////////////////////////////\nunused clutsers"<<endl;
        cout <<"U\n";
        
        for(unsigned int i=0;i<clhu->clusters.size();i++)
        {
                Managed::ManagedCluster *c = & clhu->clusters[i];
                if(c->inuse)continue;
                cout << "Cluster "<<c->id << " z " << c->z << " t " << c->t << " E " << c->e<<endl;
        }
        cout <<"V\n";
        for(unsigned int i=0;i<clhv->clusters.size();i++)
        {
                Managed::ManagedCluster *c = & clhv->clusters[i];
                if(c->inuse)continue;
                cout << "Cluster "<<c->id << " z " << c->z << " t " << c->t << " E " << c->e<<endl;
        }
        
        cout << "//////////////////////////////////////////////\n";

*/
/*
        //for now just make 2d neutrons... (don't match views..)
        
        double neut_thresh=3.0;
        
        for(int q=0;q<2;q++)
        {
                Managed::ClusterManager * ch = q==0?clhu:clhv;
                for(unsigned int i=0;i<ch->clusters.size();i++)
                {
                        Managed::ManagedCluster *c = & ch->clusters[i];
                        if(c->inuse)continue;
                        if(c->e<neut_thresh)continue;
                        
                                Particle3D * pnew = new Particle3D();
                                pnew->particletype=Particle3D::neutron;
                                double u = q==0? c->t : 0.0;
                                double v = q==1? c->t : 0.0;
                                pnew->add_to_back(u, v, c->z, c->e, -1,-1, q==0?2:3, 0);
                                pout.push_back(pnew);                   
                }
        
        }
*/
        //for now just make 2d neutrons... (don't match views..)
        
        double neut_thresh=3.0;
        
        for(int q=0;q<2;q++)
        {
                //Managed::ClusterManager * ch = q==0?clhu:clhv;
                
                std::vector<int> chidx = clhu->GetViewIndex(q==0?2:3);
                std::vector<int> chidx_other = clhu->GetViewIndex(q==1?2:3);
                
                for(unsigned int i=0;i<chidx.size();i++)
                {
                        Managed::ManagedCluster *c = & clhu->clusters[chidx[i]];
                        if(c->inuse)continue;
                        if(c->e<neut_thresh)continue;
                        
                        
                        //Managed::ClusterManager * ch_other = q==1?clhu:clhv;
                        std::vector<int>other_view;
                        for(unsigned int j=0;j<chidx_other.size();j++)
                        {
                                Managed::ManagedCluster *cother = & clhu->clusters[chidx_other[j]];
                                if(TMath::Abs(cother->z - c->z)<0.06 && !cother->inuse)
                                        other_view.push_back(j);
                        }       
                                
                        if(other_view.size()==0)
                        {       
                                Particle3D * pnew = new Particle3D();
                                pnew->particletype=Particle3D::neutron;
                                double u = q==0? c->t : vtx_u;
                                double v = q==1? c->t : vtx_v;
                                pnew->add_to_back(u, v, c->z, c->e, -1,-1, q==0?2:3, 0);
                                pout.push_back(pnew);   
                                c->inuse=1;             
                        }else{
                                //don't interpolate....
                                
                                Managed::ManagedCluster *cother0 = & clhu->clusters[chidx_other[other_view[0]]];
                                
                                double o_t = cother0->t;
                                if(other_view.size()==2) 
                                {
                                        Managed::ManagedCluster *cother1 = & clhu->clusters[chidx_other[other_view[1]]];

                                        o_t += cother1->t;
                                        o_t /=2.0;
                                }
                                
                                
                                double u = q==0? c->t : o_t;
                                double v = q==1? c->t : o_t;                            
                        
                        
                                Particle3D * pnew = new Particle3D();
                                pnew->particletype=Particle3D::neutron;
                                pnew->add_to_back(u, v, c->z, c->e, -1,-1, q==0?2:3, 0);
                                c->inuse=1;
                                pnew->add_to_back(u, v, cother0->z, cother0->e, -1,-1, q==1?2:3, 0);
                                cother0->inuse=1;                               
                                if(other_view.size()==2)
                                {
                                        Managed::ManagedCluster *cother1 = & clhu->clusters[chidx_other[other_view[1]]];
                                        pnew->add_to_back(u, v, cother1->z, cother1->e, -1,-1, q==1?2:3, 0);
                                        cother1->inuse=1;
                                }
                                
                                pout.push_back(pnew);   
                        }
        
        
                }
        
        }


}



/**determine what hits are unused
*/
void Finder::RecordLostHits(std::vector<Particle3D*> p3d)
{
        for(int i=0;i<2;i++)
        {
                Managed::ClusterManager * ch = i==0? &mypoh->clusterer : &mypoh->clusterer;
                
                for(unsigned int j=0;j<ch->clusters.size();j++)ch->clusters[j].inuse=0;
        }
        
        for(unsigned int i=0;i<p3d.size();i++)
        {
                Particle3D * p = p3d[i];
                for(unsigned int j=0;j<p->chain.size();j++)
                {
                        Managed::ClusterManager * clh = p->view[j]==2 ? &mypoh->clusterer :  &mypoh->clusterer ;
                        ChainHelper * ch = p->view[j]==2 ? &mypoh->chu :  &mypoh->chv ;
                        Chain * chain = ch->GetChain(p->chain[j]);
                        if(p->chainhit[j]>-1)
                                clh->GetCluster(chain->cluster_id[p->chainhit[j]])->inuse=1;
                }
        }
        
        for(int i=0;i<2;i++)
        {
                Managed::ClusterManager * ch = i==0? &mypoh->clusterer : &mypoh->clusterer;
                
                for(unsigned int j=0;j<ch->clusters.size();j++)
                {
                        if(ch->clusters[j].inuse)continue;
                        for(unsigned int k=0;k<ch->clusters[j].hite.size();k++)
                        {
                                mypoh->event.unused_e +=ch->clusters[j].hite[k];
                                mypoh->event.unused_strips++;
                        
                        }
                }
        }

        mypoh->event.unused_e_avg =0;
        if(mypoh->event.unused_strips>0)        
                mypoh->event.unused_e_avg = mypoh->event.unused_e / mypoh->event.unused_strips;
        
        for(int i=0;i<2;i++)
        {
                Managed::ClusterManager * ch = i==0? &mypoh->clusterer : &mypoh->clusterer;
                
                for(unsigned int j=0;j<ch->clusters.size();j++)
                {
                        if(ch->clusters[j].inuse)continue;
                        for(unsigned int k=0;k<ch->clusters[j].hite.size();k++)
                        {
                                mypoh->event.unused_e_rms += (ch->clusters[j].hite[k] - mypoh->event.unused_e_avg)*(ch->clusters[j].hite[k] - mypoh->event.unused_e_avg);
                        
                        }
                }
        }
        
        mypoh->event.unused_e_rms=sqrt(mypoh->event.unused_e_rms);
        if(mypoh->event.unused_e_avg>0)mypoh->event.unused_e_rms/=mypoh->event.unused_e_avg;
        
}



/** Finalize each Particle3D 
*
*       Finalize each Particle3D by determining final type and calculating final values based on that type
*/

void Finder::FinalizeParticles3D(std::vector<Particle3D*>pout)
{
        for(unsigned int i=0;i<pout.size();i++)
        {
                Particle3D * p = pout[i];
                if(!p)continue;
        
                p->Clean();
        
                
                if(p->entries>0)
                {               
                        double sum_z=0;
                        double sum_z_z=0;
                        double sum_u=0;
                        double sum_z_u=0;
                        double sum_v=0;
                        double sum_z_v=0;
                                

                        double n = 5 <p->entries? 5:p->entries;
                        if(n<1)continue;
                        for(int j=0;j<n;j++)
                        {
                
                                sum_z+=p->z[j];
                                sum_z_z+=p->z[j]*p->z[j];
                                sum_u+=p->u[j];
                                sum_z_u+=p->z[j]*p->u[j];
                                sum_v+=p->v[j];
                                sum_z_v+=p->z[j]*p->v[j];
        
                        }
                
                        if(n==1)//add the vertex as a point...
                        {
                        
                                sum_z+=vtx_z;
                                sum_z_z+=vtx_z*vtx_z;
                                sum_u+=vtx_u;
                                sum_z_u+=vtx_z*vtx_u;
                                sum_v+=vtx_v;
                                sum_z_v+=vtx_z*vtx_v;
                                n++;
                        }
                
                        double denom = (n*sum_z_z-(sum_z)*(sum_z)) < 1e-10 ? 0 :  (n*sum_z_z-(sum_z)*(sum_z));
                        double au=0;
                        double bu=0;
                        double av=0;
                        double bv=0;
                                
                        if(denom==0)            
                        {
                          au=std::numeric_limits<double>::infinity();
                          bu=std::numeric_limits<double>::infinity();
                          av=std::numeric_limits<double>::infinity();
                          bv=std::numeric_limits<double>::infinity();



                        }else{
                                au=(sum_u*sum_z_z-sum_z*sum_z_u)/denom;
                                bu=(n*sum_z_u-sum_z*sum_u)/denom;
                                av=(sum_v*sum_z_z-sum_z*sum_z_v)/denom;
                                bv=(n*sum_z_v-sum_z*sum_v)/denom;
                        }





                //      avg_slope = b;
                //      avg_offset = a;
                
                        double theta = TMath::ATan2(bv,bu);
                        double phi = TMath::ATan(bu*bu+bv*bv);
                
                
                        p->theta=theta;
                        p->phi=phi;
                
                }
                
                
                //for now... if everything else is muon energy like...
                p->calibrated_energy = p->sum_e * (1.46676) / 0.0241;
                
                
                //attempt em fit!
                if(p->particletype!=Particle3D::muon && p->entries>2 && p->emfit_a==0 /*already done?*/)
                {
                
                        ShwFit f;
                        
                        double startz=0;
                        
                        startz=p->z[0];
                        
                        for(int i=0;i<p->entries;i++)
                        {
                                //double planewidth=0.035;
                                
                //              f.Insert(p->e[i],(int)((p->z[i]-startz)/planewidth));
                        
                //              cout <<"inserting "<<p->e[i]<<" at plane "<<(int)((p->z[i]-startz)/planewidth)<<endl;
                
                //for now, do the simple thing and assume each hit is in the next plane...
                        
                                //                      f.Insert(p->e[i],i+1);
                        
                        
                                f.Insert3d(p->e[i],p->u[i],p->v[i],p->z[i]);
                        
                                //cout <<"inserting "<<p->e[i]<<" at plane "<<i+1<<endl;
                
                        
                        }
                        
                        //f.Fit();
                        
                        f.Fit3d();
                        
                        MSG("Finder",Msg::kDebug) <<"EM FIT! " << f.par_a << " " << f.par_b << " +- "<<f.par_b_err<<" "<<f.par_e0<<endl;
                
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
        
//      if(p)           printf("finder-- %f %d %f\n",p->cmp_chisq,p->cmp_ndf, p->peakdiff);             
                
                //      if(f.par_b<0.31 || f.par_b > 0.69) //not too emlike!
                //              if(p->particletype==Particle3D::electron)p->particletype=Particle3D::other;
                //      if(f.par_b>0.31 && f.par_b < 0.69) //force it to be emlike if not already
                //              p->particletype=Particle3D::electron;
                                
                //      if(p->particletype==Particle3D::electron && f.par_e0/25/60. < 0.8) p->particletype=Particle3D::other; //<0.8gev elec is hard to find... so thats probably not it!
                /*
                        if(TMath::Prob(p->emfit_chisq,p->emfit_chisq/p->emfit_chisq_ndf)>0.8)
                        {
                                p->particletype=Particle3D::electron;
                        }
                        else
                        {
                                p->particletype=Particle3D::other;
                        }
                */
                
                
                //      if(f.prob>0.9)  //good fit
                        if(f.conv && 
                        f.par_b - f.par_b_err*2.0 < 0.6 && f.par_b + f.par_b_err*2.0 > 0.4 //2.0 sigma from correct b value     range
                        && f.par_b_err < f.par_b //< 0.4
                        && f.par_e0/60./meupergev>0.25) //at least 0.25 gev cal e  to trust it!
                        {
                                p->particletype=Particle3D::electron;           
                        }
                        else
                        {
                                p->particletype=Particle3D::other;
                                
                //              f.Fit3dx2();
                        }
                
                }
        
                
                //proton filter
                for(unsigned int i=0;i<p->types.size();i++)
                {
                        if(p->types[i].type==ParticleType::prot)
                        {
                                //if 80% of the particle's energy is in the last two hits....
                                
                                if(p->sum_e<0.001)continue;
                                
                                if(p->entries<3)continue;
                                
                                double efrac =( p->e[p->entries-1] + p->e[p->entries-2] ) / p->sum_e;
                                
                                if(efrac > 0.6)
                                {
                                        p->particletype=Particle3D::proton;
                                        
                                        p->calibrated_energy=p->sum_e*60.;
                                        break;
                                }
                        }
                
                }
        
        
        }
}




std::vector<Particle3D> Finder::shareEnergy(std::vector<Particle3D> p3v)
{

        //now iterate over found particles
        //and try to divide up shared energy!
        
        
        std::multimap< std::pair<int, std::pair<int,int> >, Particle3D* > pmap;
        for(int i=0;i<(int)p3v.size();i++)
        for(int j=0;j<p3v[i].entries;j++)
        {
                pmap.insert(make_pair( make_pair( p3v[i].view[j],  make_pair(p3v[i].chain[j],p3v[i].chainhit[j])),&(p3v[i])) );
        
                
        }
        
        
        std::pair<int,int>last;
        int lastcount=0;
        int lastview=0; 
        std::vector<Particle3D*> shared3ds;
        std::multimap< std::pair<int, std::pair<int,int> >, Particle3D* >::reverse_iterator it1;
        for(it1=pmap.rbegin();it1!=pmap.rend();it1++)
        {
                if(lastcount==0)
                {
                        last=it1->first.second;
                        lastcount=1;
                        lastview=it1->first.first;
                        shared3ds.clear();
                        shared3ds.push_back(it1->second);
                        continue;       
                }
                
                if(lastcount==1 && ( last!=it1->first.second || lastview!=it1->first.first))
                {
                        last=it1->first.second;
                        lastview=it1->first.first;
                        shared3ds.clear();
                        shared3ds.push_back(it1->second);
                        continue;
                }
                
                if(last==it1->first.second && lastview==it1->first.first)
                {
                        lastcount++;
                        shared3ds.push_back(it1->second);
                        continue;
                }else{
                        //do what we need to do with all of the last ones...
                        
                        //find all with this key...
                        
                        //then run a function over the list of the effected particle3d's noting the chain and chain id
                        //this function will look for adjacent hits in each of the particle3ds and attempt to extrapolate
                        //energy will be divided up evenly
                        //start at the back (highest z) and work down... as most shared hits should be closer to the vertex
                        
                        for(unsigned int i=0;i<shared3ds.size();i++)
                        {
                                shared3ds[i]->SetShared(last.first, last.second);
                        }
                        
                        ShareHit(lastview,last.first,last.second,shared3ds);
                        
                        
                        
                        
                        //printf("cmp %d %d shared %d times\n",last.first, last.second,lastcount);
                        
                        
                        
                        //then start with the new one..
                        lastcount=1;
                        last=it1->first.second;
                        lastview=it1->first.first;
                        shared3ds.clear();
                        shared3ds.push_back(it1->second);
                }
                
                
        }
        
        if(lastcount>1)
        {
                for(unsigned int i=0;i<shared3ds.size();i++)
                {
                        shared3ds[i]->SetShared(last.first, last.second);
                }
                ShareHit(lastview,last.first,last.second,shared3ds);
                //printf("cmp %d %d shared %d times\n",last.first, last.second,lastcount);
        }

        return p3v;

}


std::vector<Particle3D*> Finder::SetShared(std::vector<Particle3D*> p3v)
{

        //now iterate over found particles
        //and try to divide up shared energy!
        
        
        std::multimap< std::pair<int, std::pair<int,int> >, Particle3D* > pmap;
        for(int i=0;i<(int)p3v.size();i++)
        for(int j=0;j<p3v[i]->entries;j++)
        {
                pmap.insert(make_pair( make_pair( p3v[i]->view[j],  make_pair(p3v[i]->chain[j],p3v[i]->chainhit[j])),p3v[i]) );
                p3v[i]->UnsetShared(p3v[i]->chain[j],p3v[i]->chainhit[j]);
                
        }
        
        
        std::pair<int,int>last;
        int lastcount=0;
        int lastview=0; 
        std::vector<Particle3D*> shared3ds;
        std::multimap< std::pair<int, std::pair<int,int> >, Particle3D* >::reverse_iterator it1;
        for(it1=pmap.rbegin();it1!=pmap.rend();it1++)
        {
        
                if(it1->second->ShareLocked(it1->first.second.first, it1->first.second.second)>0)continue;
        
                if(lastcount==0)
                {
                        last=it1->first.second;
                        lastcount=1;
                        lastview=it1->first.first;
                        shared3ds.clear();
                        shared3ds.push_back(it1->second);
                        continue;       
                }
                
                if(lastcount==1 && ( last!=it1->first.second || lastview!=it1->first.first))
                {
                        last=it1->first.second;
                        lastview=it1->first.first;
                        shared3ds.clear();
                        shared3ds.push_back(it1->second);
                        continue;
                }
                
                if(last==it1->first.second && lastview==it1->first.first )
                {
                        lastcount++;
                        shared3ds.push_back(it1->second);
                        continue;
                }else{
                        //do what we need to do with all of the last ones...
                        
                        //find all with this key...
                        
                        //then run a function over the list of the effected particle3d's noting the chain and chain id
                        //this function will look for adjacent hits in each of the particle3ds and attempt to extrapolate
                        //energy will be divided up evenly
                        //start at the back (highest z) and work down... as most shared hits should be closer to the vertex
                        
                        for(unsigned int i=0;i<shared3ds.size();i++)
                        {
                                shared3ds[i]->SetShared(last.first, last.second);
                        }
                        
                        //printf("cmp %d %d shared %d times\n",last.first, last.second,lastcount);
                        
                        
                        
                        //then start with the new one..
                        lastcount=1;
                        last=it1->first.second;
                        lastview=it1->first.first;
                        shared3ds.clear();
                        shared3ds.push_back(it1->second);
                }
                
                
        }
        
        if(lastcount>1)
        {
                for(unsigned int i=0;i<shared3ds.size();i++)
                {
                        shared3ds[i]->SetShared(last.first, last.second);
                }
                //printf("cmp %d %d shared %d times\n",last.first, last.second,lastcount);
        }

        return p3v;

}


void Finder::ShareHit(int view, int chain, int chainhit, std::vector<Particle3D*> shared)
{

        //for now, divide the energy based on the average of the next/previous energies
        //in a chain as a fraction of the total n/p energy in all of the shared chains....
        
        
        std::vector<double> es;
        for(unsigned int i=0;i<shared.size();i++)
        {
                int inview = view == 2 ? 3 :2; //look in the opposite view
                
                double n = shared[i]->GetNextEnergy(chain,chainhit,inview);
                double p = shared[i]->GetPreviousEnergy(chain,chainhit,inview);
                //if either is 0 then we are the end of the chain, so count the non zero one twice
                if(n==0 && p >0)n=p;
                if(p==0 && n>0)p=n;
                
                //printf("e %f\n",n+p);
                
                es.push_back(n+p);
        }
        
        double total=0;;
        double totaln=0;
        for(unsigned int i=0;i<es.size();i++)
        {
                totaln++;
                total+=es[i];
        }
        
        ChainHelper * ch = view == 2 ? &mypoh->chu : &mypoh->chv;
        
        Chain *c =ch->GetChain(chain);
        double energy = c->e[chainhit];
        
        //printf("Sharing %f\n",energy);
        
        if(total==0)return; //shared hits with no energy!
        for(unsigned int i=0;i<es.size();i++)
        {
                shared[i]->SetEnergy(chain,chainhit,view,energy*es[i]/total);
        }
        

}




std::vector<foundpath> Finder::GetPaths(ChainHelper*ch)
{


        std::vector<foundpath> p;
        for(unsigned int i=0;i<ch->finished.size();i++)
        {
                if(ch->finished[i].parentChain<0 )
                {
                        std::vector<foundpath> t = ch->FindPaths(&ch->finished[i]);
                        //printf("head chain %d has %d paths\n",ch->finished[i].myId, t.size());
                        
                        for(unsigned int j=0;j<t.size();j++)
                        p.push_back(t[j]);
        
                }
         }
         return p;
}

void Finder::DumpPaths(std::vector<foundpath> a)
{

        ostringstream os;
        os << "Dumping Paths\n";
        for(unsigned int i=0;i<a.size();i++)
        {
                os << "Path "<<i<<" :  ";
                os << " start t,z "<<a[i].start_t<<", "<<a[i].start_z;
                os << " end t,z "<<a[i].end_t<<", "<<a[i].end_z;
                os << " entries " <<a[i].entries << " energy " << a[i].energy << " muonlike " << a[i].muonlike;
                os << "\n\t Chain Path : ";
                for(unsigned int j=0;j<a[i].path.size();j++)
                        os<<a[i].path[j]<<" ";
                os <<"\n\n";
        }
        
        MSG("Finder",Msg::kDebug) << os.str();  

}
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
        
        bool pointgreater(point p1, point p2){return p1.z < p2.z;}
        

Particle3D Finder::Make3DParticle(std::vector<int>upath, std::vector<int>vpath, ChainHelper * chu, ChainHelper * chv,int multu=0, int multv=0)
{
        Particle3D p;


        
        std::vector<point> points;
        
        double start =0;
        double end  =1000;
        
        double endu=0;
        double endv=0;
                
        
        
        
        int type=0;
        
        
        for(unsigned int i=0;i<upath.size();i++)
        {
                Chain *c = chu->GetChain(upath[i]);
                //cout <<"using chain "<<upath[i]<<"\n";
                                
                if(c->particletype)
                        type = c->particletype;
                
                for(int j=0;j<c->entries;j++)
                {
                        if(c->parentChain==-1 && j==0)
                                start = start < c->z[j]? c->z[j] : start;
                        if(i == upath.size()-1 && j == c->entries-1)
                        {
                                end = end < c->z[j] ? end : c->z[j];
                        }
                        endu=endu<c->z[j]?c->z[j]:endu;

                        //don't double count children head hits
                        if(c->parentChain>-1 && j==0)
                        {
                                if(i==0)continue;
                                Chain *cpast = chu->GetChain(upath[i-1]);
                                if(cpast->cluster_id[cpast->entries-1] == c->cluster_id[0])
                                        continue; 
                        }
                        Managed::ManagedCluster * clu = mypoh->clusterer.GetCluster(c->cluster_id[j]);  
                        if(!clu)continue;               
                
                        
                        point a;
                        a.z=c->z[j];
                        a.t=c->t[j];
                        if(c->available && clu->id>0)
                                a.e=clu->e;
                        else a.e=0;
                        a.chain=c->myId;
                        if(c->available && clu->id>0)
                                a.chainhit=j;
                        else
                                a.chainhit=-1;
                        a.view=2;
                        

                        if(clu)
                        {
                                a.rms_t=clu->rms_t;
                                points.push_back(a);
                        }
                        
                        
                //      printf("adding pt u %f %f %f\n",a.z,a.t,a.e);
                        
                }
        
        }
        for(int i=0;i<(int)vpath.size();i++)
        {
                Chain *c = chv->GetChain(vpath[i]);
                
                //cout <<"using chain "<<vpath[i]<<"\n";
                
                if(c->particletype)
                        type = c->particletype;
                for(int j=0;j<c->entries;j++)
                {
                        if(c->parentChain==-1 && j==0)
                                start = start < c->z[j]? c->z[j] : start;
                        if(i == (int)vpath.size()-1 && j == c->entries-1)
                        {
                                end = end < c->z[j] ? end : c->z[j];
                        }
                        endv=endv<c->z[j]?c->z[j]:endv;
                        
                        //don't double count children head hits
                        if(c->parentChain>-1 && j==0)
                        {
                                if(i==0)continue;
                                Chain *cpast = chv->GetChain(vpath[i-1]);
                                if(cpast->cluster_id[cpast->entries-1] == c->cluster_id[0])
                                        continue; 
                        }

                        Managed::ManagedCluster * clu = mypoh->clusterer.GetCluster(c->cluster_id[j]);  
                        if(!clu)continue;               
                        

                        point a;
                        a.z=c->z[j];
                        a.t=c->t[j];
                        if(c->available && clu->id>0)
                                a.e=clu->e;
                        else
                                a.e=0;
                        a.chain=c->myId;
                        if(c->available && clu->id>0)
                                a.chainhit=j;
                        else
                                a.chainhit=-1;
                        a.view=3;
                        

                        if(clu){
                                a.rms_t=clu->rms_t;
                                points.push_back(a);
                        }

                //      printf("adding pt v %f %f %f\n",a.z,a.t,a.e);

                }
        
        }
        
        //we don't want the particle to extrapolate too far....
        double stopend = endu<endv?endu:endv;
        //but if that chain path is only used once, then we want to go all the way to the end!
        if(!multu && multv)stopend = stopend<endu?endu:stopend;
        if(!multv && multu)stopend = stopend<endv?endv:stopend;
        if(!multv && !multu)stopend = endu<endv?endv:endu;
                
        std::sort(points.begin(), points.end(),pointgreater);

/*
double largest_z=0;
double largest_idx=0;
        for(unsigned int i=0;i<points.size();i++)
        {
                if(points[i].z>largest_z)
                {
                        largest_z=points[i].z;
                        largest_idx=i;
                }
        }
*/      
        ///printf("LARGEST Z %f t %f e %f\n",largest_z,points[largest_idx].t,points[largest_idx].e);

        
        for(unsigned int i=0;i<points.size();i++)
        {
//              if( start - points[i].z > 0.045 )continue;
//              if( points[i].z - end > 0.045)continue; //within 1 planes, we want the end of the chain
        

        
        
        
        //are we at the end of usuable points?
        int done=0;
        for(unsigned int j=i;j<points.size();j++)
        {
                if(points[j].e>0)break;
        
                if(j==points.size()-1)done=1;
        }
        if(done)break;
        
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
                        
                //      printf("looking for lower value\n");
                        //do we have another lower value?
                        int lower2=-1;
                        for(int j=lower-1;j>-1;j--)
                        {
                                if(points[j].view!=myview && fabs(points[lower].z-points[j].z)>0.04)
                                {
                                        lower2=j;
                //                      printf("second lower found\n");
                                        break;
                                }
                        }
                        
                        
                        if(lower2>-1 && fabs( points[lower].z - points[lower2].z) >0)
                        {
                                double s = (points[lower].t - points[lower2].t) /  ( points[lower].z - points[lower2].z);
                        
                                double off = points[lower].t - points[lower].z*s;
                        
                                double t = s * points[i].z + off;
                                u = myview == 2 ? u : t;
                                v = myview == 3 ? v : t;
                //              printf("extrap  tz tz   %f %f    %f %f to %f\n",points[lower].t,points[lower].z,points[lower2].t,points[lower2].z,points[i].z);
                        }else{
                                u = myview == 2 ? u : points[lower].t;
                                v = myview == 3 ? v : points[lower].t;
        //                      printf("using closer point\n");
                        }
        //              printf("uv %f %f\n",u,v);
                        
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
                                double off = points[upper].t - points[upper].z*s;
                        
                                double t = s * points[i].z + off;
                                u = myview == 2 ? u : t;
                                v = myview == 3 ? v : t;
                                                                
                //              printf("ext next hit with t %f z %f  t %f z %f for ext\n",points[upper].t,points[upper].z,points[upper2].t,points[upper2].z);
                        
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
*/                              double vz =vtx_z;
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
                                
                        //      printf("---view %d u %f v %f t %f\n",myview,u,v,t);
                        //      printf("---vtx z %f t%f upper z %f t %f\n",vz,vt,points[upper].z,points[upper].t);
                        

                                u = myview == 2 ? u : t;
                                v = myview == 3 ? v : t;
//                      }       
                                

                }
                else if(upper==-1 && lower==-1) //we have an empty view!!!
                {

                        u = myview == 2 ? u : 0;
                        v = myview == 3 ? v : 0;
                }
                
                p.add_to_back(u,v,z,e,points[i].chain,points[i].chainhit,view,rms_t);
                
                //printf("adding %f %f %f --- %f\n",u,v,z,e);
        }
        
        
        if(type)
                p.particletype=(Particle3D::EParticle3DType)type;
        
        p.finalize();
        
        return p;
}






void Finder::DumpParticle3D(std::vector<Particle3D*> p3d)
{
        ostringstream s;
        
        s << "Dump of Particle3D's "<<endl;
        s << "Number of Particle3D : "<<p3d.size()<< endl<<endl;
                
                
        for (unsigned int i=0;i<p3d.size();i++)
        {
                //ostringstream s;
        
                Particle3D * p = p3d[i];
                if (p==0)continue;
        
        

        
                s << "\n---Particle3d " << i << "---\nstart, stop (u,v,z) (" << p->start_u << ", "<< p->start_v << ", " << p->start_z <<") (" <<p->end_u<<", "<< p->end_v << ", " << p->end_z <<")"<<endl;
                s << "entries "<< p->entries << " muonlike " << p->muonfrac <<endl;
                
                

                
                s<<"types: ";
                for(unsigned int j=0;j<p->types.size();j++)
                {
                        switch( p->types[j].type)
                        {
                                case    ParticleType::em:
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
                s<<"vise "<<p->sum_e<<endl;
                                
                s << "points (u,v,z,e - chain, chainhit, chainview - shared - rms_t - view) : ";
                for(int j=0;j<p->entries;j++)
                        s << "(" << p->u[j]<<", "<<p->v[j]<<", "<<p->z[j]<<", "<<p->e[j]<<" - " << p->chain[j] <<", "<<p->chainhit[j]<<", "<<p->view[j]<<" - "<<p->shared[j]<<" - "<< p->rms_t[j]<< " - " << p->view[j]<<") ";
        
        /*
        s << "points (e - shared) : ";
                for(int j=0;j<p->entries;j++)
                        s << "(" <<p->e[j]<<" - "<<p->shared[j]<<") ";
        
        */
        
        

                        
                        

                        
                        
                s<<endl<<endl;
        
                        //check to see if it is a good particle........
                for(unsigned int q=0;q<p->u.size();q++)
                {
                  if (p->e[q]<0 ) { MSG("Finder",Msg::kDebug) <<"BAD E "<<p->e[q]<<" at "<<q<<"\n"; }
                  if (p->z[q]<0  ||p->z[q]>30 || isnan(p->z[q])) { MSG("Finder",Msg::kDebug) <<"BAD Z "<<p->z[q]<<" at "<<q<<"\n"; }
                  
                  if (p->u[q]<-10||p->u[q]>10 || isnan(p->u[q])) { MSG("Finder",Msg::kDebug) <<"BAD U "<<p->u[q]<<" at "<<q<<"\n"; }
                  if (p->v[q]<-10||p->v[q]>10 || isnan(p->v[q])) { MSG("Finder",Msg::kDebug) <<"BAD U "<<p->v[q]<<" at "<<q<<"\n"; }
                }

        
                
                
        }
        
        s << "\n!!! shareholder has "<<shareholder.GetTotRemaining() << " hits still shared!\n";


        MSG("Finder",Msg::kDebug) <<s.str();

        

}

/*
enum particletype
{
        em,
        emshort,
        muon,
        prot,
        pi0,
        uniquemuon


};
*/

/*
struct ptype
{
        particletype type;
        int start;
        int stop;
};
*/


std::vector<Particle3D*> Finder::ProcessParticle3D(Particle3D * p3d)
{

        std::vector<Particle3D*> pout;

        ostringstream s;
        s<<"\nparticle types:\n";

        //std::vector<ptype> types;

        //for the muon checks
        double lowm=0.2;
        double highm=3;
        int minm=2;


        //em check
        int lastemend=-1;
        int lastemlow=-1;
        if(p3d->entries>3)
        for(int i=0;i<p3d->entries-2;i++)
        {
                double p0 = p3d->e[i];
                double p1 = p3d->e[i+1];
                double p2 = p3d->e[i+2];
                if(p0<p1 && p1>p2 && p1 > highm) 
                {

                        int low=i;
                        int high=i+2;

                        
                        for(int h=i-1;h>-1;h--)
                        {
                                if(p3d->e[h]<p3d->e[h+1])
                                {
                                        low=h;
                                }else{
                                        break;
                                }
                        }
                        for(int h=i+3;h<p3d->entries;h++)
                        {
                                if(p3d->e[h]<p3d->e[h-1])
                                {
                                        high=h;
                                }else{
                                        break;
                                }
                        }
                        if(low==lastemend)
                        {

                                s<<"pi0 like - found consecutive em like from "<<lastemlow<<" to "<<high<<endl;
                                ParticleType a;
                                a.type = ParticleType::pi0;
                                a.start=lastemlow;
                                a.stop=high;
                                p3d->types.push_back(a);
                        }
                        s<<"em like from "<< low << " to " <<high<<endl;
                        i=high; //start after this found em like structure
                        lastemend=high;
                        lastemlow=low;
                        ParticleType a;
                        a.type = ParticleType::em;
                        a.start=low;
                        a.stop=high;
                        p3d->types.push_back(a);        
                }               
        }
        
        //shortem check
        if(p3d->entries>3)
        {
                double p0 = p3d->e[0];
                double p1 = p3d->e[1];
                double p2 = p3d->e[2];  
                if(p0>p1 && p1>p2 && p0 > highm)
                {

                        s<<"shortem like"<<endl;
                        
                        int high=2;
                        for(int i=2;i<p3d->entries-1;i++)
                        {
                                if(p3d->e[i] > p3d->e[i+1])
                                {
                                        high=i+1;
                                }       
                        }
                        ParticleType a;
                        a.type = ParticleType::emshort;
                        a.start=0;
                        a.stop=high;
                        p3d->types.push_back(a);        
                }
        }       
        
        //muon check
        if(p3d->entries>1)
        for(int i=0;i<p3d->entries;i++)
        {
                int low=i;
                int high=i;
                if(p3d->e[i] > lowm && p3d->e[i] < highm)
                {
                        for(int j=i;j<p3d->entries;j++)
                        {
                                if(p3d->e[j] > lowm && p3d->e[j] < highm)
                                {
                                        high=j;
                                }else{
                                        if(j < p3d->entries-1 && high == j-1 && p3d->e[j] < highm*2 && p3d->e[j+1] > lowm && p3d->e[j+1] < highm)  //allow 1 out of range hit
                                        {
                                                high=j+1;
                                                j++;
                                                continue;
                                        }
                                        i=+1;
                                        break;
                                
                                }
                        }
                }
                
                if(high-low+1>=minm)
                {

                        s<<"muon like from "<<low<<" to "<<high<<endl;
                        ParticleType a;
                        a.type = ParticleType::muon;
                        a.start=low;
                        a.stop=high;
                        p3d->types.push_back(a);        
                }
                
                i=high+1;
        }

        //unique muon check --- should we be extracting a muon from a particle that has larger shared hits?
        double  mc=0;
        double me=0;
        int c=0;
        for(int i=0;i<p3d->entries;i++)
        {

                if(p3d->shared[i]==0)
                {
                        c++;
                        if(p3d->e[i] > lowm && p3d->e[i] < highm)
                        {
                                mc++;
                                me+=p3d->e[i];
                        }
                }
        }       

        if(mc/c>0.7)
        {
                s<<"unique muon found"<<endl;
                ParticleType a;
                a.type = ParticleType::uniquemuon;
                a.start=-1;
                a.stop=-1;
                p3d->types.push_back(a);        
        }
        //////////



        
        //prot check // really stopping of +muon
        int plow=p3d->entries-1;
        int phigh=p3d->entries-1;
        if(p3d->e[p3d->entries-1]>highm)
        for(int i=p3d->entries-2;i>-1;i--)
        {
                if(p3d->e[i]<p3d->e[i+1]) plow=i;
                else break;
        }
        if(phigh-plow>1)
        {

                s<<"proton/stopping u+ found from "<<plow<<" to "<<phigh<<endl;
                ParticleType a;
                a.type = ParticleType::prot;
                a.start=plow;
                a.stop=phigh;
                p3d->types.push_back(a);        
        }
        
        //we now have a list of types.... go through, and split up the Particle3D as needed until all Particle3Ds have a single type
        
        
        int emt=0;
        
        p3d->particletype=Particle3D::other;    
        for(unsigned int i=0;i<p3d->types.size();i++)
        {
        
        
                if(p3d->types[i].type==ParticleType::uniquemuon)
                {
                        std::pair<Particle3D*,Particle3D*> a = StripMuon(p3d);
                        if(a.first)
                        {
                                pout.push_back(a.first); //thats the stripped muon...
                                if(a.second) //make sure we removed something from p3d...
                                {
                                        std::vector<Particle3D*> b =ProcessParticle3D(a.second); //the adjusted remnant needs to be rechecked!
                                        for(unsigned int j=0;j<b.size();j++)
                                                pout.push_back(b[j]);
                                }
                                return pout;
                        }       
                }
                
        
                //see if it is definitely a muon!
                
                if(p3d->types[i].type==ParticleType::muon && p3d->types[i].stop-p3d->types[i].start+1 == p3d->entries)
                {
                        p3d->particletype=Particle3D::muon;
                        break;
                
                }

                
                if(p3d->types[i].type==ParticleType::em && p3d->types[i].stop-p3d->types[i].start+1 == p3d->entries)
                {
                        p3d->particletype=Particle3D::electron;
                        break;
                }

                if(p3d->types[i].type==ParticleType::em && (double)( p3d->types[i].stop-p3d->types[i].start) / (double) p3d->entries > 0.8) //more than half of it is em....
                {
                        p3d->particletype=Particle3D::electron;
                        break;
                }

                if(p3d->types[i].type==ParticleType::em )
                {
                        emt+= p3d->types[i].stop-p3d->types[i].start+1;
                }
        

        
        }

        if(p3d->particletype==0 && p3d->muonfrac>0.6)
        {
                p3d->particletype=Particle3D::muon;
                        
        }       

        if(p3d->particletype==0 && (double)emt / (double)p3d->entries > 40)
        {
                p3d->particletype=Particle3D::electron;
                        
        }       


        MSG("Finder",Msg::kDebug) << s.str();   
        
        pout.push_back(p3d);
        return pout;
        
}


std::pair<Particle3D*,Particle3D*> Finder::StripMuon(Particle3D * p3d)
{
        std::pair<Particle3D*,Particle3D*>  a;
        a.first=0;
        a.second=p3d;
        
        //see if all unique hits are muons hits....
        double lowm=0.0;
        double highm=3.5;
        
        int u=0;
        int m=0;
        double me=0;
        
        for(int i=0;i<p3d->entries;i++)
        {
                if(p3d->shared[i]==0)
                {
                        u++;
                        if(p3d->e[i] > lowm && p3d->e[i]<highm)
                        {
                                m++;
                                me+=p3d->e[i];
                        }
                }
        }
        
        if(u==m)
        { 
                a.second=0;  //releasing shared particles which cant be part of something else...
                //cout <<"!!!---releasing particle"<<endl;
        }
        
        
        Particle3D * c = new Particle3D();
        a.first=c;
        
        
        double eavg = me/m;
        for(int i=0;i<p3d->entries;i++)
        {
                double e = p3d->e[i] > lowm && p3d->e[i]<highm ? p3d->e[i] : eavg;
                e = e > p3d->e[i] ? p3d->e[i] : e;
        
                c->add_to_back( p3d->u[i],  p3d->v[i],  p3d->z[i], e, p3d->chain[i], p3d->chainhit[i],  p3d->view[i],0);
                c->lockshared[c->entries-1]=1;
                
                if(a.second)
                {
                        double reste=p3d->e[i]-e>0 ? p3d->e[i]-e : 0;
                        a.second->e[i]= reste;
                        
                        
                }
        }
        c->particletype=Particle3D::muon;



        if(a.second)
        {
                double sume=0;
                int num0=0;
                for(int i=0;i<a.second->entries;i++)
                {
                        sume+=a.second->e[i];
                        if(a.second->e[i]<0.01)num0++;
                }
                //cout<<"second e"<<sume<<endl;
                
                if(sume<0.0001 || (double)num0/(double)a.second->entries > 0.8)
                {
                        a.second=0;  //releasing shared particles which cant be part of something else...
                        //cout <<"!!!---releasing particle"<<endl;
                }
        }

        return a;

}






std::vector<Particle3D*> Finder::ProcessParticle3D1(Particle3D * p3d, int check_unique_muon)
{

        std::vector<Particle3D*> pout;

        ostringstream s;
        s<<"\nparticle types:\n";

        //std::vector<ptype> types;

        //for the muon checks
        double lowm=0.2;
        double highm=3;
        int minm=2;


        //em check
        int lastemend=-1;
        int lastemlow=-1;
        if(p3d->entries>3)
        for(int i=0;i<p3d->entries-2;i++)
        {
                double p0 = p3d->e[i];
                double p1 = p3d->e[i+1];
                double p2 = p3d->e[i+2];
                if(p0<p1 && p1>p2 && p1 > highm) 
                {

                        int low=i;
                        int high=i+2;

                        
                        for(int h=i-1;h>-1;h--)
                        {
                                if(p3d->e[h]<p3d->e[h+1])
                                {
                                        low=h;
                                }else{
                                        break;
                                }
                        }
                        for(int h=i+3;h<p3d->entries;h++)
                        {
                                if(p3d->e[h]<p3d->e[h-1])
                                {
                                        high=h;
                                }else{
                                        break;
                                }
                        }
                        if(low==lastemend)
                        {

                                s<<"pi0 like - found consecutive em like from "<<lastemlow<<" to "<<high<<endl;
                                ParticleType a;
                                a.type = ParticleType::pi0;
                                a.start=lastemlow;
                                a.stop=high;
                                p3d->types.push_back(a);
                        }
                        s<<"em like from "<< low << " to " <<high<<endl;
                        i=high; //start after this found em like structure
                        lastemend=high;
                        lastemlow=low;
                        ParticleType a;
                        a.type = ParticleType::em;
                        a.start=low;
                        a.stop=high;
                        p3d->types.push_back(a);        
                }               
        }
        
        //shortem check
        if(p3d->entries>3)
        {
                double p0 = p3d->e[0];
                double p1 = p3d->e[1];
                double p2 = p3d->e[2];  
                if(p0>p1 && p1>p2 && p0 > highm)
                {

                        s<<"shortem like"<<endl;
                        
                        int high=2;
                        for(int i=2;i<p3d->entries-1;i++)
                        {
                                if(p3d->e[i] > p3d->e[i+1])
                                {
                                        high=i+1;
                                }       
                        }
                        ParticleType a;
                        a.type = ParticleType::emshort;
                        a.start=0;
                        a.stop=high;
                        p3d->types.push_back(a);        
                }
        }       
        
        //muon check
        if(p3d->entries>1)
        for(int i=0;i<p3d->entries;i++)
        {
                int low=i;
                int high=i;
                if(p3d->e[i] > lowm && p3d->e[i] < highm)
                {
                        for(int j=i;j<p3d->entries;j++)
                        {
                                if(p3d->e[j] > lowm && p3d->e[j] < highm)
                                {
                                        high=j;
                                }else{
                                        if(j < p3d->entries-1 && high == j-1 && p3d->e[j] < highm*2 && p3d->e[j+1] > lowm && p3d->e[j+1] < highm)  //allow 1 out of range hit
                                        {
                                                high=j+1;
                                                j++;
                                                continue;
                                        }
                                        i=+1;
                                        break;
                                
                                }
                        }
                }
                
                if(high-low+1>=minm)
                {

                        s<<"muon like from "<<low<<" to "<<high<<endl;
                        ParticleType a;
                        a.type = ParticleType::muon;
                        a.start=low;
                        a.stop=high;
                        p3d->types.push_back(a);        
                }
                
                i=high+1;
        }

        //unique muon check --- should we be extracting a muon from a particle that has larger shared hits?

        if(check_unique_muon)
        {
        double  mc=0;
        double me=0;
        int c=0;
        for(int i=0;i<p3d->entries;i++)
        {

//              if(p3d->shared[i]==0)
//              {
                        c++;
                        if(p3d->e[i] > lowm && p3d->e[i] < highm)
                        {
                                mc++;
                                me+=p3d->e[i];
                        }
//              }
        }       

        if(mc>0.7*c)
        {
                s<<"unique muon found muonlike to all hits ratio "<<mc/c<<endl;
                ParticleType a;
                a.type = ParticleType::uniquemuon;
                a.start=-1;
                a.stop=-1;
                p3d->types.push_back(a);        
        }
        }
        //////////



        
        //prot check // can also be a  stopping of +muon
        int plow=p3d->entries-1;
        int phigh=p3d->entries-1;
        if(p3d->e[p3d->entries-1]>highm)
        for(int i=p3d->entries-2;i>-1;i--)
        {
                if(p3d->e[i]<p3d->e[i+1]) plow=i;
                else break;
        }
        if(phigh-plow>0)
        {

                s<<"proton/stopping u+ found from "<<plow<<" to "<<phigh<<endl;
                ParticleType a;
                a.type = ParticleType::prot;
                a.start=plow;
                a.stop=phigh;
                p3d->types.push_back(a);        
        }
        
        //we now have a list of types.... go through, and split up the Particle3D as needed until all Particle3Ds have a single type
        
        
        int emt=0;
        
        //p3d->particletype=0;
        p3d->particletype=Particle3D::other;
        
        
        
        std::vector<ParticleType>::iterator it;
        
        
        //cout <<"num of types for this particle "<<  p3d->types.size()<<endl;
                
        for(it=p3d->types.begin();it!=p3d->types.end();it++)
        {
        
                
                if(it->type==ParticleType::uniquemuon)
                {
                        ostringstream s;
                        s<<"ATTEMPTING TO STRIP MUON\n";
                        std::pair<Particle3D*,Particle3D*> a = StripMuon1(p3d);
                        if(a.first)
                        {
                                s << "STRIPPED MUON\n";
                                pout.push_back(a.first); //thats the stripped muon...
                                if(a.second) //make sure we removed something from p3d...
                                {
                                        s <<"  OTHERS REMAIN!\n";
                                
                                        std::vector<Particle3D*> b =ProcessParticle3D1(a.second); //the adjusted remnant needs to be rechecked!
                                        for(unsigned int j=0;j<b.size();j++)
                                                pout.push_back(b[j]);
                                }
                                //return pout;
                                break;
                        }else{  //remove the incorrectly assigned unique muon tag
                                p3d->types.erase(it);
                        
                        }       
                }
                
                MSG("Finder",Msg::kDebug)<<s;
        
        
                //see if it is definitely a muon!
                
                if(it->type==ParticleType::muon && it->stop - it->start+1 == p3d->entries)
                {
                        p3d->particletype=Particle3D::muon;
                //      break;
                
                }

                
                if(it->type==ParticleType::em && it->stop - it->start+1 == p3d->entries)
                {
                        p3d->particletype=Particle3D::electron;
                //      break;
                }

                if(it->type==ParticleType::em && (double)( it->stop - it->start) / (double) p3d->entries > 0.8) //more than half of it is em....
                {
                        p3d->particletype=Particle3D::electron;
                //      break;
                }

                if(it->type==ParticleType::em )
                {
                        emt+= it->stop - it->start+1;
                }
        

        
        }

        if(p3d->particletype==0 && p3d->muonfrac>0.6)
        {
                p3d->particletype=Particle3D::muon;
                        
        }       

        if(p3d->particletype==0 && (double)emt / (double)p3d->entries > 40)
        {
                p3d->particletype=Particle3D::electron;
                        
        }       


        MSG("Finder",Msg::kDebug) << s.str()<<endl;     
        
        pout.push_back(p3d);
        return pout;
        
}



std::pair<Particle3D*,Particle3D*> Finder::StripMuon1(Particle3D * p3d)
{
        std::pair<Particle3D*,Particle3D*>  a;
        a.first=0;
        a.second=p3d;
        
        //see if all unique hits are muons hits....
        double lowm=0.0;
        double highm=3.5;
        
        int u=0;
        int m=0;
        double me=0;
        
        for(int i=0;i<p3d->entries;i++)
        {
                if(p3d->shared[i]==0)
                {
                        u++;
                        if(p3d->e[i] > lowm && p3d->e[i]<highm)
                        {
                                m++;
                                me+=p3d->e[i];
                        }
                }
        }
        
        if(u==m)
        { 
                a.second=0;  //releasing shared particles which cant be part of something else...
                //cout <<"!!!---releasing particle"<<endl;
        }
        
        
        Particle3D * c = new Particle3D();
        a.first=c;
        
        
        double eavg = me/m;
        for(int i=0;i<p3d->entries;i++)
        {
                double e = p3d->e[i] > lowm && p3d->e[i]<highm ? p3d->e[i] : eavg;
                e = e > p3d->e[i] ? p3d->e[i] : e;
        
                c->add_to_back( p3d->u[i],  p3d->v[i],  p3d->z[i], e, p3d->chain[i], p3d->chainhit[i],  p3d->view[i],0);
                c->lockshared[c->entries-1]=1;

                shareholder.Take(p3d->chain[i],p3d->chainhit[i],e);
                
                if(a.second)
                {
                        double reste=p3d->e[i]-e>0 ? p3d->e[i]-e : 0;
                        a.second->e[i]= reste;
                        
                        shareholder.Take(p3d->chain[i],p3d->chainhit[i],e);
                        
                        
                }
        }
        c->particletype=Particle3D::muon;

        ParticleType pta;
        pta.type = ParticleType::muon;
        pta.start=0;
        pta.stop=p3d->entries;
        c->types.push_back(pta);        


        if(a.second)
        {
                double sume=0;
                int num0=0;
                for(int i=0;i<a.second->entries;i++)
                {
                        sume+=a.second->e[i];
                        if(a.second->e[i]<0.01)num0++;
                }
                MSG("Finder",Msg::kDebug)<<"second e"<<sume<<endl;
                
                if(sume<0.0001 || (double)num0/(double)a.second->entries > 0.8)
                {
                        a.second=0;  //releasing shared particles which cant be part of something else...
                        //cout <<"!!!---releasing particle"<<endl;
                }
                else  //we are keeping the particle....
                {


                        //remove the unique muon indicator....
                
                        std::vector<ParticleType>::iterator it;
                        for(it=a.second->types.begin();it!=a.second->types.end();it++)
                        {
                                if(it->type==ParticleType::uniquemuon)
                                {
                                        a.second->types.erase(it);
                                        it=a.second->types.begin(); //for now restart at the beginning... not the best way
                                        MSG("Finder",Msg::kDebug)<<"removed uniquemuon indicator from particle3d "<<endl;
                                }
                        
                        }


                        //clear particle list....
                        a.second->types.clear();

                }
        }


        return a;

}


bool p3dsharedsort(Particle3D* a, Particle3D* b)
{
        return a->numshared < b->numshared; 
}


void Finder::RemoveSharedHits(std::vector<Particle3D*> pv)
{

        MSG("Finder",Msg::kDebug) <<"\ncleaning list\n";
        for (unsigned int i=0;i<pv.size();i++)
        {
                        MSG("Finder",Msg::kDebug) << " removing with "<<pv[i]->numshared<<" shared \n";
        
                        //first clear out any hits marked as shared that are not!
                        for(int j=0;j<pv[i]->entries;j++)
                        {
                                if(pv[i]->shared[j])
                                {
                                        if(shareholder.GetNumShared(pv[i]->chain[j],pv[i]->chainhit[j])<2)
                                                pv[i]->UnsetShared(pv[i]->chain[j],pv[i]->chainhit[j]);
                                        if(shareholder.GetNumShared(pv[i]->chain[j],pv[i]->chainhit[j])==1)
                                        {
                                                int c=pv[i]->chain[j];
                                                int ch=pv[i]->chainhit[j];
                                                double e = shareholder.GetEShared(c,ch);
                                        
                                        
                                                pv[i]->e[j]=e;
                                                shareholder.Take(c,ch,e);
                                        }
                                }
                        
                        }
        
        }


        sort(pv.begin(), pv.end(), p3dsharedsort);

        MSG("Finder",Msg::kDebug) <<"\nremoving hits\n";

        for (unsigned int i=0;i<pv.size();i++)
        {               
                        shareholder.dump();

                        //first clear out any hits marked as shared that are not!
                        for(int j=0;j<pv[i]->entries;j++)
                        {
                                if(pv[i]->shared[j])
                                {
                                        if(shareholder.GetNumShared(pv[i]->chain[j],pv[i]->chainhit[j])<2)
                                                pv[i]->UnsetShared(pv[i]->chain[j],pv[i]->chainhit[j]); 
                                        if(shareholder.GetNumShared(pv[i]->chain[j],pv[i]->chainhit[j])==1)
                                        {
                                                int c=pv[i]->chain[j];
                                                int ch=pv[i]->chainhit[j];
                                                double e = shareholder.GetEShared(c,ch);
                                        
                                        
                                                pv[i]->e[j]=e;
                                                shareholder.Take(c,ch,e);
                                        }               
                                }
                        
                        }





                if (pv[i]->particletype==Particle3D::muon)
                {
                        MSG("Finder",Msg::kDebug)<<"!!muon\n";
                }
                else if (pv[i]->particletype==Particle3D::electron)
                {
                        MSG("Finder",Msg::kDebug)<<"!!elec\n";
                }
                else MSG("Finder",Msg::kDebug)<<"!!other "<<pv[i]->particletype<<"\n";
                
                
                //for now, treat all types similary.... 
                //first look for shared hit surrounded by two unshared hits, and set it to the average
                for(int j=1;j<pv[i]->entries-2;j++)
                {
                        if(pv[i]->shared[j] && !pv[i]->shared[j-1] && !pv[i]->shared[j+1])
                        {
                                int c=pv[i]->chain[j];
                                int ch=pv[i]->chainhit[j];
                                double e = shareholder.GetEShared(c,ch);
                                        
                                double totake=(pv[i]->e[j-1]+pv[i]->e[j+1])/2;
                                totake=totake>e?e:totake;
                                
                                pv[i]->e[j]=totake;
                                shareholder.Take(c,ch,totake);
                                pv[i]->UnsetShared(c,ch);
                        
                        }
                }
                
                //if that doesn't work, try the average of closest unshared hits
                for(int j=0;j<pv[i]->entries;j++)
                {
                        int foundupper=-1;
                        int foundlower=-1;
                        for(int k=j+1;k<pv[i]->entries;k++)
                        {
                                //printf("%d(%d) %d(%d)\n",j,pv[i]->shared[j],k,pv[i]->shared[k]);
                                if(pv[i]->shared[j]&&! pv[i]->shared[k])
                                {
                                        foundupper=k;
                                        break;
                                }
                        }
                        
                        for(int k=j-1;k>-1;k--)
                        {
                                if(pv[i]->shared[j]&&! pv[i]->shared[k])
                                {
                                        foundlower=k;
                                        break;
                                }
                        }
                        
                        double totake=0;
                        
                        if(foundupper>-1 && foundlower>-1)
                        {
                                totake = (pv[i]->e[foundupper]+pv[i]->e[foundlower])/2;
                        }
                        else if(foundupper>-1)
                        {
                                totake = (pv[i]->e[foundupper]);
                        }
                        else if(foundlower>-1)
                        {
                                totake = (pv[i]->e[foundlower]);
                        }else continue;
                        
                        int c=pv[i]->chain[j];
                        int ch=pv[i]->chainhit[j];
                        double e = shareholder.GetEShared(c,ch);        
                        totake=totake>e?e:totake;
                        pv[i]->e[j]=totake;
                        MSG("Finder",Msg::kDebug)<<"taking from "<<c<<" "<<ch << " with shard "<<shareholder.GetNumShared(c,ch)<<"\n";
                        shareholder.Take(c,ch,totake);
                        MSG("Finder",Msg::kDebug)<<"now has "<<shareholder.GetNumShared(c,ch)<<" shared\n";
                        pv[i]->UnsetShared(c,ch);
                        
                }
                
                pv[i]->finalize();
                
        }


}
