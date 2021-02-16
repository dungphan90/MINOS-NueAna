#include "NueAna/ParticlePID/ParticleFinder/Display/HoughView.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterManager.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedCluster.h"
#include "TPolyLine.h"
#include "TMarker.h"
#include "TArrow.h"
#include "TMath.h"

#include "Conventions/SimFlag.h"
#include <vector>
#include "MCNtuple/NtpMCStdHep.h"
#include "MCNtuple/NtpMCRecord.h"
#include "MCNtuple/NtpMCTruth.h"
#include <math.h>
#include "TruthHelperNtuple/NtpTHEvent.h"

#include <map>

#include "NueAna/ParticlePID/ParticleFinder/HoughLine.h"


using namespace std;

HoughView::HoughView():Page()
{
	viewU=0;
	viewV=0;
	padU=0;
	padV=0;
	houghviewU=0;
	houghviewV=0;

}


HoughView::~HoughView()
{}


void HoughView::BuildDisplay(TCanvas *c)
{
	myCanvas=c;
	myCanvas->cd();
	padU=new TPad("houghview_padu","houghview_padu",0,0.5,0.7,1.0);
	padV=new TPad("houghview_padv","houghview_padv",0,0,0.7,0.5);


	padU->Draw();
	padV->Draw();
	

	
}


void HoughView::DrawEvent(ParticleObjectHolder * poh, NtpStRecord * /*ntp*/,int /*ntpEvt*/)
{
	myCanvas->cd();

	padU->Clear();
	padV->Clear();

	padU->Divide(2,1);
	padV->Divide(2,1);

	mypoh=poh;
	
	int isMC=0;
	if(poh->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC)isMC=1;
	
	double minz= mypoh->event.minz;
	double minu= mypoh->event.minu;
	double minv= mypoh->event.minv;
	double maxz= mypoh->event.maxz;
	double maxu= mypoh->event.maxu;
	double maxv= mypoh->event.maxv;
	

	
	if(isMC)
	{
		if(mypoh->mctrue.vtx_z>=0)minz=minz<mypoh->mctrue.vtx_z ? minz : mypoh->mctrue.vtx_z;
		maxz=maxz>mypoh->mctrue.vtx_z ? maxz : mypoh->mctrue.vtx_z;
		minu=minu<mypoh->mctrue.vtx_u ? minu : mypoh->mctrue.vtx_u;
		maxu=maxu>mypoh->mctrue.vtx_u ? maxu : mypoh->mctrue.vtx_u;
		minv=minv<mypoh->mctrue.vtx_v ? minv : mypoh->mctrue.vtx_v;
		maxv=maxv>mypoh->mctrue.vtx_v ? maxv : mypoh->mctrue.vtx_v;
	}
	



	minz -= 0.1;
	minu -= 0.1;
	minv -= 0.1;
	maxz += 0.1;
	maxu += 0.1;
	maxv += 0.1;
	
	int nu=(int)((maxu-minu)/0.0412);
	int nv=(int)((maxv-minv)/0.0412);
	int nz=(int)((maxz-minz)/0.0354);


	if(viewU)delete viewU;
	if(viewV)delete viewV;

	if(houghviewU)delete houghviewU;
	if(houghviewV)delete houghviewV;

	viewU = new TH2D("hu","Clusters U", nz, minz,maxz,nu,minu, maxu);
	viewV = new TH2D("hv","Clusters V", nz, minz,maxz,nv,minv, maxv);
	houghviewU = new TH2D("hhu","Hough U", 80, 0,3.141592,50,-1,1);
	houghviewV = new TH2D("hhv","Hough V", 80, 0,3.141592,50,-1,1);
	viewU->SetStats(0);
	viewV->SetStats(0);
	viewU->SetDirectory(0);
	viewV->SetDirectory(0);

	houghviewU->SetStats(0);
	houghviewV->SetStats(0);
	houghviewU->SetDirectory(0);
	houghviewV->SetDirectory(0);

	if(mypoh)
	{	
	

		//psf.FindPrimaryShower(&mypoh->clusterer);
	
	
		padU->cd(1);	
		
		viewU->Draw();
		DrawClusters(2,viewU);
		DrawVertex(2);
		if(isMC)DrawTrueVertex(2);
	//	DrawChains(2);

		padU->cd(2);
		DrawHough(2,houghviewU);

		

		padV->cd(1);
		viewV->Draw();
		DrawClusters(3,viewV);
		DrawVertex(3);
		if(isMC)DrawTrueVertex(3);
	//	DrawChains(3);	
		padV->cd(2);
		DrawHough(3,houghviewV);




	}

	
	double maxz_u = viewU->GetMaximum();
	double maxz_v = viewV->GetMaximum();
	double maxz_both = maxz_u>maxz_v?maxz_u:maxz_v;
	viewU->GetZaxis()->SetRangeUser(0, maxz_both);
	viewV->GetZaxis()->SetRangeUser(0, maxz_both);
	
	padU->Update();
	padV->Update();

	myCanvas->Update();
	
}


void HoughView::DrawChains(int view)
{       
        
	ChainHelper * ch = 0;
    if(view==2)
    	ch=&mypoh->chu;
    else if(view==3)
    	ch=&mypoh->chv;
    else return;    


    std::vector<int> path =ch->FindMaxPath();

	for(unsigned int i=0;i<ch->finished.size();i++)
    {
	    if(ch->finished[i].parentChain>-1 || ch->finished[i].entries < 1)continue;
        int parentcount=2;
        Chain *todraw=ch->GetChain(ch->finished[i].myId);

		//draw the parent
	    int dsize = todraw->t.size();
        int needconnect=0;

 	    if (dsize==1)  //chain with only 1 cluster needs to draw line to parent explicitly... if not a head chain
                                                //other wise, make the single cluster, head chain, no children noticeable
         {
         //      if(todraw->parentChain<0)continue;
	         if(todraw->parentChain>-1)
    	     {
        	     dsize=2;
                 needconnect=1;
             }
         }

      	 double * tt = new double[todraw->t.size()];
         double * tz = new double[todraw->t.size()];
         if(needconnect)
         {
        	 tt[0]=ch->GetChain(todraw->parentChain)->end_t;
             tz[0]=ch->GetChain(todraw->parentChain)->end_z;
             tz[1]=todraw->z[0];
             tt[1]=todraw->t[0];
         }
         else
         {
     	     std::copy(todraw->t.begin(),todraw->t.end(),tt);
             std::copy(todraw->z.begin(),todraw->z.end(),tz);
		
             if(dsize==1)
             {
        	     tt[1]=tt[0];
                 tz[1]=tz[0]-0.03;
                 dsize=2;

             }

         }

         TPolyLine *pl = new TPolyLine(dsize,tz,tt );

         //int linecol =2+todraw->level;

        //      pl->SetLineColor(linecol);
        //      pl->SetLineColor(parentcount);
 		 pl->SetLineColor(2);

         int width=1;
         //should make this better if it is to be really used

         for(unsigned int j=0;j<path.size();j++)
         {
         	if(path[j]==todraw->myId)width+=1;
         }

         pl->SetLineWidth(width);
         //    pl->Draw();
         pl->Draw("same");
                        
         for (unsigned int k=0;k<ch->finished[i].children.size();k++)
         {
         	parentcount++;
                
            std::vector<int> c = ch->GetAllChildren(ch->GetChain(ch->finished[i].children[k]));
                        
            printf("parent %d children ", ch->finished[i].myId);
                       
            for(unsigned int j=0;j<c.size();j++)
            {               
                                
            	Chain *todraw=ch->GetChain(c[j]);

                printf("%d ",todraw->myId);

                
                int dsize = todraw->t.size();
                int needconnect=0;

                if (dsize==1)  //chain with only 1 cluster needs to draw line to parent explicitly... if not a head chain
                                                //other wise, make the single cluster, head chain, no children noticeable 
                {
                //      if(todraw->parentChain<0)continue;
                
                        if(todraw->parentChain>-1) 
                        {
                                dsize=2;
                                needconnect=1;
                        }
                }

                double * tt = new double[todraw->t.size()];
                double * tz = new double[todraw->t.size()];
                if(needconnect)
                {
                        tt[0]=ch->GetChain(todraw->parentChain)->end_t;
                        tz[0]=ch->GetChain(todraw->parentChain)->end_z;
                        tz[1]=todraw->z[0];
                        tt[1]=todraw->t[0];
                }
                else
                {
                        std::copy(todraw->t.begin(),todraw->t.end(),tt);
                        std::copy(todraw->z.begin(),todraw->z.end(),tz);

                        if(dsize==1)
                        {
                                tt[1]=tt[0];
                                tz[1]=tz[0]-0.03;
                                dsize=2;

                        }

                }

                TPolyLine *pl = new TPolyLine(dsize,tz,tt );

                //int linecol =2+todraw->level;

        //      pl->SetLineColor(linecol);
                pl->SetLineColor(parentcount);


                int width=1;
                //should make this better if it is to be really used

                for(unsigned int j=0;j<path.size();j++)
                {
                        if(path[j]==todraw->myId)width+=1;
                }

                pl->SetLineWidth(width);
        //      pl->Draw();
                pl->Draw("same");
 
                }printf("\n");
        }}      
                
}       





void HoughView::DrawClusters(int view, TH2 * h)
{

	 std::map<double, std::map<double, int>  >::iterator p_iterr;
     std::map<double, int>::iterator s_iterr;
	

	Managed::ClusterManager *cluh = &mypoh->clusterer;
/*
	     for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
         {
         for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
         {
         
         	
       		h->Fill( p_iterr->first,s_iterr->first, cluh->GetCluster(s_iterr->second)->e);  
         }
         }
		h->SetStats(0);
		h->Draw("colz");
*/

  //Managed::ManagedCluster *largest=0;
  //double last_plane=0;
    for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
    {	

        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			//Managed::ManagedCluster * clus = cluh->GetCluster(s_iterr->second);
			h->Fill( p_iterr->first,s_iterr->first, cluh->GetCluster(s_iterr->second)->e); 

	  	}
	}


	
	Managed::ClusterSaver *	cluhs = &mypoh->cluster_saver;
    for(p_iterr=cluhs->GetClusterMap(view)->begin();p_iterr!=cluhs->GetClusterMap(view)->end(); p_iterr++)
    {	

        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluhs->GetCluster(s_iterr->second);
			
			//printf("%f %f %d\n",clus->z,clus->t,clus->GetStatus());
			if(clus->GetStatus()==1)continue;//long muon clusters
			
			h->Fill( p_iterr->first,s_iterr->first, cluhs->GetCluster(s_iterr->second)->e); 

	  	}
	}
	h->SetStats(0);
	h->Draw("colz");       
	
	
         
/*
    for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
    {	
    	//new plane and a largest cluster in the previous plane?
		if(fabs(last_plane-p_iterr->first)>0.03 && largest)
		{
			TMarker *m=new TMarker(largest->z,largest->t, 25);
			m->Draw("same");
			largest=0;
			last_plane=p_iterr->first;
		}
	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluh->GetCluster(s_iterr->second);

			if(!largest)
			{
				largest=clus;
				continue;
			}
			if(largest->e<clus->e)largest=clus;
	  	}
	}
	if(fabs(last_plane-largest->z)>0.03 && largest)
	{
		TMarker *m=new TMarker(largest->z,largest->t, 25);
		m->Draw("same");
		largest=0;
	}    
*/

//to take all within xxx% of the max hit in the event!

//find max value!
double maxvale=0;
    for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
    {	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluh->GetCluster(s_iterr->second);
			if(clus->e>maxvale)
			{
				maxvale=clus->e;
			}
		}
	}


    for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
    {	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluh->GetCluster(s_iterr->second);

		//	if(clus->e>maxvale*0.1)
	//		{
						if(clus->GetStatus()==1)continue;//long muon clusters
		TMarker *m=new TMarker(clus->z,clus->t, 25);
		m->Draw("same");
			
		//	}
	  	}
	}
	
	
/*
    for(p_iterr=cluhs->GetClusterMap(view)->begin();p_iterr!=cluhs->GetClusterMap(view)->end(); p_iterr++)
    {	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluhs->GetCluster(s_iterr->second);

			if(clus->GetStatus()==1)continue;//long muon clusters
		TMarker *m=new TMarker(clus->z,clus->t, 27);
		m->Draw("same");
			

	  	}
	}
*/	
	
	

}

bool operator < (const HoughLine & left, const HoughLine & right)
{
	if(left.ncluster < right.ncluster)return true;
	return false;
}



void HoughView::DrawVertex(int view)
{

	double t=0;
	if(view==2)
		t=mypoh->event.vtx_u;
	else if(view==3)
		t=mypoh->event.vtx_v;
	else return;
		
	TMarker *m = new TMarker(mypoh->event.vtx_z, t, 3);
	m->Draw("same");

}

void HoughView::DrawTrueVertex(int view)
{

	double t=0;
	if(view==2)
		t=mypoh->mctrue.vtx_u;
	else if(view==3)
		t=mypoh->mctrue.vtx_v;
	else return;
				
	TMarker *m = new TMarker(mypoh->mctrue.vtx_z, t, 2);
	m->Draw("same");

}
  
  
void HoughView::DrawHough(int view, TH2 * /*h*/)
{

	if(view==2)padU->cd(2);
	if(view==3)padV->cd(2);
	
	double t=0;
	double z=0;
	
	if(!mypoh->psf.ran)return;
	int cp = mypoh->psf.FindFirstIntersection(view,t,z);
	if(cp)
	{
		TMarker *m=new TMarker(z,t, 29);
		m->Draw("same");
		printf("closest intersection z %f t %f\n",z,t);
	}



	TH2D * inth = view==2? mypoh->psf.intU:mypoh->psf.intV;
	if(inth)
	{
		inth->SetStats(0);
		inth->Draw("colz");
	}
	
/*	

    TH2D * hh = mypoh->psf.GetHoughMap(view);     	
	hh->SetStats(0);
	hh->Draw("colz");
*/

	if(view==2)padU->cd(1);
	if(view==3)padV->cd(1);




	//double view_xmin = view==2?viewU->GetXaxis()->GetXmin():viewV->GetXaxis()->GetXmin();
	//double view_xmax = view==2?viewU->GetXaxis()->GetXmax():viewV->GetXaxis()->GetXmax();
	

	
	


	std::vector<HoughLine>*houghlines = mypoh->psf.GetHoughLines(view);
	
	printf("# hough lines %d\n",(int)houghlines->size());
	
	int max =-1;
	if((*houghlines).size()>50)max = (*houghlines).size()-50;
	for(int i=(*houghlines).size()-1;i>max;i--)
	{
		if((*houghlines)[i].ncluster<1)continue;
//		printf("houghline with e %f clusters %d (t,z) (%f,%f) (%f,%f) chi2/ndf %f  offz %f expt %f phi %f\n", (*houghlines)[i].sum_e, (*houghlines)[i].ncluster, (*houghlines)[i].start_t, (*houghlines)[i].start_z, (*houghlines)[i].end_t, (*houghlines)[i].end_z, (*houghlines)[i].chi2/(*houghlines)[i].ncluster,(*houghlines)[i].offset_z,(*houghlines)[i].GetExpectedT((*houghlines)[i].offset_z),(*houghlines)[i].phi);
		
		/*
			double ym = (*houghlines)[i].GetExpectedT(view_xmin);
			double yp = (*houghlines)[i].GetExpectedT(view_xmax);
			
						
			TLine *l=new TLine(view_xmin,ym,view_xmax,yp);
		*/	
		
			TLine *l=new TLine((*houghlines)[i].start_z,(*houghlines)[i].start_t,(*houghlines)[i].end_z,(*houghlines)[i].end_t);
		
			if((*houghlines)[i].primary)l->SetLineWidth(2);
			l->Draw();
						
		//	printf("fit line  --%f %f - %f %f -- theta r %f %f\n",view_xmin,ym,view_xmax,yp,houghlines[i].theta,houghlines[i].r);
						
	}
	

	








return;

///the old code....
/*

	 std::map<double, std::map<double, int>  >::iterator p_iterr;
     std::map<double, int>::iterator s_iterr;
	
printf("----view %d----\n\n",view);
	 Managed::ClusterManager *cluh = view==2?&mypoh->clusterer:&mypoh->clusterer;


		double off_z = cluh->GetClusterMap(view)->begin()->first;
		double off_t = cluh->GetClusterMap(view)->begin()->second.begin()->first;


/ *
	     for(p_iterr=cluh->cluster_map.begin();p_iterr!=cluh->cluster_map.end(); p_iterr++)
         {
         for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
         {
         
         	for(double i=0;i<2*3.141592;i+=3.141592*2/99)
         	{
       			h->Fill(i, (p_iterr->first-off_z) * cos(i) + (s_iterr->first-off_t) * sin(i), 1cluh->GetCluster(s_iterr->second)->e);  
         	}
         }
         }
         
  *//*
  
  Managed::ManagedCluster *largest=0;
  double last_plane=0;
  
//to take the largest in each plane  
/ *  
    for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
    {	
    	//new plane and a largest cluster in the previous plane?
		if(fabs(last_plane-p_iterr->first)>0.03 && largest)
		{
			for(double i=0;i<3.141592;i+=3.141592/(h->GetNbinsX()-1))
         	{
         		//double val = largest->e;
				double val=1;
				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i), val);  

//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)+2./(h->GetNbinsY())., val*.5); 
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)+4./(h->GetNbinsY())., val*.5*.5);
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)-2./(h->GetNbinsY())., val*.5); 
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)-4./(h->GetNbinsY())., val*.5*.5);      
			}
			largest=0;
			last_plane=p_iterr->first;
		}
	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluh->GetCluster(s_iterr->second);
			if(!largest)
			{
				largest=clus;
				continue;
			}
			if(largest->e<clus->e)largest=clus;
	  	}
	}
	if(fabs(last_plane-largest->z)>0.03 && largest)
	{
		for(double i=0;i<3.141592;i+=3.141592/(h->GetNbinsX()-1))
		{
			//double val = largest->e;
			double val=1;
			h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i),  val);  
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)+2./(h->GetNbinsY())., val*.5); 
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)+4./(h->GetNbinsY())., val*.5*.5);
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)-2./(h->GetNbinsY())., val*.5); 
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)-4./(h->GetNbinsY())., val*.5*.5);      
		
		}

		largest=0;
	}
*/       /*
         
         
     	TH2D copy;
	h->Copy(copy);    

//to take all within xxx% of the max hit in the event!

//find max value!
double maxvale=0;
    for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
    {	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluh->GetCluster(s_iterr->second);
			if(clus->e>maxvale)
			{
				maxvale=clus->e;
			}
		}
	}


    for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
    {	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluh->GetCluster(s_iterr->second);

//			if(clus->e>maxvale*0.1)
//			{
		for(double i=0;i<3.141592;i+=3.141592/(h->GetNbinsX()-1))
		{
			//double val = clus->e;
			double val=1;
			h->Fill(i, (clus->z-off_z) * cos(i) + (clus->t-off_t) * sin(i),  val);  
			copy.Fill(i, (clus->z-off_z) * cos(i) + (clus->t-off_t) * sin(i), clus->e);
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)+2./(h->GetNbinsY())., val*.5); 
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)+4./(h->GetNbinsY())., val*.5*.5);
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)-2./(h->GetNbinsY())., val*.5); 
//				h->Fill(i, (largest->z-off_z) * cos(i) + (largest->t-off_t) * sin(i)-4./(h->GetNbinsY())., val*.5*.5);      
			
}			
//			}
	  	}
	}
	
	
       	/ *
	h->SetStats(0);
	h->Draw("colz");
         
         return;
  */
   //peak filter the hough map!
   /*
   TH2D ch;
   h->Copy(ch);
	

h->Reset();
	
	int nx=ch.GetNbinsX();
	int ny=ch.GetNbinsY();   

	while(ch.GetMaximum()>1)
	{
		int x;
		int y;
		int z;
		int m = ch.GetMaximumBin();
		ch.GetBinXYZ(m,x,y,z);
		
/ *		printf("clearing hits around %d %d\n",x,y);
		ClearHit(&ch, ch.GetMaximum(), x+1,y);
		ClearHit(&ch, ch.GetMaximum(), x,y+1);
		ClearHit(&ch, ch.GetMaximum(), x,y-1);
		ClearHit(&ch, ch.GetMaximum(), x-1,y);
		ClearHit(&ch, ch.GetMaximum(), x+1,y+1);
		ClearHit(&ch, ch.GetMaximum(), x-1,y-1);
		ClearHit(&ch, ch.GetMaximum(), x+1,y-1);
		ClearHit(&ch, ch.GetMaximum(), x-1,y+1);
*/				/*		
		SaveHitGroup(&ch,(TH2D*)h, (double)ch.GetBinContent(x,y),(double)ch.GetBinContent(x,y), x,y);
		
		
    }    
         	
	h->SetStats(0);
	h->Draw("colz");


	if(view==2)padU->cd(1);
	if(view==3)padV->cd(1);
	
	
	/ *
	for(int i=1;i<=h->GetNbinsX();i++)
	for(int j=1;j<=h->GetNbinsY();j++)
	{
		if(h->GetBinContent(i,j)>=4)
		{
			
			double theta = h->GetXaxis()->GetBinLowEdge(i);
			double r = h->GetYaxis()->GetBinLowEdge(j);
			
			double ym = (- cos(theta)/sin(theta))*(viewU->GetXaxis()->GetXmin()-off_z)+r/sin(theta) + off_t;
			double yp = (- cos(theta)/sin(theta))*(viewU->GetXaxis()->GetXmax()-off_z)+r/sin(theta) + off_t;
						
			TLine *l=new TLine(viewU->GetXaxis()->GetXmin(),ym,viewU->GetXaxis()->GetXmax(),yp);
			l->Draw();
		
		}
	}
	*/
	/*
	multimap<double, std::pair<double,double> > top_map;
	
	int tops=0;	
/ *	for(int i=1;i<=h->GetNbinsX();i++)
	for(int j=1;j<=h->GetNbinsY();j++)
	{
		if(h->GetBinContent(i,j)>1)
		{
			top_map.insert(make_pair((double)h->GetBinContent(i,j), std::pair<int,int>(i,j) ));
			tops++;
			printf("topbin found x %d y %d\n",i,j);
		}
	}
*/

//h->Rebin2D(4,4);
/*	while(h->GetMaximum()>0)
	{
		double xv=0;
		double yv=0;
		double val=0;
		int cnt =0;
		
		int x;
		int y;
		int z;
		int m = h->GetMaximumBin();
		h->GetBinXYZ(m,x,y,z);
		
		GetPeakAreaAverage(xv, yv,val, cnt, x, y, (TH2D*)h);
		xv/=(double)cnt;
		yv/=(double)cnt;
		
		val=copy.GetBinContent(copy.FindBin(xv,yv));
		
		top_map.insert(make_pair((double)val, std::pair<double,double>(xv,yv)) );
		tops++;
	//	printf("topbin found x %f y %f\n",xv,yv);

	
	}
	
	
	
	printf("printing top %d\n",tops);
	
	multimap<double, std::pair<double,double> >::reverse_iterator it = top_map.rbegin();
	
	double view_xmin = view==2?viewU->GetXaxis()->GetXmin():viewV->GetXaxis()->GetXmin();
	double view_xmax = view==2?viewU->GetXaxis()->GetXmax():viewV->GetXaxis()->GetXmax();
	
	std::vector<HoughLine> houghlines;
	
	for(unsigned int k=0;k<tops;k++)
	{
		//int i=it->second.first;
		//int j=it->second.second;
		
			//double theta = h->GetXaxis()->GetBinCenter(i);
			//double r = h->GetYaxis()->GetBinCenter(j);
		double theta=it->second.first;
		double r = it->second.second;	
			
//			double ym = (- cos(theta)/sin(theta))*(view_xmin-off_z)+r/sin(theta) + off_t;
//			double yp = (- cos(theta)/sin(theta))*(view_xmax-off_z)+r/sin(theta) + off_t;


			HoughLine hl(theta, r, off_t, off_z);
	
			PrimaryShowerFinder::LoadCloseHits(&hl,&mypoh->clusterer,view,0.03);
			
			houghlines.push_back(hl);

		it++;
	}
	
	sort(&houghlines[0],&houghlines[houghlines.size()],CompareLength);
	
	//should remove duplicates --- (caused by the same clusters being matched to different hough map peaks)
	std::vector<HoughLine> houghlinestemp;
	for(int i=0;i<houghlines.size();i++)
	{
		HoughLine * l = &houghlines[i];
		HoughLine * r=0;
		if(i<houghlines.size()-1)r = &houghlines[i+1];
		if(houghlines[i].ncluster>1)houghlinestemp.push_back(houghlines[i]);		
		if(l && r && l->start_z == r->start_z && l->start_t == r->start_t && l ->end_z == r->end_z && l->end_t == r->end_t)
		{

			while (i<houghlines.size())
			{
				i++;
				HoughLine * l = &houghlines[i];
				HoughLine * r = &houghlines[i+1];
				if(! (l->start_z == r->start_z && l->start_t == r->start_t && l ->end_z == r->end_z && l->end_t == r->end_t)) break;
			}
		}
	
	}
	
	houghlines = houghlinestemp;
	
	sort(&houghlines[0],&houghlines[houghlines.size()],CompareForwardAndClusters);
	
	int max =-1;
	if(houghlines.size()>10)max = houghlines.size()-10;
	for(int i=houghlines.size()-1;i>max;i--)
	{
		if(houghlines[i].ncluster<1)break;
		printf("houghline with e %f clusters %d (t,z) (%f,%f) (%f,%f) chi2/ndf %f\n", houghlines[i].sum_e, houghlines[i].ncluster, houghlines[i].start_t, houghlines[i].start_z, houghlines[i].end_t, houghlines[i].end_z, houghlines[i].chi2/houghlines[i].ncluster);
		
			double ym = houghlines[i].GetExpectedT(view_xmin);
			double yp = houghlines[i].GetExpectedT(view_xmax);
			
						
			TLine *l=new TLine(view_xmin,ym,view_xmax,yp);
			
			if(i==houghlines.size()-1)l->SetLineWidth(2);
			l->Draw();
						
		//	printf("fit line  --%f %f - %f %f -- theta r %f %f\n",view_xmin,ym,view_xmax,yp,houghlines[i].theta,houghlines[i].r);
						
	}
	
	
*/
}


void HoughView::ClearHit(TH2D * his, double below_val, int curx, int cury)
{
	if(curx<1 || cury < 1 || curx> his->GetNbinsX() || cury >his->GetNbinsY())return;
	
//	printf("looking in %d %d\n",curx,cury);
	
	double thisval = his->GetBinContent(curx,cury);
	if(thisval>=below_val)return;
	
	
//	printf("val %f >= %f ?\n",thisval,below_val);
	ClearHit(his,below_val,curx,cury+1);
	ClearHit(his,below_val,curx+1,cury);
	ClearHit(his,below_val,curx,cury-1);
	ClearHit(his,below_val,curx-1,cury);
	ClearHit(his,below_val,curx+1,cury+1);
	ClearHit(his,below_val,curx-1,cury-1);
	ClearHit(his,below_val,curx+1,cury-1);
	ClearHit(his,below_val,curx-1,cury+1);

//	printf("recursion done\n");

	his->SetBinContent(curx,cury,0);	
	
}
 
void HoughView::SaveHitGroup(TH2D * his, TH2D * saver, double save_val, double with_val, int curx, int cury)
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

//assuming that the region is bounded by bins with 0 as their content
void HoughView::GetPeakAreaAverage(double &x, double &y,double &val, int & cnt, int curx, int cury, TH2D * hist)
{

	double thisval = hist->GetBinContent(curx,cury);

	if(thisval==0)return;
	
	val+=thisval;
	x+=hist->GetXaxis()->GetBinCenter(curx);
	y+=hist->GetYaxis()->GetBinCenter(cury);
	cnt++;
	hist->SetBinContent(curx,cury,0);

	
	GetPeakAreaAverage(x, y, val,cnt, curx+1, cury, hist);
	GetPeakAreaAverage(x, y, val,cnt, curx, cury+1, hist);
	GetPeakAreaAverage(x, y, val,cnt, curx, cury-1, hist);
	GetPeakAreaAverage(x, y, val,cnt, curx-1, cury, hist);
	GetPeakAreaAverage(x, y, val,cnt, curx+1, cury+1, hist);
	GetPeakAreaAverage(x, y, val,cnt, curx-1, cury-1, hist);
	GetPeakAreaAverage(x, y, val,cnt, curx+1, cury-1, hist);
	GetPeakAreaAverage(x, y, val,cnt, curx-1, cury+1, hist);

		
}
