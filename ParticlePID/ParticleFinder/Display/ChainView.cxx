#include "NueAna/ParticlePID/ParticleFinder/Display/ChainView.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterManager.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterSaver.h"
#include "TPolyLine.h"
#include "TMarker.h"
#include "TBox.h"

#include "TArrow.h"
#include "TMath.h"

#include "Conventions/SimFlag.h"
#include <vector>
#include "MCNtuple/NtpMCStdHep.h"
#include "MCNtuple/NtpMCRecord.h"
#include "MCNtuple/NtpMCTruth.h"

#include "TruthHelperNtuple/NtpTHEvent.h"

using namespace std;

ChainView::ChainView():Page()
{
	viewU=0;
	viewV=0;
	padU=0;
	padV=0;


	info6=new TLatex();
	info7=new TLatex();
	info8=new TLatex();
	info9=new TLatex();
	info10=new TLatex();
	info11=new TLatex();
	info12=new TLatex();
	info13=new TLatex();
}


ChainView::~ChainView()
{}


void ChainView::BuildDisplay(TCanvas *c)
{
	myCanvas=c;
	padV=new TPad("chainview_padv","chainview_padu",0,0,0.7,0.5);
	padU=new TPad("chainview_padu","chainview_padv",0,0.5,0.7,1.0);
	padMC=new TPad("chainview_padmc","chainview_padmc",0.71,0.0,1.0,0.8);
	padInfo=new TPad("chainview_padinfo","chainview_padinfo",0.71,0.8,1.0,1.0);
	myCanvas->cd();
	padU->Draw();
	padV->Draw();
	padMC->Draw();
	padInfo->Draw();
	
}


void ChainView::DrawEvent(ParticleObjectHolder * poh,  NtpStRecord *ntp,int ntpEvt)
{
	padU->Clear();
	padV->Clear();
	padMC->Clear();
	padInfo->Clear();
	
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

	viewU = new TH2D("hu","Clusters U", nz, minz,maxz,nu,minu, maxu);
	viewV = new TH2D("hv","Clusters V", nz, minz,maxz,nv,minv, maxv);
	viewU->SetStats(0);
	viewV->SetStats(0);
	viewU->SetDirectory(0);
	viewV->SetDirectory(0);

	if(mypoh)
	{	
		padU->cd();	
		viewU->Draw();
		DrawClusters(2,viewU);
		DrawVertex(2);
		if(isMC)DrawTrueVertex(2);
		DrawChains(2);

		

		padV->cd();
		viewV->Draw();
		DrawClusters(3,viewV);
		DrawVertex(3);
		if(isMC)DrawTrueVertex(3);
		DrawChains(3);	

	}
	
	if(ntp)
	{
		padMC->cd();
		if(isMC)DrawInteractionDiagram(ntp,ntpEvt);
		padInfo->cd();
		if(isMC)DrawBasicInfo(ntp,ntpEvt);
	}
	
	double maxz_u = viewU->GetMaximum();
	double maxz_v = viewV->GetMaximum();
	double maxz_both = maxz_u>maxz_v?maxz_u:maxz_v;
	viewU->GetZaxis()->SetRangeUser(0, maxz_both);
	viewV->GetZaxis()->SetRangeUser(0, maxz_both);
	
	padU->Update();
	padV->Update();
	padMC->Update();
	padInfo->Update();
	myCanvas->Update();
	
}


void ChainView::DrawChains(int view)
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
                        
           // printf("parent %d children ", ch->finished[i].myId);
                       
            for(unsigned int j=0;j<c.size();j++)
            {               
                                
            	Chain *todraw=ch->GetChain(c[j]);

                //printf("%d ",todraw->myId);

                
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
						if(!ch->GetChain(todraw->parentChain))
                		{
                			cout<<"ERROR!!!!!! bad chain in hitview.cxx line 272\n";
                			return;
                		
                		}
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
 
                }//printf("\n");
        }}      
                
}       





void ChainView::DrawClusters(int view, TH2 * h)
{

	 std::map<double, std::map<double, int>  >::iterator p_iterr;
     std::map<double, int>::iterator s_iterr;
	

	Managed::ClusterManager *cluh = view==2?&mypoh->clusterer:&mypoh->clusterer;

	//maxe should hold the minimum e cutoff to restrict one hit per plane
	double maxe=0;
	double lastplane=0;

	double maxpe=0;
    double secmaxpe=0;
	for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
    {

		if(TMath::Abs(p_iterr->first-lastplane)>0.02){maxe = maxe<secmaxpe?secmaxpe:maxe;  //printf("maxe %f\n",maxe);
		}

         if(TMath::Abs(p_iterr->first-lastplane)>0.02){
         	maxpe=0;
         	secmaxpe=0;
         	//printf("newplane\n");
         }
   
        //printf("plane %f\n",p_iterr->first);
         for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
         {
         
         	double ce = cluh->GetCluster(s_iterr->second)->e;
         	if(ce>secmaxpe)
         	{
         		if(ce>maxpe)
         		{
         			secmaxpe=maxpe+0.001; //need to be slightly larger than second largest peak to make sure we cut it out
         			maxpe=ce;
         		}else secmaxpe=ce+0.001;
         		
         		
         	}
		
         	//printf("t %f  e %f\n",s_iterr->first,ce);
         	//printf("updating maxe %f secmaxpe %f\n",maxpe,secmaxpe);
			
         }

         lastplane=p_iterr->first;
   }
   maxe = maxe<secmaxpe?secmaxpe:maxe;
	//printf("maxe %f\n",maxe);

	     for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
         {
         for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
         {
         
         	if(!cluh->GetCluster(s_iterr->second))continue;
       		
       		h->Fill( p_iterr->first,s_iterr->first, cluh->GetCluster(s_iterr->second)->e);  


       		
         }
         }
         
         
         	
	Managed::ClusterSaver *	cluhs = &mypoh->cluster_saver;
    for(p_iterr=cluhs->GetClusterMap(view)->begin();p_iterr!=cluhs->GetClusterMap(view)->end(); p_iterr++)
    {	

        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			//Managed::ManagedCluster * clus = cluhs->GetCluster(s_iterr->second);
			h->Fill( p_iterr->first,s_iterr->first, cluhs->GetCluster(s_iterr->second)->e); 

	  	}
	}
         
         
         printf("CLUSTER SAVER HAS %d HITS\n",(int)cluhs->clusters.size());
         
		h->SetStats(0);
		h->Draw("colz");

	     for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
         {
         for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
         {
         
       	/*	if( cluh->GetCluster(s_iterr->second)->e>maxe)
       		{
       		    TMarker *m=new TMarker(p_iterr->first,s_iterr->first, 25);
       			m->Draw("same");
       		
       		}
*/

			if( cluh->GetCluster(s_iterr->second)->GetStatus()!=1)
       		{

       		    TMarker *m=new TMarker(p_iterr->first,s_iterr->first, 25);
       			m->Draw("same");
       		
       		}
       		
       		//draw the size of the cluster....
       		

       		
       		
         }
         }

/*
    for(p_iterr=cluhs->GetClusterMap(view)->begin();p_iterr!=cluhs->GetClusterMap(view)->end(); p_iterr++)
    {	

        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{

       			Managed::ManagedCluster * c = cluh->GetCluster(s_iterr->second);
       			TBox *b = new TBox(c->zmin,c->tmin,c->zmax,c->tmax);
       			b->Draw("same");
		}
	}
	

    for(p_iterr=cluh->GetClusterMap(view)->begin();p_iterr!=cluh->GetClusterMap(view)->end(); p_iterr++)
    {	

        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{

       			Managed::ManagedCluster * c = cluh->GetCluster(s_iterr->second);
       			
       			TBox *b = new TBox(c->zmin-0.01,c->tmin,c->zmax+0.01,c->tmax);
       			
       			char txt[200];
       			sprintf(txt,"(%f,%f) (%f %f) e %f\n",c->zmin,c->tmin,c->zmax,c->tmax,c->e);
       			b->SetToolTipText(txt);
       			b->Draw("1same");
		}
	}	
	*/
/*
    for(p_iterr=cluhs->GetClusterMap(view)->begin();p_iterr!=cluhs->GetClusterMap(view)->end(); p_iterr++)
    {	
        for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
    	{
			Managed::ManagedCluster * clus = cluhs->GetCluster(s_iterr->second);


		TMarker *m=new TMarker(clus->z,clus->t, 27);
		m->Draw("same");
			

	  	}
	}
*/	
	


	//show cluster saver hits
/*
	Managed::ClusterSaver *cs = &mypoh->cluster_saver;

	for(p_iterr=cs->cluster_map.begin();p_iterr!=cs->cluster_map.end(); p_iterr++)
    {

         for(s_iterr=p_iterr->second.begin();s_iterr!=p_iterr->second.end(); s_iterr++)
         {
         
         	if(cs->GetCluster(s_iterr->second)->view!=view)continue;
         	double ce = cs->GetCluster(s_iterr->second)->e;
         	h->Fill( p_iterr->first,s_iterr->first, ce);  
			TMarker *m=new TMarker(p_iterr->first,s_iterr->first, 27);
       		m->Draw("same");
         }

   	}
*/



}


void ChainView::DrawVertex(int view)
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

void ChainView::DrawTrueVertex(int view)
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
  
  
  



void ChainView::DrawInteractionDiagram(NtpStRecord *str,int index)
{
	cout <<"int "<<&str<<" "<<index<<endl;

	if(!str || index<0)return;
	
	cout <<"interaction\n";

  //modified by steve cavanaugh aug 9, 2006 to take into account double particles in mc with IstHEP=2 
  //with modifications for readability
  

/*
	TPad *tp1=new TPad("tp1","tp1",0.0,0.0,1,0.75);
		
	tp1->Draw();



	tp1->cd();
*/	
	
	 
    std::vector<const NtpMCTruth* > thmc = str->GetMCTruths();
    if(thmc.size()<1)return;
     std::vector<const NtpTHEvent* > thevt = str->GetTHEvents();    
    index =  thevt[index]->neumc;

	


  
  vector<const NtpMCStdHep *> hep;
 // if (foundST) hep = SntpHelpers::GetStdHepArray(index, st);
 // else hep = SntpHelpers::GetStdHepArray(index, mc);

	hep = str->GetMCStdHeps(index);




  Int_t nStdHep = int(hep.size());
  Int_t *indicesToUse = new Int_t[nStdHep];
  Int_t *parent = new Int_t[nStdHep];
  Int_t *parent1 = new Int_t[nStdHep];
  Int_t incomingNeutrino = -1;
  Int_t cnt = 0;



  for(int i=0;i<nStdHep;i++){    
    if(hep[i]->mc==index) {
      indicesToUse[cnt] = i;



      //parent[i] = hep[i]->parent[0];

      //in the case where we have more than 1 event per snarl we might 
      //have indices which refer to the actual index, not the one for the particular event.... 
      //so make sure we match this up properly
      for(int j=0;j<nStdHep;j++){
        if ((Int_t)hep[j]->index==hep[i]->parent[0]){
          parent[i]=j;
          break;
        }
      }

      
      for(int j=0;j<nStdHep;j++){
        if ((Int_t)hep[j]->index==hep[i]->parent[1]){
          parent1[i]=j;
          break;
        }
      }
	



      if(hep[i]->IstHEP==0){
	if(abs(hep[i]->IdHEP)==12||abs(hep[i]->IdHEP)==14
	   ||abs(hep[i]->IdHEP)==16) {
	  incomingNeutrino=i;


	}
	parent[i] = -1; //don't draw arrows to initial state particles
	parent1[i] = -1;
      }
      cnt++;
    }
    else{
      parent[i] = -1;
      parent1[i] = -1;
    }
  }


  




  //make arrows and markers
  TArrow *arrow[1000];
  TMarker *marker[1000];
  for(int i=0;i<nStdHep;i++) {
    arrow[i] = new TArrow(0,0,0,0,.01,"|>"); 
    arrow[i]->SetLineWidth(1);
    arrow[i]->SetLineStyle(3);
    arrow[i]->SetFillColor(39);
    arrow[i]->SetLineColor(39);
    marker[i] = new TMarker(0,0,24);
    marker[i]->SetMarkerSize(.2);
    marker[i]->SetMarkerColor(38);
  }

  //now loop through valid stdhep entries and fill variables  
  Double_t Available[5] = {0.9,0.7,0.7,0.1};  
  
  //double max_a_3=0.7;
  
  //cout<<"nStdHep: "<<nStdHep<<"\n";

  for(int i=0;i<cnt;i++){
    
    int toUse = indicesToUse[i];
    if(hep[toUse]->IstHEP==999){
      parent[i]=-1;
      continue;
    }



    //    cout<<toUse<<"  "<<hep[toUse]->IstHEP<<"   "<<nStdHep<<endl;
    //cout<<parent[toUse]<<"   "<<parent1[toUse]<<endl;



    
    
    //sc - avoid any double taus' charms' etc, but set the proper parent/child 
    //(remove the particle between true parent and child)

    if((hep[toUse]->child[0]==hep[toUse]->child[1]) &&
       (hep[toUse]->IstHEP==2 || hep[toUse]->IstHEP==3 || hep[toUse]->IstHEP==14) && 
       hep[toUse]->child[0]!=-1){
      int j;
      int child = hep[toUse]->child[0];
      for(j=0;j<cnt;j++){
        if ((Int_t)hep[j]->index==child)break;
      }

      //we've lost arrow tail info, so fix it
      if (j<nStdHep){
	arrow[j]->SetX1(arrow[i]->GetX1());
	arrow[j]->SetY1(arrow[i]->GetY1());
      
	parent[i]=-1;
      
	continue;
      }else{
	cout<<"Error! - linking failed, index out of range! \n";
      }
    }

    //sc - avoid any "later decays" not comming from IstHEP=1 decays .. .these are children of double taus, etc

    if((hep[toUse]->IstHEP==205||hep[toUse]->IstHEP==1205)&&
       (hep[parent[toUse]]->IstHEP!=1 || hep[parent1[toUse]]->IstHEP!=1)){
      parent[i]=-1;
      //cout <<"breaking out of loop for particle "<<i<<"\n";
      continue;  
    }



    
    Double_t mom = sqrt(hep[toUse]->p4[0]*hep[toUse]->p4[0] + 
		       hep[toUse]->p4[1]*hep[toUse]->p4[1] + 
		       hep[toUse]->p4[2]*hep[toUse]->p4[2]);
    double x=0.,y=0.;
    int col=0;
    char text[256];

    //set x,y    
    if(hep[toUse]->IstHEP==0) {
      x = 0.05;
      y=Available[0]; Available[0] -= 0.05;
    }
    else if(hep[toUse]->IstHEP==2) {
      x = 0.15;
      y=Available[1]; Available[1] -= 0.05;   
    }
    else if(hep[toUse]->IstHEP==11) {
      x = 0.05;
      y=Available[0]; Available[0] -= 0.05; 
    }
    else if(hep[toUse]->IstHEP==3||hep[toUse]->IstHEP==14) {   //sc - allow for both IstHEP=3 and 14 intermediate decays
      x = 0.3;
      y=Available[1]; Available[1] -= 0.05;
    }
    else if(hep[toUse]->IstHEP==1){
      x = 0.55;
      y=Available[2]; Available[2] -= 0.05;
    }
    else if((hep[toUse]->IstHEP==205||hep[toUse]->IstHEP==1205)){
      x = 0.8;
      y=Available[3]; Available[3] += 0.05; //seems to be opposite ordered....
    }
 
    //set colour and label (and override y in special cases)
    if(abs(hep[toUse]->IdHEP)==12) { //nue
      if(parent[toUse]==incomingNeutrino)  y = 0.9; 
      sprintf(text,"#nu_{e}"); col = 3;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#bar{#nu}_{e}");
    }
    else if(abs(hep[toUse]->IdHEP)==14) { //numu
     if(parent[toUse]==incomingNeutrino)  y = 0.9;  
     sprintf(text,"#nu_{#mu}"); col = 4;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#bar{#nu}_{#mu}");
    }
    else if(abs(hep[toUse]->IdHEP)==16) { //nutau
     if(parent[toUse]==incomingNeutrino)  y = 0.9;  
     sprintf(text,"#nu_{#tau}"); col = 105;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#bar{#nu}_{#tau}"); 
    }    
    else if(abs(hep[toUse]->IdHEP)==11) { //e
      if(parent[toUse]==incomingNeutrino) y = 0.9;          
      sprintf(text,"e^{-}"); col = 3;
      if(hep[toUse]->IdHEP<0) sprintf(text,"e^{+}");
    }
    else if(abs(hep[toUse]->IdHEP)==13) { //mu
      if(parent[toUse]==incomingNeutrino) y = 0.9;               
      sprintf(text,"#mu^{-}"); col = 4;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#mu^{+}");
    }
    else if(abs(hep[toUse]->IdHEP)==15) { //tau
      //tau will decay so it is easier to read if it is not on the top line...
      // if(parent[toUse]==incomingNeutrino) y = 0.9;                  
      sprintf(text,"#tau^{-}"); col = 105;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#tau^{+}"); 
    }
    else if(hep[toUse]->IdHEP==22) { //photon
      sprintf(text,"#gamma"); col = 9;      
    }
    else if(hep[toUse]->IdHEP>1000000000) { //nucleus
      y = 0.8;
      sprintf(text,"nucleus(%i,%i)",int((hep[toUse]->IdHEP-1e9)/1e6),
	      int((hep[toUse]->IdHEP-1e9 - 1e6*int((hep[toUse]->IdHEP-1e9)
						  /1e6))/1e3)); 
      col = 15;
    }
    else if(hep[toUse]->IdHEP==2112){ 
      sprintf(text,"neutron"); col = 28;
    }
    else if(hep[toUse]->IdHEP==2212){
      sprintf(text,"proton"); col = 2;
    }
    else if(abs(hep[toUse]->IdHEP)==211) {
      sprintf(text,"#pi^{+}"); col = 6;
      if(hep[toUse]->IdHEP<0) sprintf(text,"#pi^{-}");
    }
    else if(hep[toUse]->IdHEP==111) {
      sprintf(text,"#pi^{0}"); col = 7;
    }
    else if(hep[toUse]->IdHEP==130) {
      sprintf(text,"K^{0}_{L}"); col = 31;
    }
    else if(hep[toUse]->IdHEP==310) {
      sprintf(text,"K^{0}_{S}"); col = 31;
    }
    else if(hep[toUse]->IdHEP==311) {
      sprintf(text,"K^{0}"); col = 31;
    }
    else if(abs(hep[toUse]->IdHEP)==321) {
      sprintf(text,"K^{+}"); col = 31;
      if(hep[toUse]->IdHEP<0) sprintf(text,"K^{-}"); col = 31;
    }
    else if(hep[toUse]->IdHEP==28) {
      sprintf(text,"Geantino"); col = 46;
      if(hep[toUse]->IdHEP<0) sprintf(text,"K^{-}"); col = 31;
    }
    else {
      sprintf(text,"ID: %i",hep[toUse]->IdHEP); col=43;
    }

    sprintf(text,"%s [%.1f GeV/c]",text,mom);
    
    arrow[toUse]->SetX2(x-0.02);   
    arrow[toUse]->SetY2(y-0.02);   
    marker[toUse]->SetX(x-0.02);
    marker[toUse]->SetY(y-0.02);

    for(int j=0;j<nStdHep;j++){
   
      if(parent[j]==toUse){
	arrow[j]->SetX1(x-0.02);
	arrow[j]->SetY1(y-0.02);
	//cout<<"writing arrows from "<<j << " to "<<parent[j]<<"\n";
      }
    }
 
    TLatex *tex = new TLatex(x,y,text);
    char texname[256];
    sprintf(texname,"tex%i",i);
    tex->SetName(texname);
    tex->SetTextSize(0.03);
    tex->SetTextColor(col);
    tex->Draw();
  }

  for(int i=0;i<nStdHep;i++){
    if(parent[i]==-1){
      delete arrow[i];
      delete marker[i];
    }
    else {
    
    //	if(arrow[i]->GetX2()== 0.78 )
   // 	{
   // 		arrow[i]->SetY2(arrow[i]->GetY2()+0.5);//max_a_3-Available[3]);
    //		marker[i]->SetY(marker[i]->GetY()+0.5);//max_a_3-Available[3]);
    //	}
    
    
      arrow[i]->Draw();
      marker[i]->Draw();
      
    }
  }
  




  Double_t minAvail = 0;
  for(int i=0;i<4;i++){
    if(Available[i]<minAvail) minAvail = Available[i];
  }

  delete [] indicesToUse;
  delete [] parent;



}




void ChainView::DrawBasicInfo(NtpStRecord *str,int index)
{


  info6->Clear();
  info7->Clear();
  info8->Clear();
  info9->Clear();
  info10->Clear();
  info11->Clear();
  info12->Clear();
  info13->Clear();

  //char text1[100];
  //char text2[100];
  //char text3[100];
  //char text31[100];
  //char text4[100];
  char text5[100];
  char text6[100];
  char text7[100];
  char text8[100];
  char text9[100];
  char text10[100];
  //char text11[100];
  //char text12[100];

  	if(!str || index <0)return;  
    
    sprintf(text5,"Run: %d, Snarl: %d, Event: %d",str->GetHeader().GetRun(),str->GetHeader().GetSnarl(),index);
   


  std::vector<const NtpMCStdHep* > stdhep = str->GetMCStdHeps();


double energy =0;

if(stdhep.size()>0)
{
	printf("STDHEPS %d %f %f %f %f\n",(int)stdhep.size(),stdhep[0]->p4[0],stdhep[0]->p4[1],stdhep[0]->p4[2],stdhep[0]->p4[3]);
	energy = stdhep[0]->p4[3];
}else return;


        int mcinfoc=1;

        double textsize=0.12;

if(stdhep.size()==1)
{

 sprintf(text6,"Single Particle Energy %f",energy);



    info6->SetText(0.05,1-1*(textsize*1.1),text5);
    info6->SetTextSize(textsize);
    info6->SetTextColor(mcinfoc);
	info6->Draw();

    info7->SetText(0.05,1-2*(textsize*1.1),text6);
    info7->SetTextSize(textsize);
    info7->SetTextColor(mcinfoc);
        info7->Draw();

return; 
 }
    
    std::vector<const NtpMCTruth* > thmc = str->GetMCTruths();
     std::vector<const NtpTHEvent* > thevt = str->GetTHEvents();
    
    
    
    
    if(thmc.size()<1)return;
    
    const NtpMCTruth * mctruth = thmc[thevt[index]->neumc];
    

 	char tmptxt[100];
 	tmptxt[0]=0;//sprintf(tmptxt,"");
    if (mctruth->iaction==0){
      sprintf(text6,"NC event %s",tmptxt);
    }
    else{
      if (abs(mctruth->inunoosc)==12&&abs(mctruth->inu)==12){
	sprintf(text6,"beam #nu_{e} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==12&&abs(mctruth->inu)==14){
	sprintf(text6,"beam #nu_{e}#rightarrow#nu_{#mu} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==12&&abs(mctruth->inu)==16){
	sprintf(text6,"beam #nu_{e}#rightarrow#nu_{#tau} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==14&&abs(mctruth->inu)==14){
	sprintf(text6,"#nu_{#mu}#rightarrow#nu_{#mu} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==14&&abs(mctruth->inu)==12){
	sprintf(text6,"#nu_{#mu}#rightarrow#nu_{e} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==14&&abs(mctruth->inu)==16){
	sprintf(text6,"#nu_{#mu}#rightarrow#nu_{#tau} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==16&&abs(mctruth->inu)==14){
	sprintf(text6,"#nu_{#tau}#rightarrow#nu_{#mu} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==16&&abs(mctruth->inu)==12){
	sprintf(text6,"#nu_{#tau}#rightarrow#nu_{e} CC %s",tmptxt);
      }
      if (abs(mctruth->inunoosc)==16&&abs(mctruth->inu)==16){
	sprintf(text6,"#nu_{#tau}#rightarrow#nu_{#tau} CC %s",tmptxt);
      }
    }
    if (mctruth->iresonance==1001 && mctruth->iaction==1){
      sprintf(text7,"quasi-elastic");
    } 
    if (mctruth->iresonance==1001 && mctruth->iaction==0){
      sprintf(text7,"elastic");
    } 
    if (mctruth->iresonance==1002){
      sprintf(text7,"resonance production");
    }
    if (mctruth->iresonance==1003){
      sprintf(text7,"DIS");
    }
    if (mctruth->iresonance==1004){
      sprintf(text7,"coherent production");
    }
    sprintf(text8,"E_{#nu}: %.1fGeV, E_{shw}: %.1fGeV",TMath::Abs(mctruth->p4neu[3]),TMath::Abs(mctruth->p4shw[3]));
    sprintf(text9,"Y: %.2f, emfrac: %.2f",mctruth->y,mctruth->emfrac);
    sprintf(text10,"E_{#mu}: %.1fGeV,  E_{el}: %.1fGeV",TMath::Abs(mctruth->p4mu1[3]),TMath::Abs(mctruth->p4el1[3]));

/*
    if (foundST) {
      sprintf(text11,"E_{#pi^{0}}_tot: %.2fGeV, E_{#pi^{0}}_neugen: %.2fGeV",stdhepinfo.epi0_total,stdhepinfo.epi0_neugen);
      sprintf(text12,"E_{#pi^{0}}_abs: %.2fGeV, E_{#pi^{0}}_intranuke: %.2fGeV",stdhepinfo.epi0_abs,stdhepinfo.epi0_intranuke);
    }
*/

//	int mcinfoc=1;
	
//	double textsize=0.12;

    info6->SetText(0.05,1-1*(textsize*1.1),text5);
    info6->SetTextSize(textsize);
    info6->SetTextColor(mcinfoc);
    

    info7->SetText(0.05,1-2*(textsize*1.1),text6);
    info7->SetTextSize(textsize);
    info7->SetTextColor(mcinfoc);
    
    info8->SetText(0.05,1-3*(textsize*1.1),text7);
    info8->SetTextSize(textsize);
    info8->SetTextColor(mcinfoc);

    info9->SetText(0.05,1-4*(textsize*1.1),text8);
    info9->SetTextSize(textsize);
    info9->SetTextColor(mcinfoc);
    
    info10->SetText(0.05,1-5*(textsize*1.1),text9);
    info10->SetTextSize(textsize);
    info10->SetTextColor(mcinfoc);

    info11->SetText(0.05,1-6*(textsize*1.1),text10);
    info11->SetTextSize(textsize);
    info11->SetTextColor(mcinfoc);

/*
    info12->SetText(0.05,0.14,text11);
    info12->SetTextSize(textsize);
    info12->SetTextColor(mcinfoc);

    info13->SetText(0.05,0.08,text12);
    info13->SetTextSize(textsize);
    info13->SetTextColor(mcinfoc);
  */  
    
  //  padr->cd(1);

	info6->Draw();
    info7->Draw();
    info8->Draw();
    info9->Draw();
    info10->Draw();
    info11->Draw();
   // info12->Draw();
   // info13->Draw();
    
}


 
