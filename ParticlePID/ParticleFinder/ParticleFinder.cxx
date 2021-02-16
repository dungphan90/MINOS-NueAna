////////////////////////////////////////////////////////////////////////
/// $Id: ParticleFinder.cxx,v 1.11 2009/09/23 00:47:12 scavan Exp $
///
/// (Document me!)
///
/// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////

#include "TTree.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TMath.h"
#include "TROOT.h"

#include "NueAna/ParticlePID/ParticleFinder/ParticleFinder.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro

#include "CandNtupleSR/NtpSREvent.h"
#include "CandNtupleSR/NtpSRStrip.h"
#include "CandNtupleSR/NtpSRPulseHeight.h"

#include "MCNtuple/NtpMCStdHep.h"
#include "MCNtuple/NtpMCTruth.h"
#include "StandardNtuple/NtpStRecord.h"

#include "TruthHelperNtuple/NtpTHEvent.h"


#include "NueAna/ParticlePID/ParticleFinder/ParticleObject.h"
#include "NueAna/ParticlePID/ParticleFinder/Particle.h"

#include "NueAna/ParticlePID/ParticleFinder/ParticleEvent.h"

#include "Conventions/SimFlag.h"

#include "DataUtil/MCInfo.h"

#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedCluster.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ManagedHit.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterManager.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/HitManager.h"
#include "NueAna/ParticlePID/ParticleFinder/Managed/ClusterSaver.h"

#include "Calibrator/CalMIPCalibration.h"
#include "DatabaseInterface/DbiResultPtr.h"


#include "MuonRemoval/NtpMRRecord.h"
#include "MuonRemoval/NtpMREvent.h"

#include "NueAna/NueAnaTools/SntpHelpers.h"


#include "TClonesArray.h"

JOBMODULE(ParticleFinder, "ParticleFinder",
          "finds particles");


CVSID("$Id: ParticleFinder.cxx,v 1.11 2009/09/23 00:47:12 scavan Exp $");

TF1 * ParticleFinder::osceq=0;
SKZPWeightCalculator * ParticleFinder::skzpCalc=0;
//POT * ParticleFinder::pot=0; 

//......................................................................

ParticleFinder::ParticleFinder()
{
        lastrun=-1;

        recoCnt=0;
        pecutlevel=0;
        MIPperPE=60.;
        DoMRCC=0;

//      if(!pot)pot = new POT();
	if(!skzpCalc)skzpCalc=new SKZPWeightCalculator("DetXs",true);
  
  

        if(!osceq)osceq = new TF1("f2","sin(1.267*[0]*[1]/x)*sin(1.267*[0]*[1]/x)",0.,120.);

 


        
  ///
  /// (Document me!)
  ///
}
//......................................................................

ParticleFinder::~ParticleFinder()
{
  ///
  /// (Document me!)
  ///

  
  //if(pot)delete pot;
  //pot=0;
  if(skzpCalc)delete skzpCalc;
  skzpCalc=0;
  if(osceq)delete osceq;osceq=0;
  //kb =0;
}

//......................................................................

void ParticleFinder::BeginJob()
{
  ///
  /// (Document me!)
  ///
 
  
}

//......................................................................

void ParticleFinder::EndJob()
{
  ///
  /// (Document me!)
  ///

/*
  //MSG("ParticleFinder",Msg::kInfo)
  cout<<"Number of POT in this job: "<<pot->pot<<endl;
        
   TDirectory *savedir = gDirectory;
        
   TFile *fpf = dynamic_cast<TFile *>(gROOT->GetListOfFiles()->FindObject(kPOTTreeName.c_str()));
   if(fpf){
     fpf->cd();
     TTree *pottree = new TTree("pottree","pottree");
     pottree->Branch("ParticlePOT",&pot);
     pottree->Fill();
     pottree->Write();
     savedir->cd();
   }
   else{
     MSG("ParticleFinder",Msg::kError)<<"Could not find the file to write the pottree to, there will be no pottree"<<endl;
   }

*/

}

//......................................................................

JobCResult ParticleFinder::Reco(MomNavigator* mom)
{


   bool foundMR=false;
   bool foundST=false;

  NtpStRecord *str=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord","Primary"));
  NtpMRRecord *mr = 0;
  NtpStRecord *oldstr = 0;

  VldContext vc;  
   
  if(str){
  	foundST=true;	
	vc=str->GetHeader().GetVldContext();
     	ReleaseType::Release_t release = str->GetRelease();
     	//check for MR:
     	mr =  dynamic_cast<NtpMRRecord *>(mom->GetFragment("NtpMRRecord"));
     	if(mr){
       		if(ReleaseType::IsBirch(release)) {
	 	  		printf("please don't run with birch (particle finder)\n");exit(1);
	 	  		//oldstr=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord",
				//			    "NtpStRecordOld"));	
     		} else {
	    		oldstr=str;
	    		str=dynamic_cast<NtpStRecord *>(mom->GetFragment("NtpStRecord",
							  "MuonRemoved"));
     		}
     		if(oldstr) foundMR = true;
     	}
  }

  if(!foundST && !foundMR)
  {
      MSG("ParticleFinder",Msg::kWarning)<<"Got Nothing from mom"<<endl;
  
      return JobCResult::kFailed;

  }

  //std::vector<ParticleObjectHolder *>* hf;
  std::vector<TObject* > hft =( mom->GetFragmentList("ParticleObjectHolder"));

MSG("ParticleFinder",Msg::kDebug) << "Snarl " << str->GetHeader().GetSnarl() << endl;

  int make_new_holder=1;


  if(hft.size()>0)
  {
      make_new_holder=0;
  }
 

MSG("ParticleFinder",Msg::kDebug) << "size of found array " << hft.size()<<endl; 

  int ismc=str->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC;


  std::vector<const NtpSREvent* > evts = str->GetEvents();
  std::vector<const NtpSRStrip* > stp = str->GetStrips();

        std::vector<const NtpMCStdHep* > stdhep;
    std::vector<const NtpTHEvent* > thevt;
        std::vector<const NtpMCTruth* > mc;

        if(ismc)
        {
                stdhep = oldstr ? oldstr->GetMCStdHeps() : str->GetMCStdHeps();
                thevt = oldstr ? oldstr->GetTHEvents() : str->GetTHEvents();
                mc = oldstr ? oldstr->GetMCTruths() : str->GetMCTruths();
        }
        
        
  MSG("ParticleFinder",Msg::kDebug) << "Events " << evts.size() <<endl;
  if(evts.size()<1)
        {

/*              if(make_new_holder)
                {
                                RecCandHeader ntphdr(str->GetHeader().GetVldContext(),str->GetHeader().GetRun(),
                                        str->GetHeader().GetSubRun(),str->GetHeader().GetRunType(),str->GetHeader().GetErrorCode(),
                                        str->GetHeader().GetSnarl(),str->GetHeader().GetTrigSrc(),str->GetHeader().GetTimeFrame(),
                                        str->GetHeader().GetRemoteSpillType(),-1);
                
                
                        ParticleObjectHolder * h = new ParticleObjectHolder(ntphdr); 
                        
                
                        mom->AdoptFragment(h);
                }
*/
                return JobCResult::kPassed;  
        }

        ParticleBeamMon *bmon=0;
        MSG("ParticleFinder",Msg::kDebug) << "!!! "<<evts.size()<<" events, "<<  hft.size()<<" currently there"<<endl;

for(unsigned int ievt=0;ievt<evts.size();ievt++)
{


	int bestmatch=ievt;

	if(mr)
	{


	   bestmatch=-1;
	   //Int_t best_rmmu = i;
	   Float_t best_com = 0;
	   for(int xi=0;xi<mr->mrhdr.nmrevt;xi++){	  
	     NtpMREvent *ev = SntpHelpers::GetMREvent(xi,mr);
	     if(ev && ev->best_event==(int)ievt && ev->best_complete>best_com) {
	        best_com = ev->best_complete;
	      // best_rmmu = ev->orig_event;
		bestmatch = ev->orig_event;       
	//	       cout<<"found mr match "<<i<<" "<<ev->best_event<<" "<<xi<<endl;
	     }
	   }
	}



//int ievt=0;

        
//      if(ievt!=4 ||str->GetHeader().GetSnarl()!=2800034)continue;
//if(ievt!=0)continue;

        ParticleObjectHolder *n=0;


    //    if(make_new_holder || hft.size() < ievt+1)
     //   {
                        RecCandHeader ntphdr(str->GetHeader().GetVldContext(),str->GetHeader().GetRun(),
                                        str->GetHeader().GetSubRun(),str->GetHeader().GetRunType(),str->GetHeader().GetErrorCode(),
                                        str->GetHeader().GetSnarl(),str->GetHeader().GetTrigSrc(),str->GetHeader().GetTimeFrame(),
                                        str->GetHeader().GetRemoteSpillType(),ievt);
                n = new ParticleObjectHolder(ntphdr);

                
    //    }else{
     //           n=dynamic_cast<ParticleObjectHolder *>(hft[ievt]);
    //    }


        MSG("ParticleFinder",Msg::kDebug)<<"Generating event for run "<<str->GetHeader().GetRun()<<" subrun "<<str->GetHeader().GetSubRun()<<" snarl "<<str->GetHeader().GetSnarl()<<" event "<<ievt<<" matchingevt "<<bestmatch<<"\n";
        int mc_idx=-1;
        if(ismc && bestmatch>-1)mc_idx=thevt[bestmatch]->neumc;
                


if(ievt==0 || !bmon)
{

        bmon = new ParticleBeamMon(ntphdr);
                mom->AdoptFragment(bmon);
        

                bma.fname=kPOTTreeName;
                bma.ana(n,bmon);

                
                FillPot(bmon);
}

MSG("ParticleFinder",Msg::kDebug) << "partices in current evt "<< ievt << " array " << n->particles->GetEntries() << endl; 

        if(n->particles->GetEntries()>0)return JobCResult::kPassed;//continue; //dont rerun it


        finder.Reset();
        finder.SetMRCC(DoMRCC);
        finder.SetPOH(n);
        finder.SetMEUperGeV(MEUperGeV);







        if(ismc)
        {
        int beststdhep = -1;
	if(bestmatch>-1)beststdhep = thevt[bestmatch]->neustdhep;
        
        if(beststdhep>-1)//we may not have a match!
        {       
                double vx = stdhep[beststdhep]->vtx[0];
                double vy = stdhep[beststdhep]->vtx[1];
                double vz = stdhep[beststdhep]->vtx[2];
        
                finder.SetTrueVertex((vx+vy)/sqrt(2.0),(vy-vx)/sqrt(2.0),vz);

                n->mctrue.vtx_u=(vx+vy)/sqrt(2.0);
                n->mctrue.vtx_v=(-vx+vy)/sqrt(2.0);
                n->mctrue.vtx_z=vz;
        
                //want to save stdhep array...
                int thismc = stdhep[beststdhep]->mc;

                for(unsigned int i=beststdhep;i<stdhep.size() && stdhep[i]->mc==thismc;i++)
                {
                //      new(&(n->mctrue.stdhep[i-beststdhep])) NtpMCStdHep();
                //      *((NtpMCStdHep *)(n->mctrue.stdhep->AddrAt(i-beststdhep)))=*(stdhep[i]);
                        
                //      NtpMCStdHep *me = new NtpMCStdHep();
                //      n->mctrue.stdhep->AddAtAndExpand(me,i-beststdhep);
                //      *me = *(stdhep[i]);
        
                        NtpMCStdHep a=*(stdhep[i]);
        
                        a.parent[0]=a.parent[0]-beststdhep;
                        a.parent[0]=a.parent[1]-beststdhep;
                        a.child[0]=a.child[0]-beststdhep;
                        a.child[0]=a.child[1]-beststdhep;
                        a.index=i;
                        n->mctrue.stdhep.push_back(a);
        
        
                        if(mc_idx>-1)n->mctrue.visible_energy=mc[mc_idx]->tphu+mc[mc_idx]->tphv;
                }
        //      printf("saved %d stdheps\n",n->mctrue.stdhep.size());
        
        }
        }




        
        n->event.sr_vtx_u=evts[0]->vtx.u;
        n->event.sr_vtx_v=evts[0]->vtx.v;
        n->event.sr_vtx_z=evts[0]->vtx.z;       
        
        //set this to what we are actually filling based on sigcors
        n->event.visenergy=-1;//evts[ievt]->ph.mip;
                
        if(ismc && mc_idx>-1)
        {

                n->mctrue.iaction=mc[mc_idx]->iaction;
                n->mctrue.inu=mc[mc_idx]->inu;
                n->mctrue.iresonance=mc[mc_idx]->iresonance;
                n->mctrue.nuenergy=mc[mc_idx]->p4neu[3];

                n->mctrue.inunoosc=mc[mc_idx]->inunoosc;

                n->mctrue.flux=mc[mc_idx]->flux;
                
        /*      n->mctrue.osc_dm2=0.0024;
                n->mctrue.osc_L=735.; 
                n->mctrue.osc_sinth23=0.5;
                n->mctrue.osc_sin2th13=0.15;
        */      
                n->mctrue.type=(n->mctrue.iaction>0)*(abs(n->mctrue.inu)/2-5);
                //and for beam ves
                if(n->mctrue.inu==n->mctrue.inunoosc && TMath::Abs(n->mctrue.inu)==12)n->mctrue.type=4;
                
                
        }
        
        
        
        ////fluxweights/////
        
        if(ismc)
        {
        VldContext evt_vldc = n->GetHeader().GetVldContext();


        int det=evt_vldc.GetDetector(); 
        BeamType::BeamType_t beam=bmon->beamtype;
    NtpMCFluxInfo fi= n->mctrue.flux;
        double pt = sqrt(fi.tpx*fi.tpx+fi.tpy*fi.tpy);
        double pz = 1.*fi.tpz;
        int tptype = fi.tptype;
        int zbeam = BeamType::ToZarko(beam);
        n->mctrue.totbeamweight = skzpCalc->GetBeamWeight(det,zbeam,tptype,pt,pz,n->mctrue.nuenergy,n->mctrue.inu);
        }else{
                n->mctrue.totbeamweight = 1;
        }
        //////


        /////oscillation/////////
        if(ismc && n->GetHeader().GetVldContext().GetDetector()==Detector::kFar)
        {
        int otype=0;
        if(n->mctrue.iaction==1)
        {
                if(abs(n->mctrue.inu)==12)
                {
                        if(abs(n->mctrue.inunoosc)==12) otype=4;
                        else otype=2;
                }       
                if(abs(n->mctrue.inu)==14)
                {
                        if(abs(n->mctrue.inunoosc)==12) otype=6;
                        else otype=1;
                }       
                if(abs(n->mctrue.inu)==16)
                {
                        if(abs(n->mctrue.inunoosc)==12) otype=5;
                        if(abs(n->mctrue.inunoosc)==14) otype=3;
                }       
        }


        n->mctrue.osc_L=735;
        n->mctrue.osc_dm2=0.0024;
        n->mctrue.osc_sinth23=sin(3.141592/4.0);
        n->mctrue.osc_sin2th13=0.15;
        osceq->SetParameters(n->mctrue.osc_dm2,n->mctrue.osc_L);
        n->mctrue.oscprob=OscillationProb(osceq, otype, n->mctrue.nuenergy, n->mctrue.osc_sinth23* n->mctrue.osc_sinth23, n->mctrue.osc_sin2th13); 

//      printf("type %d oscprob %f e %f\n",n->mctrue.type,n->mctrue.oscprob, n->mctrue.nuenergy);

        }
        ///////
        


        Managed::HitManager hitmanager_u;
        Managed::ClusterManager clustermanager_u;
        clustermanager_u.SetHitManager(&hitmanager_u);

//      Managed::HitManager hitmanager_v;
//      Managed::ClusterManager clustermanager_v;
//      clustermanager_v.SetHitManager(&hitmanager_v);


        MSG("ParticleFinder",Msg::kDebug) << "-event " << ievt << "  with " << evts[ievt]->nstrip << " strips" <<endl; 

        double strip_e=0;


        //Setup and give context
        DbiResultPtr<CalMIPCalibration> fResPtr;
        fResPtr.NewQuery(str->GetHeader().GetVldContext(),10);

        //Get Scale
        const CalMIPCalibration* mipcal = fResPtr.GetRow(0);
        double sigcorpermip= mipcal->GetScale();

        if(recoCnt<1)printf("Using SigCorPerMip %f\n",sigcorpermip);

/*
        double sigcorpermip=0;
        if(n->GetHeader().GetVldContext().GetDetector()==Detector::kFar)sigcorpermip=kb->Get("FarSigCor/MIP");
        if(n->GetHeader().GetVldContext().GetDetector()==Detector::kNear)sigcorpermip=kb->Get("NearSigCor/MIP");
*/      
        //double pecutlevel=kb->Get("PECutLevel");
        //double MIPperPE=kb->Get("MIPperPE");


        //store the frontmost n hits in each view for containment

        double front_t[2][nfront];
        double front_z[2][nfront];

        for(int v=0;v<2;v++)
        for(int i=0;i<nfront;i++)
        {
                front_t[v][i]=0;
                front_z[v][i]=10000;
        }

        //Load processor with event strip data
        for(int i=0;i<evts[ievt]->nstrip;i++)
        {
        
                double strip_z = stp[evts[ievt]->stp[i]]->z;
                double strip_t = stp[evts[ievt]->stp[i]]->tpos;
                int view = stp[evts[ievt]->stp[i]]->planeview;
        
                if(! (view==2 || view ==3))continue;//we don't know what to do with any other view
                
                int forwardpos=nfront;
                for(;forwardpos>0;forwardpos--)
                        if(front_z[view-2][forwardpos-1]<strip_z)break;
                        
                //shift positions
                for(int k=nfront-1;k>=forwardpos;k--)
                {
                        front_z[view-2][k]=front_z[view-2][k-1];
                        front_t[view-2][k]=front_t[view-2][k-1];
                }
                
                if(forwardpos<nfront)
                {
                        front_z[view-2][forwardpos]=strip_z;
                        front_t[view-2][forwardpos]=strip_t;
                }
                
                        
        
                //pe cut????
                if((stp[evts[ievt]->stp[i]]->ph0.sigcor+stp[evts[ievt]->stp[i]]->ph1.sigcor)/sigcorpermip<pecutlevel*MIPperPE)continue;
        

                int plane = stp[evts[ievt]->stp[i]]->plane;
                int strip = stp[evts[ievt]->stp[i]]->strip;
                double energy = (stp[evts[ievt]->stp[i]]->ph0.sigcor+stp[evts[ievt]->stp[i]]->ph1.sigcor)/sigcorpermip;
                double tpos = stp[evts[ievt]->stp[i]]->tpos;
                double z = stp[evts[ievt]->stp[i]]->z;
                int planeview = stp[evts[ievt]->stp[i]]->planeview;
        
                finder.AddStrip(plane, strip, energy, tpos, z, planeview);
                strip_e+=energy;


                hitmanager_u.InsertHit(planeview, plane, strip, z, tpos, energy);
//              printf("adding strip %d %d %f\n",plane,strip,energy);
        }

//      printf("stripenergy %f\n",strip_e);
        /*
        printf("front strips\n");
        for(int v=0;v<2;v++)
        {
                printf("view %d\n",v+2);
        for(int i=0;i<nfront;i++)
        {
                if(front_z[v][i]<10000)printf("(z %f t %f)",front_z[v][i],front_t[v][i]);
        }
                printf("\n");
        }
        */
        
        CheckContainment(n, front_z, front_t);  
        //printf("forward contained: %d\n",n->event.contained);

        


                
                



        //clustermanager_u.DumpClusters();
//      clustermanager_v.MakeClusters(0.2,0.05,0.01);
        //clustermanager_v.DumpClusters();


        
        Managed::ClusterSaver * cs=&n->cluster_saver;

        
        clustermanager_u.SetClusterSaver(cs);
//      clustermanager_v.SetClusterSaver(cs);

        finder.SetClusterManagerU(clustermanager_u);
        finder.SetClusterManagerV(clustermanager_u);
        
n->event.visenergy=strip_e;

        MSG("ParticleFinder",Msg::kDebug) << finder.GetStrips() <<" strips added"<<endl; 

        //process event, returns particleobject 

        
                

        
        
        finder.Process(*n);


        
        //printf("p3d size %d\n",n->particles3d1->GetEntries());

        n->event.vtx_u=finder.vtx_u;
        n->event.vtx_v=finder.vtx_v;
        n->event.vtx_z=finder.vtx_z;    
        

        double x=(n->event.vtx_u-n->event.vtx_v)/sqrt(2.0);
        double y=(n->event.vtx_u+n->event.vtx_v)/sqrt(2.0);
        double z=n->event.vtx_z;
        
        //check for fiducial

        if(n->GetHeader().GetVldContext().GetDetector() ==Detector::kNear)
                if(IsInsideNearFiducial_Nue_Standard(x,y,z,ismc)==1) 
                        n->event.inFiducial=1;  
 
        if(n->GetHeader().GetVldContext().GetDetector()==Detector::kFar)        
                if(IsInsideFarFiducial_Nue_Standard(x,y,z,ismc)==1)
                        n->event.inFiducial=1;





                n->SetName(outObjectName.c_str());
                mom->AdoptFragment(n);

        recoCnt++;

}





/////////////////////
//a check to see whats there...
/*
 hft =( mom->GetFragmentList("ParticleObjectHolder"));
n=dynamic_cast<ParticleObjectHolder *>(hft[0]);

printf("#poh is %d   , #p3d %d\n",hft.size(),n->particles3d1->GetEntries());
*/



/////////////////////



  return JobCResult::kPassed; // kNoDecision, kFailed, etc.
}


//......................................................................





void ParticleFinder::CheckContainment(ParticleObjectHolder *p, double  front_z[2][nfront], double front_t[2][nfront])
{

        if(p->GetHeader().GetVldContext().GetDetector()!=Detector::kNear)       
        {
                p->event.contained=1;
                return;
        }       
        
        //we are in the near... make sure that the event starts in the fully instrumented region........

        //here we use a fidual volume like check....
        //actual dimensions taken from viewing strip.tpos vs strip.z is as folows:
        // u view:  -0.2<t<2.3  0<z<7
        // v view:  -2.3<t<0.2  0<z<7


        double SuperModule1Beg = 0.1;//1.01080;  //Data and MC values (according to DataUtil/infid.h on 10/02/07
        double SuperModule1End = 7;//4.99059;
        double radialInner = 0;
        double radialOuter = 1.;//0.8;
        
        double uCenter = 1.1513;
        double vCenter = -0.9537;
        
        
        int contained=1;
        for(int view=0;view<2;view++)
        for(int c=0;c<nfront;c++)
        {
        
                if(front_z[view][c]>=10000)continue;//not used.... its a short event
                
                
                if(front_z[view][c]<SuperModule1Beg || front_z[view][c]>SuperModule1End)
                {
                        contained=0;
                        break;
                }
                
                double tc=0;
                
                if(view==0)
                {
                        tc=front_t[view][c]-uCenter;
                }else{
                        tc=front_t[view][c]-vCenter;
                }
                
                tc=fabs(tc);//because its a radial cut and 0 is the center...
        
                if(tc < radialInner || tc > radialOuter)
                {
                        //printf("failing contained view %d st %f t %f z %f\n",view, tc,front_t[view][c], front_z[view][c]);
                        contained=0;
                        break;
                }
                
        }

        p->event.contained=contained;

////
/*
Double_t SuperModule1Beg = 1.01080;  //Data and MC values (according to DataUtil/infid.h on 10/02/07
  Double_t SuperModule1End = 4.99059;

    fX = 1./sqrt(2.0)*(fU - fV);
    fY = 1./sqrt(2.0)*(fU + fV);

        fu = 1./sqrt(2.)*(fx+fy);
        fv = 1./sqrt(2.)*(fy-fx);

  Double_t radialInner = 0;
  Double_t radialOuter = 0.8;
  Double_t xCenter = 1.4885;
  Double_t yCenter = 0.1397;

        uCenter = 1.1513;
        vCenter = -0.9537;


*/
//////
        
                

}







const Registry& ParticleFinder::DefaultConfig() const
{
  ///
  /// Supply the default configuration for the module
  ///
  static Registry r; // Default configuration for module

  // Set name of config
  std::string name = this->GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  //r.Set("POTTreeFileName","pottree.root");
  

  // Set values in configuration
  r.UnLockValues();

  r.Set("DoMRCC",0);
  r.Set("OutObjectName","Normal");
  r.Set("PECutLevel",2);
  r.Set("MIPperPE",0.1);
  r.Set("MEUperGeV",25.);
  r.LockValues();


  




  return r;
}

//......................................................................

void ParticleFinder::Config(const Registry& r)
{
  ///
  /// Configure the module given the Registry r
  ///
  
      const char* tmps;
        int tmpi;
        double tmpd;
  //  if(r.Get("POTTreeFileName",tmps)){kPOTTreeName=tmps;}
        if(r.Get("OutObjectName",tmps))
        {
                outObjectName=tmps;
        //      cout << "Setting Particle Finder output object name to "<<outObjectName<<"\n";
        }
        if(r.Get("MEUperGeV",tmpd)){MEUperGeV=tmpd;}
        if(r.Get("DoMRCC",tmpi)){DoMRCC=tmpi;}
        if(r.Get("PECutLevel",tmpd)){pecutlevel=tmpd;}
        if(r.Get("MIPperPE",tmpd)){MIPperPE=tmpd;}


}

//......................................................................

void ParticleFinder::Reset()
{
  ///
  /// (Document me!)
  ///
}




void ParticleFinder::FillPot(ParticleBeamMon *bmon)
{
                int ismc=bmon->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC;
                VldContext evt_vldc = bmon->GetHeader().GetVldContext();


                if(!ismc)
                {
                        //pot->pot+=bmon->GetPot();
                }else{
                
                        std::string relName = bmon->GetTitle();

                        std::string mcinfo = "";
                        if(bmon->GetHeader().GetVldContext().GetSimFlag() == SimFlag::kMC){
                        mcinfo = "Daikon";
        //              string temp = mchdr.geninfo.codename;
        //              if(temp.size() != 0){   mcinfo = temp;  }
                        }
                                        
                        //ReleaseType::Release_t rel = ReleaseType::MakeReleaseType(relName, mcinfo);
                
                        
                        //near for every snarl
/*                      if(evt_vldc.GetDetector()==Detector::kNear){
              double thispot = MCInfo::GetMCPoT(Detector::kNear, bmon->beamtype, rel);
              pot->pot+=thispot;
           // printf("near pot %f bt %d r %d\n",thispot,bmon->beamtype,rel);
            }
*/                      
                        //far only once per run
/*                      if(evt_vldc.GetDetector()==Detector::kFar && 
                bmon->GetHeader().GetRun()!=lastrun){
                        double thispot=6.5e8;//
                        //dogwood0 isn't showing up as daikon....
                        //MCInfo::GetMCPoT(Detector::kFar, bmon->beamtype, rel);
                        pot->pot+=thispot;
                        //printf("far pot %f bt %d r %d\n",thispot,bmon->beamtype,rel);
                }
  */  
                }
                
                
                
/*              pot->nsnarls++;
                if(bmon->GetHeader().GetRun()!=lastrun)
                {
                        pot->nruns++;
                        pot->beamtype = bmon->beamtype;
                }
*/              lastrun=bmon->GetHeader().GetRun();
                

}


////////////////////////////////////////////////////////////////////////





//Oscillation code from Tingjun  7/15/08
double ParticleFinder::OscillationProb(TF1* f2, int ntype, double NuE, double sinth23, double sin2th13) {

 // sinth23 : sine square theta 23
 // sin2th13 : sine square 2 theta 13
 //
 double OscProb = 0 ;
 double NumuToNutau ;
 double NumuToNue ;
 double NueSurvival ;
 double NumuSurvival ;
 double NueToNutau  ;
 double NueToNumu;

 NumuToNutau = 4.*sinth23*(1.-sinth23)*pow(1-sin2th13/4,2) ;
 NumuToNutau *= f2->Eval(TMath::Abs(NuE)) ;

 NumuToNue = sinth23*sin2th13*f2->Eval(TMath::Abs(NuE)) ;
 
//printf("ting p = %f  sinth23 %f sin2th13 %f * l/e %f\n",NumuToNue,sinth23,sin2th13,f2->Eval(TMath::Abs(NuE))) ;

//printf("ting l/e term %f\n",f2->Eval(TMath::Abs(NuE))) ;

 NueSurvival = 1.- sin2th13*f2->Eval(TMath::Abs(NuE)) ;
 //NueSurvival = 1.;

 NumuSurvival = 1. - NumuToNutau - NumuToNue ;
 //NumuSurvival = 1.;

 NueToNutau = (1.-sinth23)*sin2th13*f2->Eval(TMath::Abs(NuE)) ;

 NueToNumu = NumuToNue;

 if (ntype==0) OscProb = 1 ;
 else if (ntype==1) OscProb = NumuSurvival ;
 else if (ntype==2) OscProb = NumuToNue ;
 else if (ntype==3) OscProb = NumuToNutau ;
 else if (ntype==4) OscProb = NueSurvival ;
 else if (ntype==5) OscProb = NueToNutau ;
 else if (ntype==6) OscProb = NueToNumu;
//  cout<<"Period = "<<f2->Eval(TMath::Abs(NuE))<<endl ;
//  cout<<"Oscillation probability = "<<OscProb<<endl ;

 return OscProb ;
} 




int ParticleFinder::IsInsideNearFiducial_Nue_Standard(double x, double y, double z, bool /*isMC*/)
{
  Double_t SuperModule1Beg = 1.01080;  //Data and MC values (according to DataUtil/infid.h on 10/02/07
  Double_t SuperModule1End = 4.99059;

  Double_t radialInner = 0;
  Double_t radialOuter = 0.8;
  Double_t xCenter = 1.4885;
  Double_t yCenter = 0.1397;

  Bool_t zContained = false;
  Bool_t xyContained = false;

  Double_t r = TMath::Sqrt((x-xCenter)*(x-xCenter) + (y-yCenter)*(y-yCenter));

  if( z >= SuperModule1Beg && z <=SuperModule1End)
     zContained = true;
  if( r >= radialInner && r <= radialOuter)
     xyContained = true;

  Int_t retVal = 0;
  if(zContained && xyContained) retVal = 1;
  if(!zContained) retVal = -1;
  if(!xyContained) retVal -= 2;

  return retVal;
}

int ParticleFinder::IsInsideFarFiducial_Nue_Standard(double x, double y, double z, bool isMC){

//cout << "vtx "<< x <<" "<<y<<" "<<z<<" "<<isMC<<endl;

  Double_t SuperModule1Beg =  0.49080;   // These are data values
  Double_t SuperModule2Beg = 16.27110;
  Double_t SuperModule1End = 14.29300;
  Double_t SuperModule2End = 27.98270;

  if(isMC){
    SuperModule1Beg =  0.47692;   // These are mc values
    SuperModule2Beg = 16.26470;
    SuperModule1End = 14.27860;
    SuperModule2End = 27.97240;
  }

  Double_t radialInner = 0.50;
  Double_t radialOuter = TMath::Sqrt(14.0);
  Bool_t zContained = false;
  Bool_t xyContained = false;

  Double_t r = TMath::Sqrt(x*x + y*y);

  if( (z >= SuperModule1Beg && z <=SuperModule1End) ||
      (z >= SuperModule2Beg && z <=SuperModule2End) )
     zContained = true;

  if( r >= radialInner && r <= radialOuter)
     xyContained = true;

  Int_t retVal = 0;
  if(zContained && xyContained) retVal = 1;
  if(!zContained) retVal = -1;
  if(!xyContained) retVal -= 2;

  return retVal;  //  1 contained, -1 out of bounds z
                  //  -2 oob xy, -3 oob both
}
