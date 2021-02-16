////////////////////////////////////////////////////////////////////////
//
// FILL_IN: [Document your code!!]
//
// scavan@fas.harvard.edu
////////////////////////////////////////////////////////////////////////
#include <cmath>
#include <fstream>
#include <map>
#include <sstream>
#include <string>

#include "TH1F.h"
#include "TFile.h"
#include "Conventions/Detector.h"
#include "NueAna/NueRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro

#include "AnalysisNtuples/ANtpDefaultValue.h"                                   
#include "TClass.h"
#include "TDataType.h"
#include "TDataMember.h"
#include "TRealData.h"
#include "PETrimmer.h"
#include "StandardNtuple/NtpStRecord.h"

#include "CandNtupleSR/NtpSRRecord.h"

#include "TClonesArray.h"

#include "CandNtupleSR/NtpSRStrip.h"  
#include "CandNtupleSR/NtpSREvent.h" 
#include "CandNtupleSR/NtpSRShower.h" 
#include "CandNtupleSR/NtpSRTrack.h" 

#include "TObject.h"

#include "VertexFinder/NtpVtxFinder/NtpVtxFinder.h"

JOBMODULE(PETrimmer, "PETrimmer",
          "Filter out strips below a given PE");
          
CVSID("$Id: PETrimmer.cxx,v 1.15 2009/08/12 15:22:46 scavan Exp $");

//......................................................................

PETrimmer::PETrimmer():
  pecut(0.0),
  updateEventEnergy(1)
{}

//......................................................................

PETrimmer::~PETrimmer()
{}

//......................................................................
void PETrimmer::BeginJob()
{

}


JobCResult PETrimmer::Reco(MomNavigator* mom)
{
   bool foundST=false;
                                                                               
//   VldContext vc;
   //first ask mom for a NtpStRecord

//iterate over all NtpStRecords in MOM

std::vector<TObject*> records = mom->GetFragmentList("NtpStRecord");


for(unsigned int i=0;i<records.size();i++)
{

   NtpStRecord *str=dynamic_cast<NtpStRecord *>(records[i]); //(mom->GetFragment("NtpStRecord","Primary"));
   if(str){
     foundST=true;
  //   vc=str->GetHeader().GetVldContext();
   }else{
   	continue;
   }
   //cout<<"Running over "<<str->GetName()<<endl;

   TrimRecord(str);

}

   return JobCResult::kPassed; // kNoDecision, kFailed, etc.

}


void PETrimmer::TrimRecord(NtpStRecord *str)
{
   DumpIt(str);



	MSG("PETrimmer",Msg::kDebug)<< "Cutting out strips below "<<pecut<<" pe"<<endl;
	
	std::vector<const NtpSRStrip* > stp = str->GetStrips();
	
	if(stp.size()<1)return; //no strips!

	std::vector< NtpSREvent* > evts;
	std::vector< NtpSRShower* > shw;
	std::vector< NtpSRTrack* > trk;


MSG("PETrimmer",Msg::kDebug)<< "Snarl " << str->GetHeader().GetSnarl() << endl;



	TClonesArray* tc = str->evt;	
	for(int i=0;i<tc->GetEntries();i++)
		evts.push_back((NtpSREvent*)tc->At(i));
		
	TClonesArray* ts = str->shw;	
	for(int i=0;i<ts->GetEntries();i++)
		shw.push_back((NtpSRShower*)ts->At(i));
		
	TClonesArray* tt = str->trk;	
	for(int i=0;i<tt->GetEntries();i++)
		trk.push_back((NtpSRTrack*)tt->At(i));
		

   for(int i=0;i<(int)evts.size();i++)
        {       
	  MSG("PETrimmer",Msg::kDebug)<< "Event "<<i<<" strips "<<evts[i]->nstrip<<" idx "<<evts[i]->index<<endl;
	};




///////last shower strip bad index check.....
        for(int i=0;i<(int)shw.size();i++)
          {
       


                //found a missing strip in the mrcc far data...  hopefully this will fix it!
                if(shw[i]->stp[shw[i]->nstrip-1]<0)
                {
                        int nstrip = shw[i]->nstrip;
                        MSG("PETrimmer",Msg::kError)<<"Shower "<<i<<" is missing the last strip!!!! adjusting strip size from "<< nstrip <<" to "<<(nstrip-1)<<endl;
                        shw[i]->nstrip=shw[i]->nstrip-1;
			shw[i]->nstpcnt=shw[i]->nstrip;
                }

          }

	

































      

		

		
//	float *pepass = new float[stp.size()];	
//	float *stripe_trk = new float[3 * stp.size()];
//	float *stripe_shw = new float[3 * stp.size()];


        int pepass[100000];  //will store if that strip passes the pe cut
        float stripe_trk[3 * 100000]; //energy of strip as appearing in track objects - will hold sigcor, mip, gev
        float stripe_shw[3 * 100000]; //energy of strip as appearing in shower object
	

	//these are arrays that will store 1 if that index in shower or event is to be kept
	//hopfully there wont be a snarl with more than 1000 showers or events
	int goodshower[1000]; //if its bad (a view with no hits).. we will remove them from the event!
	for(int i=0;i<1000;i++)goodshower[i]=1;

	int goodevent[1000];
        for(int i=0;i<1000;i++)goodevent[i]=1;

        int goodtrack[1000];
        for(int i=0;i<1000;i++)goodtrack[i]=1;


	for(int i=0;i<(int)stp.size();i++)
	{
		pepass[i]=stp[i]->ph0.pe+stp[i]->ph1.pe - pecut > 0.0001 ? 1: 0;

		for(int j=0;j<3;j++)
		{
			stripe_trk[j + 3*i]=0;
			stripe_shw[j + 3*i]=0;
		}
	}
		
	



	for(int i=0;i<(int)shw.size();i++)
	{
		MSG("PETrimmer",Msg::kDebug)<< "Shw "<<i<<" nstrip " << shw[i]->nstrip<< " stpcnt "<<shw[i]->nstpcnt<<endl;
		NtpSRShower *sh = shw[i];
		goodshower[i]=CleanShower(sh,stp,pepass,stripe_shw);
	}	

	
	
	for(int i=0;i<(int)trk.size();i++)
	{
		MSG("PETrimmer",Msg::kDebug)<< "Trk "<<i<<" nstrip " << trk[i]->nstrip<< " stpcnt "<<trk[i]->nstpcnt<<endl;
		NtpSRTrack *tr = trk[i];
		goodtrack[i]=CleanTrack(tr,stp,pepass,stripe_trk);
	}	

	
	for(int i=0;i<(int)evts.size();i++)
	{
		MSG("PETrimmer",Msg::kDebug)<< "Event "<<i<<" nstrip " << evts[i]->nstrip<< " stpcnt "<<evts[i]->nstpcnt<<endl;
		NtpSREvent* evt = evts[i];
		goodevent[i]=CleanEvent(evt,shw,trk,stp,pepass,stripe_trk,stripe_shw,goodshower,goodtrack);
//                 MSG("PETrimmer",Msg::kInfo)<<"Status: " << goodevent[i]<<endl;


		if(!goodevent[i])continue;

		NtpVtxFinder finder(i,str);
		if(finder.FindVertex() >0)
		{
			evt->vtx.x=finder.VtxX();
			evt->vtx.y=finder.VtxY();
			evt->vtx.z=finder.VtxZ();
			evt->vtx.u=finder.VtxU();
			evt->vtx.v=finder.VtxV();
			evt->vtx.t=finder.VtxT();
			evt->vtx.plane=finder.VtxPlane();
		}

                int thisVtxPlane=evt->vtx.plane;
		int invtxwindow=0;
 		for(int j=0;j<evt->nstpcnt;j++)
		{
			if(evt->stp[j]<0)continue;
			if(pepass[evt->stp[j]])
            if(stp[evt->stp[j]]->plane >= thisVtxPlane && stp[evt->stp[j]]->plane < (thisVtxPlane+5))
            {
            	invtxwindow=1;
                break;
            }   
        }              

        if(!invtxwindow)
        {
			float evtpes=0;
            for(int j=0;j<evt->nstpcnt;j++)
            {
            	if(evt->stp[j]<0)continue;
                if(pepass[evt->stp[j]])
                	evtpes+= (stp[evt->stp[j]]->ph0.pe + stp[evt->stp[j]]->ph1.pe);
            }
           	MSG("PETrimmer",Msg::kWarning)<<"removing "<<evtpes<<" pe event that has no strips in vtx window!"<<endl;

        	goodevent[i]=0;
        }

	}


      
     
	//remove bad events from the ntpstrecord..

	std::vector<NtpSREvent *>todelrec;

	//int tot=evts.size();
	int off=0;
        for(int i=0;i<(int)evts.size();i++)
        {
        	if(!goodevent[i+off])
                {
		
			MSG("PETrimmer",Msg::kDebug)<< "Event removed from snarl "<<i+off<<endl;
	
			todelrec.push_back(evts[i]);
			evts.erase(evts.begin()+i);
			i--;	
			off++;
			//str->evthdr.nevent--; //don't adjust this... we are not resizing the event array due to truth matching
		}

        }

	
	//clean up primary tracks, showers... 
	//for now, if we deleted a primary, set it to -1 (don't try to find a replacement!)

	for(int i=0;i<(int)evts.size();i++)
	{
           if(evts[i]->primtrk > -1){
              bool isOK = false;
              int longest = 0;
              int longestTrackIndex = -1;

              for(int j=0;j<(int)evts[i]->ntrack && !isOK; j++)
              {
                 int index = evts[i]->trk[j];
                 if(index == -1) continue;
                 if(index == evts[i]->primtrk) isOK = true;

                 NtpSRTrack* track = trk[index];
                 int length =  TMath::Abs(track->plane.end - track->plane.beg);
                 if(length > longest){
                    longest = length;  longestTrackIndex = j;
                 }
              }
              if(!isOK) {  evts[i]->primtrk = longestTrackIndex; 
  //               cout<<"Assigning new primary track"<<longest<<endl;
              }
           }

	   //check prim shw
	   if(evts[i]->primshw>-1) {
	      bool isOK = false;

              //code ruthlessly borrowed from CandEventHandle::SetPrimaryShower

              double biggestE = -1;
              double closeE = -1;
              int primShwIndex = -1;
              int largestShwIndex = -1;
              int closeShwIndex = -1;
              double dzLarge = 0;
              double dzClose = 1e6;
              NtpSRTrack *primaryTrack = 0;
              if(evts[i]->primtrk != -1) primaryTrack = trk[evts[i]->trk[evts[i]->primtrk]];

	      for(int j=0;j<(int)evts[i]->nshower && !isOK;j++)
	      {
                 int index = evts[i]->shw[j];
                 if(index == -1) continue;
		 if(index==evts[i]->primshw) isOK = true;
                 
                 NtpSRShower* shower = shw[index];
                 if(shower->ph.mip > biggestE){
                    largestShwIndex = j;
                    biggestE = shower->ph.mip;
                    if(primaryTrack)  dzLarge = fabs(shower->vtx.z - primaryTrack->vtx.z); 
                 }
                 if(primaryTrack){
                     if(fabs(shower->vtx.z - primaryTrack->vtx.z) < dzClose){
			 dzClose = fabs(shower->vtx.z - primaryTrack->vtx.z);
                         closeShwIndex = j;
                         closeE = shower->ph.mip;
                      }
                 }
              }
              if(!isOK){  //the primary shower has died, we must elect a new leader
                float eratio = 1;
                if(biggestE > 0) eratio = closeE/biggestE;

//                cout<<primShwIndex<<"  "<<eratio<<"  "<<dzLarge<<"  "<<dzClose<<endl;
                if(closeShwIndex > -1 && (dzLarge - dzClose) > 1 && eratio > 0.2)
                   primShwIndex = closeShwIndex;
                else
                   primShwIndex = largestShwIndex;

                if(primShwIndex > -1 && primaryTrack){
                  NtpSRShower *PShw = shw[evts[i]->shw[primShwIndex]];
                  if(PShw == 0) cout<<"Failed to load shower - VERY BAD"<<endl;
                  if(fabs(PShw->vtx.z - primaryTrack->vtx.z) > 0.5 &&
                          PShw->shwph.linNCgev < 2.0){
	 	    primShwIndex = -1;
                  }
                }
                evts[i]->primshw = primShwIndex;
//                cout<<"New primary shower assigned!"<<primShwIndex<<"  "<<largestShwIndex
//                    <<"  "<<closeShwIndex<<"  "<<endl;
              }
//              cout<<"Event Primary shw: "<<evts[i]->primshw<<"  "<<evts[i]->nshower<<endl;                  
	   } 
        }

	//I am not directly delete from the tclonesarray, as it seems to cause problems....
	//instead, clear the tclones array, and put the pointers for the good events back in
//check that events make sense....
//str->evt->Clear();


for(int i=0;i<(int)todelrec.size();i++)
	str->evt->Remove(todelrec[i]);  //does not shift entires...


////!!!!!! do not compress... this will mess up the truth matching.....
//.......   str->evt->Compress();




/*
   for(int i=0;i<(int)str->evt->GetEntries();i++)
     {
       ( (NtpSREvent*)str->evt->At(i))->index=i;
     }
*/
   //check that we saved them
   evts.clear();
   for(int i=0;i<tc->GetEntries();i++)
 	{
		NtpSREvent * e = ((NtpSREvent*)str->evt->At(i));
		if(!e)continue;
		e->index=i;
	  	evts.push_back((NtpSREvent*)str->evt->At(i));
	}


   for(int i=0;i<(int)evts.size();i++)
   {      
	  MSG("PETrimmer",Msg::kDebug)<< "Event "<<i<<" strips "<<evts[i]->nstrip<<" idx "<<evts[i]->index<<endl;
   };


      

    
   
  
 









//for debuggins...


	for(int i=0;i<(int)evts.size();i++)
        {
        std::ostringstream s;
		s<<" showers ";
		for(int j=0;j<evts[i]->nshower;j++)
			s<<" " << evts[i]->shw[j] <<" ";
                s<<" tracks ";
		for(int j=0;j<evts[i]->ntrack;j++)
                        s<<" " << evts[i]->trk[j] <<" ";


        MSG("PETrimmer",Msg::kDebug)<< "Event done    strips "<<evts[i]->nstrip<< s.str() << endl;

	}
	
//	delete[] pepass;
//	delete[] stripe_trk;
//	delete[] stripe_shw;


//	DumpIt(str);

   

}


        bool showergreater(std::pair<float, NtpSRShower *>  p1, std::pair<float, NtpSRShower *>  p2){return p1.first < p2.first;}
        bool trackgreater(std::pair<float, NtpSRTrack *>  p1, std::pair<float, NtpSRTrack *>  p2){return p1.first < p2.first;}
 


int PETrimmer::CleanEvent(NtpSREvent * evt, std::vector< NtpSRShower* > shw, std::vector< NtpSRTrack* > trk, std::vector<const NtpSRStrip* > stp, int *pepass, float *stripe_trk, float* stripe_shw,int*goodshower,int*goodtrack)
{
     int goodevent=1;

		//remove bad showers from this event....

		std::map<int,int> shwstpmap;

		//mark bad showers as bad, record which strips are effected
		for(int i=0;i<evt->nshower;i++)
		{
			if(!goodshower[evt->shw[i]])
			{
                            //cout<<"Removing shw: "<<i<<endl;	
				for(int j=0;j<shw[evt->shw[i]]->nstrip;j++)
                                shwstpmap[shw[evt->shw[i]]->stp[j]]=1; 					


				evt->shw[i]=-1;
				
			}	
		}
		
		//remove bad showers from event
		for(int i=0;i<evt->nshower;i++)
		{
			if(evt->shw[i]<0)
			{
				for(int j=i+1;j<evt->nshower;j++)
					evt->shw[j-1]=evt->shw[j];
				evt->shw[evt->nshower-1]=-1;
				evt->nshower--;
				i--;
			}
		}




                //mark bad trackss as bad, record which strips are effected
                for(int i=0;i<evt->ntrack;i++)
                {      
			if(evt->trk[i]<0)continue; 
                        if(!goodtrack[evt->trk[i]])
                        {       
                           //cout<<"removing track"<<i<<endl;
                                for(int j=0;j<trk[evt->trk[i]]->nstrip;j++)
                                shwstpmap[trk[evt->trk[i]]->stp[j]]=1;

                         
                                evt->trk[i]=-1;
                         
                        }
                }
                
                //remove bad tracks from event
                for(int i=0;i<evt->ntrack;i++)
                {       
                        if(evt->trk[i]<0)
                        {       
                                for(int j=i+1;j<evt->ntrack;j++)
                                        evt->trk[j-1]=evt->trk[j];
                        	evt->trk[evt->ntrack-1]=-1; 
			       evt->ntrack--;
                                i--;
                        }
                }









		//which strips are in the tracks? removing these strips due to removed showers has no effect on event energy
		for(int i=0;i<evt->ntrack;i++)
		{
			for(int j=0;j<trk[evt->trk[i]]->nstrip;j++)
			{
				shwstpmap[trk[evt->trk[i]]->stp[j]]=0;
			
		
			}
		}

		//if we have strips shared between good and bad showers, we are not removing those strips from the event
                for(int i=0;i<evt->nshower;i++)
                {
                       for(int j=0;j<shw[evt->shw[i]]->nstrip;j++)
			{
                                shwstpmap[shw[evt->shw[i]]->stp[j]]=0;
	

			}
		}
		


		std::map<int,int>::iterator iter;
		for(iter=shwstpmap.begin();iter!=shwstpmap.end();iter++)
		{
			
			
			if(iter->second==1)
			{
			  int j=evt->nstrip;
			  for(int i=0;i<evt->nstrip;i++) //find the start index in evt->stp for this bad strip
			    {
			      if(evt->stp[i]==iter->first)
				{
				  j=i;
				  break;
				}
			    }
				for(int i=j;i<evt->nstrip-1;i++)
				{
				  evt->stp[i]=evt->stp[i+1]; //remove that strip from the event strip list
				  
				}
				//mark the last strip, so we dont get confused
				evt->stp[evt->nstrip-1]=-1;
				evt->nstrip--; //adjust both .... we are loosing access to the end of the array
				evt->nstpcnt--;
			}	

		}
                //cout<<"Remaining strips: "<<evt->nstrip<<endl;
		if(evt->nstrip<1)return 0; //mark this event for deletion	


		int cutcount =0;
	
		//temporary vectors for saving good hits - better than erasing single entries in arrays at a time!
		std::vector<float> stpph0sigmap;
		std::vector<float> stpph0mip;
		std::vector<float> stpph0gev;
		std::vector<float> stpph1sigmap;
		std::vector<float> stpph1mip;
		std::vector<float> stpph1gev;
		std::vector<int> stpid;

		if(updateEventEnergy){
		//we are recomputing the energies... clear out those variables...
		evt->ph.raw =0;
		evt->ph.pe =0;
		evt->ph.sigcor =0;
		evt->ph.siglin =0;
		evt->ph.sigmap =0;
		//evt->ph.mip =0; //don't change this
		evt->ph.gev =0;		
		}

		int us=0;
		int vs=0;
		for(int j=0;j<evt->nstpcnt;j++)
		{

	
			if(evt->stp[j]<0)continue;



			
			if(pepass[evt->stp[j]])
			{
			

                                 //printf("stp %d plane %d\n",stp[evt->stp[j]],stp[evt->stp[j]]->plane);


				//keep it!
				stpid.push_back(evt->stp[j]);
/*				stpph0sigmap.push_back(evt->stpph0sigmap[j]);
				stpph1sigmap.push_back(evt->stpph1sigmap[j]);
				stpph0gev.push_back(evt->stpph0gev[j]);
				stpph1gev.push_back(evt->stpph1gev[j]);
				stpph0mip.push_back(evt->stpph0mip[j]);
				stpph1mip.push_back(evt->stpph1mip[j]);
*/




				if(updateEventEnergy){

			//recompute the energy!!!!!!!
				evt->ph.raw += (stp[evt->stp[j]]->ph0.raw + stp[evt->stp[j]]->ph1.raw);
				evt->ph.pe += (stp[evt->stp[j]]->ph0.pe + stp[evt->stp[j]]->ph1.pe);
				evt->ph.sigcor += (stp[evt->stp[j]]->ph0.sigcor + stp[evt->stp[j]]->ph1.sigcor);
				evt->ph.siglin += (stp[evt->stp[j]]->ph0.siglin + stp[evt->stp[j]]->ph1.siglin);

				//give preference to track calibrated strips???
				//this puts preference for all tracks before all showers
				//should it primary track, primary shower, other stuff, etc...???
						
				evt->ph.sigmap += stripe_trk[3*evt->stp[j]]>0 ? stripe_trk[3*evt->stp[j]] :stripe_shw[3*evt->stp[j]];
			//	evt->ph.mip += stripe_trk[1+3*evt->stp[j]]>0 ? stripe_trk[1+3*evt->stp[j]] :stripe_shw[1+3*evt->stp[j]];
//				evt->ph.gev += stripe_trk[1+3*evt->stp[j]]>0 ? stripe_trk[1+3*evt->stp[j]] :stripe_shw[1+3*evt->stp[j]];

//trying a new method of energy recalculation..... store it in gev - that is recalculated anyways!
//give priority to showers, in order of energy, and then to tracks, in order of energy
//this is to get energies for nue like events!
				std::vector<std::pair<float, NtpSRShower *> >showers;
                                std::vector<std::pair<float, NtpSRTrack *> >tracks;

				for(int k=0;k<evt->nshower;k++)
					showers.push_back(make_pair((float)shw[evt->shw[k]]->ph.mip,shw[evt->shw[k]]));
				for(int k=0;k<evt->ntrack;k++)	
					tracks.push_back(make_pair((float)trk[evt->trk[k]]->ph.mip,trk[evt->trk[k]]));


				std::sort(showers.begin(),showers.end(),showergreater);
				std::sort(tracks.begin(),tracks.end(),trackgreater);		

				float ees[100000];
				for(int k=0;k<100000;k++)ees[k]=0;

				for(unsigned int k=0;k<showers.size();k++)
				for(int j=0;j<showers[k].second->nstrip;j++)
				{
					int sidx = showers[k].second->stp[j];
					if(sidx>0 && sidx < 100000)				
					{
						if(ees[sidx]==0)ees[sidx]=showers[k].second->stpph0mip[j]+showers[k].second->stpph1mip[j];
					}
				}

                                for(unsigned int k=0;k<tracks.size();k++)
                                for(int j=0;j<tracks[k].second->nstrip;j++)
                                {
                                        int sidx = tracks[k].second->stp[j];
                                        if(sidx>0 && sidx < 100000)
                                        {
                                                if(ees[sidx]==0)ees[sidx]=tracks[k].second->stpph0mip[j]+tracks[k].second->stpph1mip[j];
                                        }
                                }

				float summeu=0;
				for(int k=0;k<evt->nstrip;k++)
					if(ees[evt->stp[k]] > 0 && ees[evt->stp[k]] < 100000)
						summeu+=ees[evt->stp[k]];

				evt->ph.gev=summeu;


				}

			         if(stp[evt->stp[j]]->planeview==2)
                                        us++;
                                if(stp[evt->stp[j]]->planeview==3)
                                        vs++;


			}else{
				cutcount++;
			}
		}


		MSG("PETrimmer",Msg::kDebug)<<"Cut "<<cutcount<<" strips from event!"<<endl;
                MSG("PETrimmer",Msg::kDebug)<<"Shower us "<<us << " vs "<< vs <<endl;

		if(us==0 || vs==0)
		{
			goodevent=0;
			return goodevent;
		}
		








		if (cutcount==0)return goodevent;  //we didn't cut anything, so continue	

		//rewrite event stp vectors
		for(int j=0;j<evt->nstpcnt;j++)
		{
			evt->stp[j]=-1;
	/*		evt->stpph0sigmap[j]=0;
			evt->stpph1sigmap[j]=0;
			evt->stpph0gev[j]=0;
			evt->stpph1gev[j]=0;
			evt->stpph0mip[j]=0;
			evt->stpph1mip[j]=0;*/
		}


		//deleteing these arrays and putting new ones in works on my mac, but not on minos machines.... why?

//		delete[] evt->stp;
/*		delete[] evt->stpph0sigmap;
		delete[] evt->stpph1sigmap;
		delete[] evt->stpph0gev;
		delete[] evt->stpph1gev;
		delete[] evt->stpph0mip;
		delete[] evt->stpph1mip;*/
		
//		evt->stp = new Int_t[stpid.size()];
/*		evt->stpph0sigmap = new Float_t[stpid.size()];
		evt->stpph1sigmap = new Float_t[stpid.size()];
		evt->stpph0gev = new Float_t[stpid.size()];
		evt->stpph1gev = new Float_t[stpid.size()];
		evt->stpph0mip = new Float_t[stpid.size()];
		evt->stpph1mip = new Float_t[stpid.size()];	//these arrays dont seem to be in the ntpst file!
		
*/		
		int planes[1000];
		for(int i=0;i<1000;i++)planes[i]=0;	
		
		evt->nstpcnt=stpid.size();
		evt->nstrip=stpid.size();  ///incorrect... we did not change array size.. workaround for other mods...
		for(int j=0;j<(int)stpid.size();j++)
		{
			if(stpid[j]<0 || stpid[j]>(int)stp.size()-1)continue;
			
			int p=stp[stpid[j]]->plane;
			if(p>-1 && p < 1000) planes[p]=stp[stpid[j]]->planeview; //2 or 3 for u or v

//		        evt->stp[j]=-1;
                        if(stpid[j]>-1 && stpid[j] < (int)stp.size())   
				evt->stp[j]=stpid[j];
/*			evt->stpph0sigmap[j]=stpph0sigmap[j];
			evt->stpph1sigmap[j]=stpph1sigmap[j];
			evt->stpph0gev[j]=stpph0gev[j];
			evt->stpph1gev[j]=stpph1gev[j];
			evt->stpph0mip[j]=stpph0mip[j];
			evt->stpph1mip[j]=stpph1mip[j];*/
		}


		int totplanes =0;
		int begu=1000;
		int begv=1000;
		int endu=0;
		int endv=0; 
		int nu=0;
		int nv=0;
		for(int i=0;i<1000;i++)
		{
			if(planes[i]>0)totplanes++;
			if(planes[i]==2)
			{
				nu++;
				begu=begu<i? begu:i;
				endu=endu<i? i :endu;
			}else if(planes[i]==3)
                        {
                                nv++;
                                begv=begv<i? begv:i;
				endv=endv<i?i:endv;
                        }
		}
		evt->plane.n=totplanes;
		evt->plane.nu=nu;
		evt->plane.nv=nv;
		evt->plane.beg=begu<begv?begu:begv;
		evt->plane.begu=begu;
		evt->plane.begv=begv;
		evt->plane.end=endu<endv?endv:endu;
		evt->plane.endu=endu;
		evt->plane.endv=endv;


	return goodevent;

}





int PETrimmer::CleanShower(NtpSRShower * shw, std::vector<const NtpSRStrip* > stp, int *pepass, float* stripe)
{
		int goodshower=1;
		int cutcount =0;
	
		//temporary vectors for saving good hits - better than erasing single entries in arrays at a time!
		std::vector<float> stpph0sigmap;
		std::vector<float> stpph0mip;
		std::vector<float> stpph0gev;
		std::vector<float> stpph1sigmap;
		std::vector<float> stpph1mip;
		std::vector<float> stpph1gev;
		std::vector<int> stpid;

			int us=0;
                        int vs=0;	
		for(int j=0;j<shw->nstpcnt;j++)
		{
		

	

			if(pepass[shw->stp[j]])
			{
			
				//keep it!
				stpid.push_back(shw->stp[j]);
				stpph0sigmap.push_back(shw->stpph0sigmap[j]);
				stpph1sigmap.push_back(shw->stpph1sigmap[j]);
				stpph0gev.push_back(shw->stpph0gev[j]);
				stpph1gev.push_back(shw->stpph1gev[j]);
				stpph0mip.push_back(shw->stpph0mip[j]);
				stpph1mip.push_back(shw->stpph1mip[j]);
				
				stripe[shw->stp[j]*3]= shw->stpph0sigmap[j]+shw->stpph1sigmap[j];
				stripe[1+shw->stp[j]*3]= shw->stpph0mip[j]+shw->stpph1mip[j];
				stripe[2+shw->stp[j]*3]= shw->stpph0gev[j]+shw->stpph1gev[j];

				if(stp[shw->stp[j]]->planeview==2)
					us++;
                                if(stp[shw->stp[j]]->planeview==3)
                                        vs++;


			}else{
				cutcount++;
			
				//adjust total energy of shower
				shw->ph.raw -= (stp[shw->stp[j]]->ph0.raw + stp[shw->stp[j]]->ph1.raw);
				shw->ph.pe -= (stp[shw->stp[j]]->ph0.pe + stp[shw->stp[j]]->ph1.pe);
				shw->ph.sigcor -= (stp[shw->stp[j]]->ph0.sigcor + stp[shw->stp[j]]->ph1.sigcor);
				shw->ph.siglin -= (stp[shw->stp[j]]->ph0.siglin + stp[shw->stp[j]]->ph1.siglin);
				shw->ph.sigmap -= (shw->stpph0sigmap[j] + shw->stpph1sigmap[j]);
				shw->ph.mip -= (shw->stpph0mip[j] + shw->stpph1mip[j]);
				shw->ph.gev -= (shw->stpph0gev[j] + shw->stpph1gev[j]);
	
			}

			
		}

                MSG("PETrimmer",Msg::kDebug)<<"Shower us "<<us << " vs "<< vs <<endl;
		if(us ==0 || vs==0)goodshower=0; //require && here as a test

		MSG("PETrimmer",Msg::kDebug)<<"Cut "<<cutcount<<" strips from shower!"<<endl;


		if (cutcount==0)return goodshower;  //we didn't cut anything, so continue	

		//rewrite event stp vectors
		for(int j=0;j<shw->nstpcnt;j++)
		{
			shw->stp[j]=-1;
			shw->stpph0sigmap[j]=0;
			shw->stpph1sigmap[j]=0;
			shw->stpph0gev[j]=0;
			shw->stpph1gev[j]=0;
			shw->stpph0mip[j]=0;
			shw->stpph1mip[j]=0;
		}


/*		delete[] shw->stp;
		delete[] shw->stpph0sigmap;
		delete[] shw->stpph1sigmap;
		delete[] shw->stpph0gev;
		delete[] shw->stpph1gev;
		delete[] shw->stpph0mip;
		delete[] shw->stpph1mip;
		
		shw->stp = new Int_t[stpid.size()];
		shw->stpph0sigmap = new Float_t[stpid.size()];
		shw->stpph1sigmap = new Float_t[stpid.size()];
		shw->stpph0gev = new Float_t[stpid.size()];
		shw->stpph1gev = new Float_t[stpid.size()];
		shw->stpph0mip = new Float_t[stpid.size()];
		shw->stpph1mip = new Float_t[stpid.size()];	//these arrays dont seem to be in the ntpst file!
		
*/		
		
		
		shw->nstpcnt=stpid.size();
		shw->nstrip=stpid.size();  ///incorrect... we did not change array size.. workaround for other mods...

		for(int j=0;j<(int)stpid.size();j++)
		{
			shw->stp[j]=-1;
			if(stpid[j]>-1 && stpid[j] < (int)stp.size())	
				shw->stp[j]=stpid[j];
			shw->stpph0sigmap[j]=stpph0sigmap[j];
			shw->stpph1sigmap[j]=stpph1sigmap[j];
			shw->stpph0gev[j]=stpph0gev[j];
			shw->stpph1gev[j]=stpph1gev[j];
			shw->stpph0mip[j]=stpph0mip[j];
			shw->stpph1mip[j]=stpph1mip[j];
		}

	return goodshower;

}





int PETrimmer::CleanTrack(NtpSRTrack * trk, std::vector<const NtpSRStrip* > stp,  int *pepass, float* stripe)
{
		int cutcount =0;
int goodtrack=1;	
		//temporary vectors for saving good hits - better than erasing single entries in arrays at a time!
		std::vector<float> stpph0sigmap;
		std::vector<float> stpph0mip;
		std::vector<float> stpph0gev;
		std::vector<float> stpph1sigmap;
		std::vector<float> stpph1mip;
		std::vector<float> stpph1gev;
		std::vector<int> stpid;


           		  int us=0;
                        int vs=0;
           		 	
		for(int j=0;j<trk->nstpcnt;j++)
		{
		

			if(trk->stp[j]>0)	
			if(pepass[trk->stp[j]])
			{
			
				//keep it!
				stpid.push_back(trk->stp[j]);
				stpph0sigmap.push_back(trk->stpph0sigmap[j]);
				stpph1sigmap.push_back(trk->stpph1sigmap[j]);
				stpph0gev.push_back(trk->stpph0gev[j]);
				stpph1gev.push_back(trk->stpph1gev[j]);
				stpph0mip.push_back(trk->stpph0mip[j]);
				stpph1mip.push_back(trk->stpph1mip[j]);
				
				
				stripe[trk->stp[j]*3]= trk->stpph0sigmap[j]+trk->stpph1sigmap[j];
				stripe[1+trk->stp[j]*3]= trk->stpph0mip[j]+trk->stpph1mip[j];
				stripe[2+trk->stp[j]*3]= trk->stpph0gev[j]+trk->stpph1gev[j];


		                if(stp[trk->stp[j]]->planeview==2)
                                        us++;
                                if(stp[trk->stp[j]]->planeview==3)
                                        vs++;


				

			}else{
				cutcount++;	

				//adjust total energy of track // we dont use these values of the track (usually...)
				trk->ph.raw -= (stp[trk->stp[j]]->ph0.raw + stp[trk->stp[j]]->ph1.raw);
				trk->ph.pe -= (stp[trk->stp[j]]->ph0.pe + stp[trk->stp[j]]->ph1.pe);
				trk->ph.sigcor -= (stp[trk->stp[j]]->ph0.sigcor + stp[trk->stp[j]]->ph1.sigcor);
				trk->ph.siglin -= (stp[trk->stp[j]]->ph0.siglin + stp[trk->stp[j]]->ph1.siglin);
				trk->ph.sigmap -= (trk->stpph0sigmap[j] + trk->stpph1sigmap[j]);
				trk->ph.mip -= (trk->stpph0mip[j] + trk->stpph1mip[j]);
				trk->ph.gev -= (trk->stpph0gev[j] + trk->stpph1gev[j]);


			}
		}

 		MSG("PETrimmer",Msg::kDebug)<<"Track us "<<us << " vs "<< vs <<endl;
                if(us ==0 || vs==0)goodtrack=0;
           
		MSG("PETrimmer",Msg::kDebug)<<"Cut "<<cutcount<<" strips from track!"<<endl;
		if(!goodtrack)MSG("PETrimmer",Msg::kDebug)<<"Track marked as bad!"<<endl;

		if (cutcount==0 )return goodtrack;  //we didn't cut anything, so continue	
		
		//rewrite event stp vectors
		for(int j=0;j<trk->nstpcnt;j++)
		{
			trk->stp[j]=-1;
			trk->stpph0sigmap[j]=0;
			trk->stpph1sigmap[j]=0;
			trk->stpph0gev[j]=0;
			trk->stpph1gev[j]=0;
			trk->stpph0mip[j]=0;
			trk->stpph1mip[j]=0;
		}

/*
		delete[] trk->stp;
		delete[] trk->stpph0sigmap;
		delete[] trk->stpph1sigmap;
		delete[] trk->stpph0gev;
		delete[] trk->stpph1gev;
		delete[] trk->stpph0mip;
		delete[] trk->stpph1mip;
		
		trk->stp = new Int_t[stpid.size()];
		trk->stpph0sigmap = new Float_t[stpid.size()];
		trk->stpph1sigmap = new Float_t[stpid.size()];
		trk->stpph0gev = new Float_t[stpid.size()];
		trk->stpph1gev = new Float_t[stpid.size()];
		trk->stpph0mip = new Float_t[stpid.size()];
		trk->stpph1mip = new Float_t[stpid.size()];	//these arrays dont seem to be in the ntpst file!
*/		
		
		
		
		trk->nstpcnt=stpid.size();
		trk->nstrip=stpid.size();  ///incorrect... we did not change array size.. workaround for other mods...

		for(int j=0;j<(int)stpid.size();j++)
		{
		        trk->stp[j]=-1;
                        if(stpid[j]>-1 && stpid[j] < (int)stp.size())   	
				trk->stp[j]=stpid[j];
			trk->stpph0sigmap[j]=stpph0sigmap[j];
			trk->stpph1sigmap[j]=stpph1sigmap[j];
			trk->stpph0gev[j]=stpph0gev[j];
			trk->stpph1gev[j]=stpph1gev[j];
			trk->stpph0mip[j]=stpph0mip[j];
			trk->stpph1mip[j]=stpph1mip[j];
		}



	return goodtrack; //track is bad only if it is empty


}



////////////////////////////////////////////////////////////////////////
void PETrimmer::EndJob()
{
}

const Registry& PETrimmer::DefaultConfig() const
{
//======================================================================
// Supply the default configuration for the module
//======================================================================
   MSG("PETrimmer",Msg::kDebug)<<"In Trimmer::DefaultConfig"<<endl;

  
  static Registry r;
  std::string name = this->JobCModule::GetName();
  name += ".config.default";
  r.SetName(name.c_str());

  r.UnLockValues();
  r.Set("PECut", 0.0);
  r.Set("updateEventEnergy",1);               

  r.LockValues();             
               
                                                                              
  return r;
}

void PETrimmer::Config(const Registry& r)
{
//======================================================================
// Configure the module given the Registry r
//======================================================================
  MSG("PETrimmer",Msg::kDebug)<<"In Trimmer::Config"<<endl;
                                                                               
  bool islocked = this->GetConfig().ValuesLocked();
  if (islocked) this->GetConfig().UnLockValues();
 
 
  double fmps;
  if(r.Get("PECut", fmps)) {pecut = fmps;}

  int imps;
  if(r.Get("updateEventEnergy", imps)) {updateEventEnergy = imps;}
  if (islocked) this->GetConfig().LockValues();

}



void PETrimmer::DumpIt(NtpStRecord * str)
{
	std::vector< NtpSREvent* > evts;
	std::vector< NtpSRShower* > shw;
	std::vector< NtpSRTrack* > trk;

	TClonesArray* tc = str->evt;	
	for(int i=0;i<tc->GetEntries();i++)
		evts.push_back((NtpSREvent*)tc->At(i));
		
	TClonesArray* ts = str->shw;	
	for(int i=0;i<ts->GetEntries();i++)
		shw.push_back((NtpSRShower*)ts->At(i));
		
	TClonesArray* tt = str->trk;	
	for(int i=0;i<tt->GetEntries();i++)
		trk.push_back((NtpSRTrack*)tt->At(i));

	std::vector<const NtpSRStrip* > stp = str->GetStrips();
		
	
  //iterate over showers, print(stripindx, pe)

	
	for(int i=0;i<(int)shw.size();i++)
	  {
	    ostringstream os;

	    os <<"Shower "<<i<<": ";
	    for(int j=0;j<shw[i]->nstrip;j++)
	      {

		if(i<0 || i > (int)shw.size()-1)
                {
                        os << "bad shower index!";
			continue;
		}
                int strip = shw[i]->stp[j]; 
                if(strip<0)
		{
                	os<<"("<<shw[i]->stp[j]<<",NO STRIP HERE....) ";                
			MSG("PETrimmer",Msg::kError)<<"Shower "<<i<<" is missing a strip at index "<<j<< " strip index "<<shw[i]->stp[j]<<endl;
		}
		else
		os<<"("<<shw[i]->stp[j]<<","<<stp[shw[i]->stp[j]]->ph0.pe+stp[shw[i]->stp[j]]->ph1.pe<<") ";
	      }
	    os <<"\n";		
	    MSG("PETrimmer",Msg::kDebug)<<os.str()<<endl;

	
	  }

	for(int i=0;i<(int)trk.size();i++)
	  {
	    ostringstream os;

	    os <<"Track "<<i<<": ";
	    for(int j=0;j<trk[i]->nstrip;j++)
	      {
                if(i<0 || i > (int)trk.size()-1)
                {
                        os << "bad track index!";
			continue;
		}	
                int strip = trk[i]->stp[j];
                if(strip<0)
                os<<"("<<trk[i]->stp[j]<<",NO STRIP HERE...) ";
                else
		os<<"("<<trk[i]->stp[j]<<","<<stp[trk[i]->stp[j]]->ph0.pe+stp[trk[i]->stp[j]]->ph1.pe<<") ";
	      }
	    os <<"\n";		
	    MSG("PETrimmer",Msg::kDebug)<<os.str()<<endl;
	  }


	for(int i=0;i<(int)evts.size();i++)
	  {
	    ostringstream os;

	    os <<"Event "<<i<<": ";
	    for(int j=0;j<evts[i]->nstrip;j++)
	      {

                if(i<0 || i > (int)evts.size()-1)
		{
			os << "bad event index!";
			continue;
                }
		int strip = evts[i]->stp[j];
                if(strip<0)
                 os<<"("<<evts[i]->stp[j]<<",NO STRIP HERE...) ";
                else
		os<<"("<<evts[i]->stp[j]<<","<<stp[evts[i]->stp[j]]->ph0.pe+stp[evts[i]->stp[j]]->ph1.pe<<") ";
	      }
	    os <<"\n";		


	    os<<"\tShowers: ";
	    for(int j=0;j<evts[i]->nshower;j++)
	      os<<evts[i]->shw[j]<<" ";

	    os<<"\tTracks: ";
	    for(int j=0;j<evts[i]->ntrack;j++)
	      os<<evts[i]->trk[j]<<" ";




	    MSG("PETrimmer",Msg::kDebug)<<os.str()<<endl;
	  }


}
