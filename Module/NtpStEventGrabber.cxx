#include "TFile.h"
#include "Conventions/Detector.h"
#include "NueAna/NueRecord.h"
#include "MessageService/MsgService.h"
#include "MinosObjectMap/MomNavigator.h"
#include "MessageService/MsgService.h"
#include "JobControl/JobCModuleRegistry.h" // For JOBMODULE macro
#include "HistMan/HistMan.h"
#include <fstream>
#include <iostream>
#include "AnalysisNtuples/ANtpDefaultValue.h"
#include "TClass.h"
#include "TDataType.h"
#include "TDataMember.h"
#include "TRealData.h"
#include <string>
#include "NtpStEventGrabber.h"
#include "NueAna/NueAnalysisCuts.h"


#include "TLeaf.h"
#include <map>
//......................................................................

NtpStEventGrabber::NtpStEventGrabber()
{
        outfile=0;
        outtree=0;
	outtreeMRCC=0;
        counter=0;
	doMRCC=0;
        ntpst=new NtpStRecord();
	ntpmr=new NtpMRRecord();
}

//......................................................................

NtpStEventGrabber::~NtpStEventGrabber()
{}

void NtpStEventGrabber::Gather(string filelist, string runlist, string outfile, int doMRCC)
{
	counter=0;
	LoadFileList(filelist);
	LoadRunList(runlist);

	this->outfile=new TFile(outfile.c_str(),"RECREATE");
	outtree=new TTree("NtpSt","NtpSt");
	outtree->Branch("NtpSt",&ntpst);	

	this->doMRCC=doMRCC;
	if(doMRCC)
	{
		outtreeMRCC=new TTree("NtpMR","NtpMR");
		outtreeMRCC->Branch("NtpMR",&ntpmr);
	}

	DoGather();
	
	this->outfile->cd();
	outtree->Write();
	if(doMRCC)outtreeMRCC->Write();
	this->outfile->Close();

	printf("got %d of %lu\n",counter, (long unsigned int)this->runlist.size());

}

void NtpStEventGrabber::DoGather()
{
	for(int i=0;i<(int)runlist.size();i++)
	{
		int run=runlist[i];
		
		string file="";
		for(int j=0;j<(int)filerunlist.size();j++)
		{
			if(filerunlist[j]==run)
			{
				file=filelist[j];
				break;
			}
		}

		TChain * c = new TChain("NtpSt");
		TChain *m=0;
		if(doMRCC)m = new TChain("NtpMR");

		int a = file.find("/pnfs/minos/");

		if(a>-1)file=file.substr(11); //cut out /pnfs/minos/

		file="dcap://fndca1.fnal.gov:24125/pnfs/fnal.gov/usr/minos/"+file;
		printf("%s\n",file.c_str());
		c->Add(file.c_str());

		if(doMRCC)
		{
			m->Add(file.c_str());
			m->SetBranchAddress("NtpMRRecord",&ntpmr);	
		}

//TChain::BuildIndex is broken in this version of root (all run/snarl entries are read as 0)
//   so do it manually

		c->SetBranchStatus("*",0);
		c->SetBranchStatus("fHeader.*",1);

		map<std::pair<int,int>, int> index;

		for(int j=0;c->GetEntry(j);j++)
		{
			int r=(int)c->GetLeaf("fHeader.fRun")->GetValue();
			int s=(int)c->GetLeaf("fHeader.fSnarl")->GetValue();
			index[std::pair<int,int>(r,s)]=j;
		}
		printf("index built\n");

		c->SetBranchStatus("*",1);
		c->ResetBranchAddresses();
		c->SetBranchAddress("NtpStRecord",&ntpst);

               while(1)
                {
                        int found=index[std::pair<int,int>(runlist[i],snarllist[i])];
                        if(found>0)
                        {
                	 	c->GetEntry(found);
		                outtree->Fill();

				if(doMRCC)
				{
					m->GetEntry(found);
					outtreeMRCC->Fill();
				}

                                counter++;

                                printf("just got snarl %d  (%d of %lu)\n",ntpst->GetHeader().GetSnarl(),counter,(long unsigned int)runlist.size());
 
                                //do we have more from this file?
                                if(runlist[i+1]==runlist[i])
                                        i++;
                                else
                                        break;


                        }else break;                    
                }




/*		c->BuildIndex("fHeader.fRun","fHeader.fSnarl");

		while(1)
		{
			int found=c->GetEntryWithIndex(runlist[i],snarllist[i]);
			if(found>0)
			{
				outtree->Fill();
				counter++;

                                printf("just got snarl %d  (%d of %d)\n",ntpst->GetHeader().GetSnarl(),counter,runlist.size());
 
                                //do we have more from this file?
                                if(runlist[i+1]==runlist[i])
                                        i++;
                                else
                                        break;


			}else break;			
		}


*/	
//		SetFastBranchStatus(c);
/*
		for(int k=0;c->GetEntry(k);k++)
		{
			if(ntpst->GetHeader().GetSnarl()==snarllist[i])
			{
  //                              c->SetBranchStatus("*",1);
//				c->ResetBranchAddresses();
//				c->SetBranchAddress("NtpSt",&ntpst);
//				c->GetEntry(k);
				//storeit....
				outtree->Fill();
				counter++;

				printf("just got snarl %d  (%d of %d)\n",ntpst->GetHeader().GetSnarl(),counter,runlist.size());
		//		SetFastBranchStatus(c);

				//do we have more from this file?
				if(runlist[i+1]==runlist[i])
					i++;
				else
					break;

			
			}
		}
*/




		delete c;

	}

}

void NtpStEventGrabber::SetFastBranchStatus(TChain *c)
{
	c->SetBranchStatus("*",1);
	c->SetBranchStatus("evthdr.*",0);
	c->SetBranchStatus("vetohdr.*",0);
        c->SetBranchStatus("crhdr.*",0);
        c->SetBranchStatus("dmxstatus.*",0);
        c->SetBranchStatus("detstatus.*",0);
        c->SetBranchStatus("calstatus.*",0);
        c->SetBranchStatus("dataquality.*",0);
        c->SetBranchStatus("mchdr.*",0);
        c->SetBranchStatus("photon.*",0);
        c->SetBranchStatus("detsim.*",0);
        c->SetBranchStatus("vetostp.*",0);
        c->SetBranchStatus("vetoexp.*",0);
        c->SetBranchStatus("deadchips.*",0);
        c->SetBranchStatus("stp.*",0);
        c->SetBranchStatus("slc.*",0);
        c->SetBranchStatus("clu.*",0);
        c->SetBranchStatus("shw.*",0);
        c->SetBranchStatus("trk.*",0);
        c->SetBranchStatus("evt.*",0);
        c->SetBranchStatus("mc.*",0);
        c->SetBranchStatus("stdhep.*",0);
        c->SetBranchStatus("digihit.*",0);
        c->SetBranchStatus("thstp.*",0);
        c->SetBranchStatus("thslc.*",0);
        c->SetBranchStatus("thshw.*",0);
        c->SetBranchStatus("thtrk.*",0);
        c->SetBranchStatus("thevt.*",0);
}

void NtpStEventGrabber::LoadFileList(string fl)
{
	filelist.clear();
	filerunlist.clear();

	string f="";

	ifstream file(fl.c_str());
	if(file.is_open())
	{
		while(!file.eof())
		{
			getline(file,f);

			if(f.length()<2)continue;
			filelist.push_back(f);

			//extract run 
			int run=0;

			run = atoi(f.substr(f.find_last_of("/")+2,8).c_str());
			filerunlist.push_back(run);

//			printf("%d %s\n",run,f.c_str());
//
		}
		file.close();
	}

	printf("loaded %lu file entries\n",(long unsigned int)filelist.size());

}



void NtpStEventGrabber::LoadRunList(string rl)
{

	runlist.clear();
	snarllist.clear();
	eventlist.clear();

        string f="";

        ifstream file(rl.c_str());
        if(file.is_open())
        {
                while(!file.eof())
                {
                        getline(file,f);

			if(f.length()<2)continue;
                    	int p1 = f.find(" ");
			int p2 = f.find(" ",p1+1);

			runlist.push_back(atoi(f.substr(0,p1).c_str()));
			snarllist.push_back(atoi(f.substr(p1,p2).c_str()));
			eventlist.push_back(atoi(f.substr(p2).c_str()));
                }
                file.close();
        }

        printf("loaded %lu request entries\n",(long unsigned int)runlist.size());


}



