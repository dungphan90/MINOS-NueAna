#include "NueAna/ParticlePID/ParticleAna/ParticleTrimmer.h"
#include "TTree.h"



#include <stdio.h>

using std::cout;


ParticleTrimmer::ParticleTrimmer()
{
	chain_precord=new TChain("PA");
	chain_pot=new TChain("pottree");
		
	precord=new PRecord();
	pot=new POT();
	totalpot=new POT();
	
	chain_precord->SetBranchAddress("PRecord",&precord);
	chain_pot->SetBranchAddress("ParticlePOT",&pot);	
	totalpot->Reset();
	
	out=0;
	
}

ParticleTrimmer::~ParticleTrimmer()
{
	out->Close();
	
}

void ParticleTrimmer::AddFiles(string file)
{


	int added = chain_precord->Add(file.c_str());
	chain_pot->Add(file.c_str());
	cout <<"Adding "<<added<<": "<<file<<"\n";
}

void ParticleTrimmer::RunTrimmer()
{





	string file="";
	 if(chain_precord->GetEntries() > 0){
        chain_precord->GetEntry(0);
        file = chain_precord->GetFile()->GetName();
     }
     string minifile = file.substr(file.find_last_of("/")+1, file.find_last_of(".root")-file.find_last_of("/") - 5);
     minifile += "-Trim.root";
  
     outfile = minifile;

  if(outfile == "-Trim.root"){
     cout<<"No input file found\n"; 
     return;
  }
                                                                         
  cout<<"Setting output to "<<outfile<<"\n";
  
  out=new TFile(outfile.c_str(),"RECREATE");

	trimPOT();
	trimPRecord();

}
		
		
void ParticleTrimmer::trimPOT()
{
	printf("\nTrimming POTs\n");

	nearFileTypeCount=0;
	farFileTypeCount[0]=0;
	farFileTypeCount[1]=0;
	farFileTypeCount[2]=0;
	farFileTypeCount[3]=0;
	farFileTypeCount[4]=0;	
	
  TTree *tree = new TTree("pottree","pottree");
  TTree::SetBranchStyle(1);
  TBranch* br = tree->Branch("ParticlePOT", &totalpot );
  br->SetAutoDelete(kFALSE);
  
  int ent = (int)chain_pot->GetEntries();
  for(int i=0;i<ent;i++)
  {
      if(i%10000==0) cout << 100*float(i)/float(ent)
                       << "% done\n";
  	chain_pot->GetEntry(i);
  	
	string file="";
	file = chain_pot->GetFile()->GetName();
    
    string minifile = file.substr(file.find_last_of("/")+1, file.find_last_of(".root")-file.find_last_of("/") - 5);
    
    string filetype = minifile.substr(minifile.find_first_of("PO-")+3,4);
    printf("file type %s\n",filetype.c_str());	
  	
  	
  	string det = filetype.substr(0,1);
  	
  	string beam = filetype.substr(3,4);
  	
  	printf("%s %s\n",det.c_str(),beam.c_str());
  	
  	if(det=="f")
  	{
  		int beamtype=-1;
  		sscanf(beam.c_str(),"%d",&beamtype);
  		farFileTypeCount[beamtype]++;
  	}else if(det=="n" || det=="N" || det=="F") //far data is treated as near... (no adjustment)
  	{
  		nearFileTypeCount++;
  	}
  	
  	
  	printf("%d %f\n",i,pot->pot);
  	totalpot->pot+=pot->pot;
  	totalpot->nruns+=pot->nruns;
  	totalpot->nsnarls+=pot->nsnarls;
  	if(i==0)totalpot->beamtype=pot->beamtype;
  }         
 
 	int farcnt=0;
	for(int i=0;i<5;i++) 
 		farcnt+=farFileTypeCount[i];
 		
 	if(farcnt && nearFileTypeCount)
 	{
 		printf("attempting to trim files from both detectors!\n");
 		exit(1);
 	}
 	
 	//don't be fancy.. require the same number of each far type
 	if(farFileTypeCount[0]!=farFileTypeCount[3] || farFileTypeCount[0] != farFileTypeCount[4] 
 			|| farFileTypeCount[3] != farFileTypeCount[4])
 	{
 		printf("need same number of each file type in far!\n");
 		exit(1);
 	}
 	
 	if(farFileTypeCount[0])totalpot->pot/=farFileTypeCount[0];
  
  tree->Fill();
  out->cd();
  tree->Write();
  
                        
}

void ParticleTrimmer::trimPRecord()
{
	printf("\nTrimming PRecords...\n");
	
	
	//calculate the adjustment for totbeamweight based on the number of files
	
		double eventTypeAdj[5];
		for(int i=0;i<6;i++)eventTypeAdj[i]=1.;
		 
		if(nearFileTypeCount)
		{
			//do nothing...
		}else if(farFileTypeCount[0])
		{
			eventTypeAdj[0]=1./3.;
		}else{
			printf("no files found?\n");
			exit(1);
		}
	
	
	//
	
  out->cd();
  TTree *tree = new TTree("PA","PA");
  TTree::SetBranchStyle(1);
  TBranch* br = tree->Branch("PRecord", &precord );
  br->SetAutoDelete(kFALSE);
  
  int ent = (int)chain_precord->GetEntries();

	printf("%d total entries\n",ent);

	int totalpassed=0;

  for(int i=0;i<ent;i++)
  {


    if(i%10000==0)
	{ cout << 100*float(i)/float(ent)
                       << "% done\n";
  		tree->AutoSave("SaveSelf");
	}

	chain_precord->GetEntry(i);
  	
  	//some trimming
  	
  	if(precord->event.inFiducial!=1)continue;
  	if(precord->event.contained!=1)continue;
  	
  	if(precord->mctrue.type<0 || precord->mctrue.type>4)continue;
  	
  	precord->mctrue.totbeamweight*=eventTypeAdj[precord->mctrue.type];
  	
  	tree->Fill();
	totalpassed++;
  }         
  

   tree->AutoSave("SaveSelf");
printf("kept %d of %d entries\n",totalpassed,ent);
	//dups trees?
  //tree->Write();
  


}




