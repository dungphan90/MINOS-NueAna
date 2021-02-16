#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"

#include "TObjArray.h"
#include "TChainElement.h"

#include "TH1.h"

#include <string>
#include <vector>

#include "DirectoryHelpers.C"
using namespace DirectoryHelpers;

using namespace std;




//Stage1 takes in a TFile that has stage0 histograms already populated..



class Stage1
{

	public:
		Stage1();
		
				
		void SetInputFile(string i){infileString=i;};
		void SetOutputFile(string o){outfileString=o;};
		void SetOverwrite(int i=0){overwrite=i;};
		
		void MakeEfficiencies();
		void MakeEfficiencies(TDirectory *saveDir, TDirectory *dirNum, TDirectory *dirDen);

		void MakePurities();
		void MakePurities(TDirectory *saveDir, TDirectory *dir);
	
		void PrintConfigs();
		
		void Run();
	
	private:
		string infileString;
		string outfileString;
		TFile * infile;
		TFile * outfile;
		int overwrite;
		
		void PrepareFiles();
		
		//get a directory 0: requires it to exist 1: creates  2: creates clean (removes existing content)
		//TDirectory * GetDirectory(TFile *f, string dir, int create=0);
	//	TDirectory * GetDirectory(TDirectory *f, string dir, int create=0);
		//void Tokenize(const string& str, vector<string>& tokens,
          //            const string& delimiters = "/");
};

void Stage1::MakeEfficiencies()
{
	//call MakeEfficiencies for specific directories....

	TDirectory * saveDir = 0;
	TDirectory * dden = 0;
	TDirectory * dnum = 0;

	TDirectory *workingDir = GetDirectory(infile,"/stage1/efficiencies/",2); //clear out old info...
	workingDir->SaveSelf();

	for(int det=0;det <2;det++)
	for(int horn=0;horn< (det==0?1:2);horn++)
	for(int mrcc=0;mrcc<2;mrcc++)
	for(int mc=0;mc<2;mc++)
	for(int rec=0;rec<2;rec++)
	for(int rn=0;rn<4;rn++)
	{

	
		string dirp="";
		
		switch(det)
		{
			case 0:dirp+="far/";break;
			case 1:dirp+="near/";break;
		};

		switch(mc)
		{
			case 0:dirp+="MC/";break;
			case 1:dirp+="data/";break;
		};

		switch(mrcc)
		{
			case 0:dirp+="standard/";break;
			case 1:dirp+="MRCC/";break;
		};

		if(det==1)
		switch(horn)
		{
			case 0:dirp+="horn_on/";break;
			case 1:dirp+="horn_off/";break;
		};

		switch(rec)
		{
			case 0:dirp+="normal/";break;
			case 1:dirp+="ParticlePID/";break;
		};
		
		switch(rn)
		{
			case 0:dirp+="Run1/";break;
			case 1:dirp+="Run2/";break;
			case 2:dirp+="Run3/";break;
			case 3:dirp+="All/";break;
		};
		
		
		
		
		string outdir="/stage1/efficiencies/"+dirp;
		
		string indir="/stage0/"+dirp;
		
		printf("in  %s\n",indir.c_str());
		printf("out %s\n",outdir.c_str());
		
	


		dden = GetDirectory(infile,indir+"fiducial");
		dnum = GetDirectory(infile,indir+"presel");	
		if(dden && dnum)saveDir = GetDirectory(outfile, outdir+"/presel__over__fiducial",2);
		MakeEfficiencies(saveDir, dnum, dden);

		if(rec==0)
		{
			dden = GetDirectory(infile,indir+"/fiducial");
			dnum = GetDirectory(infile,indir+"/ann11");	
			if(dden && dnum)saveDir = GetDirectory(outfile, outdir+"/ann11__over__fiducial",2);
			MakeEfficiencies(saveDir, dnum, dden);


			dden = GetDirectory(infile,indir+"/fiducial");
			dnum = GetDirectory(infile,indir+"/ann11_firebird");	
			if(dden && dnum)saveDir = GetDirectory(outfile, outdir+"/ann11_firebird__over__fiducial",2);
			MakeEfficiencies(saveDir, dnum, dden);
		}
	
	
		if(rec==1)
		{
			for(int pid=0;pid<6;pid++)
			{
				string pidN="";
				switch(pid)
				{
					case 0:pidN="pidA";break;
					case 1:pidN="pidB";break;
					case 2:pidN="pidC";break;
					case 3:pidN="pidD";break;
					case 4:pidN="pidE";break;
					case 5:pidN="pidF";break;
				};
			
				
				dden = GetDirectory(infile,indir+"/fiducial");
				dnum = GetDirectory(infile,indir+pidN);	
				if(dden && dnum)saveDir = GetDirectory(outfile, outdir+"/"+pidN+"__over__fiducial",2);
				MakeEfficiencies(saveDir, dnum, dden);
			
			}
		
		
		}
		
	}
	
	workingDir->SaveSelf();
}

void Stage1::MakeEfficiencies(TDirectory *saveDir, TDirectory *dirNum, TDirectory *dirDen)
{
	if(!dirNum || !dirDen || !saveDir)
	{
		//printf("cant read directory (makeefficiencies)\n");
		return;
	}

	//iterate over TH1's in dirNum

	dirNum->ReadAll();
	TList *l = dirNum->GetList();
	for(int i=0;i<l->GetEntries();i++)
	{
		TObject *t = l->At(i);
		if(t->InheritsFrom("TH1"))
		{
			
	
			TObject *td = dirDen->FindObjectAny(t->GetName());  //use the Any to ignore key cycles..
			if(!td)continue;

			string tname=t->GetName();	
			string cmp="recoE";
			int match=0;
			if(tname.substr(0,cmp.length())==cmp)match=1;
			if(!match)continue;

	
			TH1 * tnum = (TH1*)t->Clone();
			TH1 * tden = (TH1*)td->Clone();
			
			tden->SetDirectory(0);
			
			tnum->Divide(tden);
			
			tnum->SetDirectory(saveDir);
			
		
		}

	}
}


void Stage1::MakePurities()
{
	//call MakeEfficiencies for specific directories....

	TDirectory * saveDir = 0;
	TDirectory * dir = 0;	

	TDirectory * workingDir = GetDirectory(infile,"/stage1/purities/",2); //clear out old info...
	workingDir->SaveSelf();

	for(int det=0;det <2;det++)
	for(int horn=0;horn< (det==0?1:2);horn++)
	for(int mrcc=0;mrcc<2;mrcc++)
	for(int mc=0;mc<2;mc++)
	for(int rec=0;rec<2;rec++)
	for(int rn=0;rn<4;rn++)
	{

	
		string dirp="";
		
		switch(det)
		{
			case 0:dirp+="far/";break;
			case 1:dirp+="near/";break;
		};

		switch(mc)
		{
			case 0:dirp+="MC/";break;
			case 1:dirp+="data/";break;
		};

		switch(mrcc)
		{
			case 0:dirp+="standard/";break;
			case 1:dirp+="MRCC/";break;
		};

		if(det==1)
		switch(horn)
		{
			case 0:dirp+="horn_on/";break;
			case 1:dirp+="horn_off/";break;
		};

		switch(rec)
		{
			case 0:dirp+="normal/";break;
			case 1:dirp+="ParticlePID/";break;
		};
		
				
		switch(rn)
		{
			case 0:dirp+="Run1/";break;
			case 1:dirp+="Run2/";break;
			case 2:dirp+="Run3/";break;
			case 3:dirp+="All/";break;
		};
		
		
		string outdir="/stage1/purities/"+dirp;
		
		string indir="/stage0/"+dirp;
		
		//printf("in  %s\n",indir.c_str());
		//printf("out %s\n",outdir.c_str());
		

		dir = GetDirectory(infile,indir+"fiducial");		
		if(dir)
		{
			saveDir = GetDirectory(outfile,outdir+"fiducial",2);
			MakePurities(saveDir, dir);
		}
		
		
		dir = GetDirectory(infile,indir+"presel");
		if(dir)
		{
			saveDir = GetDirectory(outfile,outdir+"presel",2);
			MakePurities(saveDir, dir);
		}
		
		
		if(rec==0)
		{
			dir = GetDirectory(infile,indir+"ann11");
			if(dir)
			{
				saveDir = GetDirectory(outfile, outdir+"ann11",2);
				MakePurities(saveDir, dir);
			}
			
			
			dir = GetDirectory(infile,indir+"ann11_firebird");
			if(dir)
			{			
				saveDir = GetDirectory(outfile, outdir+"ann11_firebird",2);
				MakePurities(saveDir, dir);
			}
		}
		
		if(rec==1)
		{
			for(int pid=0;pid<6;pid++)
			{
				string pidN="";
				switch(pid)
				{
					case 0:pidN="pidA";break;
					case 1:pidN="pidB";break;
					case 2:pidN="pidC";break;
					case 3:pidN="pidD";break;
					case 4:pidN="pidE";break;
					case 5:pidN="pidF";break;
				};
			
				dir = GetDirectory(infile,indir+pidN);
				if(dir)
				{
					saveDir = GetDirectory(outfile, outdir+pidN,2);
					MakePurities(saveDir, dir);
				}		
			
			}
		
		
		}
	}

	workingDir->SaveSelf();

}

void Stage1::MakePurities(TDirectory *saveDir, TDirectory *dir)
{
	//iterate over TH1's in dir

	if(!dir)
	{
		//printf("cant read directory!\n");
		return;
	}


	dir->ReadAll();
	TList *l = dir->GetList();
	for(int i=0;i<l->GetEntries();i++)
	{
		TObject *t = l->At(i);
		if(t->InheritsFrom("TH1"))
		{
		
			//make sure it is an "all" variable
			vector<string> parts;
			Tokenize(t->GetName(),parts,"_");
			string endp = parts[(int)parts.size()-1];
			if(endp=="sig" || endp=="nc" || endp=="cc" || endp=="tau" || endp=="beam")continue;
			

			for(int i=0;i<5;i++)
			{	
				string numn="";
				switch(i)
				{
					case 0:numn=(string)t->GetName()+"_nc";break;
					case 1:numn=(string)t->GetName()+"_cc";break;
					case 2:numn=(string)t->GetName()+"_sig";break;
					case 3:numn=(string)t->GetName()+"_tau";break;
					case 4:numn=(string)t->GetName()+"_beam";break;
				};
				
				TObject *td = dir->FindObjectAny(numn.c_str());  //use the Any to ignore key cycles..
				if(!td)continue;

			
				string tname=t->GetName();	
				string cmp="recoE";
				int match=0;
				if(tname.substr(0,cmp.length())==cmp)match=1;
				if(!match)continue;
			
			
			
			
				TH1 * tden = (TH1*)t->Clone();
				TH1 * tnum = (TH1*)td->Clone();
			
				tden->SetDirectory(0);
				tnum->SetDirectory(0);

				tnum->Divide(tden);
			
				string newname = numn+"__over__"+t->GetName();
				tnum->SetName(newname.c_str());
				tnum->SetDirectory(saveDir);

			}
		
		}

	}
}

void Stage1::PrintConfigs()
{
	printf("configuration:\n");
	printf("input file %s\n",infileString.c_str());
	printf("output file %s\n",outfileString.c_str());
	
}



Stage1::Stage1()
{
	infileString="";
	outfileString="";
	overwrite=0;
}

void Stage1::PrepareFiles()
{
	if(infileString=="")
	{
		printf("specify an input file....\n");
		exit(1);
	}

	if(outfileString=="" || infileString==outfileString)
	{
		infile=TFile::Open(infileString.c_str(),"UPDATE","",9);
		outfile=infile;
	}else{
		infile=TFile::Open(infileString.c_str());
		if(!overwrite)
			outfile=TFile::Open(outfileString.c_str(),"UPDATE","",9);
		else
			outfile=TFile::Open(outfileString.c_str(),"RECREATE","",9);
	}
	
	if(!infile)
	{
		printf("error opening infile\n");
		exit(1);
	}
	
	if(!outfile)
	{
		printf("error opening outfile\n");
		exit(1);
	}

}

void Stage1::Run()
{

	PrepareFiles();
	
	
	printf("ME\n");
	MakeEfficiencies();

	printf("MP\n");
	MakePurities();

	printf("Done\n");
	
	outfile->Write();
	outfile->Close();
	
}


