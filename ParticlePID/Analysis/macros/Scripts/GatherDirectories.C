#include "TDirectory.h"
#include "TFile.h"
#include "TObject.h"
#include <string>

void GatherDirectories(){

	TFile *fout = TFile::Open("merged.root","RECREATE");
	
		
	DoFile("far_mc_standard_normal/out.root",fout);
        DoFile("far_mc_standard_ParticlePID/out.root",fout);
        DoFile("far_mc_mrcc_normal/out.root",fout);
        DoFile("far_mc_mrcc_ParticlePID/out.root",fout);


 //       DoFile("near_mc_standard_hornOFF_normal/out.root",fout);
 //       DoFile("near_mc_standard_hornOFF_ParticlePID/out.root",fout);
        DoFile("near_mc_standard_hornON_normal/out.root",fout);
        DoFile("near_mc_standard_hornON_ParticlePID/out.root",fout);
 //       DoFile("near_mc_mrcc_hornOFF_normal/out.root",fout);
        DoFile("near_mc_mrcc_hornON_normal/out.root",fout);
        DoFile("near_mc_mrcc_hornON_ParticlePID/out.root",fout);


   //     DoFile("near_data_standard_hornOFF_normal/out.root",fout);
  //      DoFile("near_data_standard_hornOFF_ParticlePID/out.root",fout);
        DoFile("near_data_standard_hornON_normal/out.root",fout);
        DoFile("near_data_standard_hornON_ParticlePID/out.root",fout);
  //      DoFile("near_data_mrcc_hornOFF_normal/out.root",fout);
        DoFile("near_data_mrcc_hornON_normal/out.root",fout);
        DoFile("near_data_mrcc_hornON_ParticlePID/out.root",fout);

	fout->Write();
	fout->Close();

}


void DoFile(string file, TFile * fout)
{
	TFile *f = TFile::Open(file.c_str());
    	if(f)
    	{
    		CopyDir(fout->GetDirectory("/"),f->GetDirectory("/"));
		f->Close();
	}
	
	fout->Flush();
}

void CopyDir(TDirectory *dw, TDirectory *di)
{
	di->ReadAll();
	TList *l = di->GetList();
	for(int i=0;i<l->GetEntries();i++)
	{
		TObject *t = l->At(i);
		if(t->InheritsFrom("TDirectory"))
		{
			TDirectory * tw = dw->GetDirectory(t->GetName());
			if(!tw)
			{
				dw->mkdir(t->GetName());
				tw=dw->GetDirectory(t->GetName());
	

			}
			CopyDir(tw,di->GetDirectory(t->GetName()));
		} else if(t->InheritsFrom("TH1")){

			dw->cd();
			t->Clone();
		
		}

	}

}






