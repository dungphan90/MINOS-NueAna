#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "NueAna/NueRecord.h"

using namespace std;

void PlotParA(const char* fname)
{
   TFile *f=new TFile(fname);
   TTree *tree = (TTree *)f->Get("ana_nue");

   NueRecord *nr = new NueRecord();

   tree->SetBranchAddress("NueRecord",&nr);

   float NENTRIES = tree->GetEntries();

   TFile out("testnue-nojc.root","RECREATE");
   TH1F *hpara=new TH1F("hpara","hpara",100,0,100);

   cout<<"Running over "<<NENTRIES<<endl;
   for(int z=0;z<NENTRIES;z++){
      if(z%1000==0){
	 cout<<"On entry "<<z<<endl;
      }

      tree->GetEntry(z);

      hpara->Fill(nr->shwfit.par_a);
   }
			   
   out.Write();
   out.Close();
}
