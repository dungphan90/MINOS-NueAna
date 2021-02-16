#include "TH1D.h"
#include "TFile.h"


void T()
{
	TFile *fout = TFile::Open("out.root","RECREATE");
	fout->cd();
	TH1D * h = new TH1D("h","h",100,0,1);
	h->Write();
	
	fout->Close();
}
