#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"

#include "TObjArray.h"
#include "TChainElement.h"

#include "TH1.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TLegend.h"


#include <string>
#include <vector>

#include <math.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include "TStyle.h"


#include "DirectoryHelpers.C"
using namespace DirectoryHelpers;

using namespace std;




//Stage2 takes in a TFile that has stage1 histograms already populated..
//It will make "final" plots as well as output text and latex


class Stage2
{

	public:
		Stage2();
		
				
		void SetInputFile(string i){infileString=i;};
		void SetOutputFile(string o){outfileString=o;};
		void SetOverwrite(int i=0){overwrite=i;};
		
		void SetPlotDirectory(string pd){plotDirectory=pd;};
		
		void CountComponentsByCutLevelFar(double sys=0.07);
		void CountComponentsByCutLevelFar(double sys, TDirectory *base);
		void CountComponentsByCutLevelFar(double sys, ostringstream &sraw, ostringstream &slatex, TDirectory *d);

				
		void CountComponentsByCutLevelNear(double sys=0.07);
		void CountComponentsByCutLevelNear(double sys, TDirectory *data, TDirectory *mc);
		void CountComponentsByCutLevelNear(double sys, ostringstream &sraw, ostringstream &slatex, TDirectory *data, TDirectory *mc);

		
		void DrawParsFar();
		void DrawParsFar(TDirectory *d, double sigscale=1);
		TCanvas * DrawParsFar(double sigscale, string localPlotDirectory, string newpath, TDirectory *data, TDirectory *mc, string var, double xlow=0, double xhigh=0, int rebin=1, double lx=0.7, double ly=0.7, char * tname="", string desc="");
		

		void DrawParsCompareNear();
		void DrawParsCompareNear(TDirectory *data, TDirectory *mc);
		TCanvas * DrawParsCompareNear(string localPlotDirectory, string newpath, TDirectory *data, TDirectory *mc, string var, double xlow=0, double xhigh=0, int rebin=1, double lx=0.7, double ly=0.7, char * tname="", string desc="");

		void DrawSuperFOMs(double sys=0.07);
		void DrawSuperFOMs(TDirectory *d, string pid, double sys);

		void DoMRCCPredictionNear();
		void DoMRCCPredictionNear(TDirectory *dataMRCC, TDirectory *mcMRCC, TDirectory *data, TDirectory *mc);
		
		void DoPurEffPlot();
		
		void MRCCDoubleRatioPlot();
		void MRCCDoubleRatioPlot(TDirectory * data, TDirectory *mc, TDirectory * mrccdata, TDirectory * mrccmc, string title);
		
		void DrawRecoEByRun();
		void DrawRecoEByRun(TDirectory * data, TDirectory *mc, string cut, string title="");
		
		
		
		void PrintConfigs();
		
		void Run();
		void SaveAndClose();
	
	private:
		string infileString;
		string outfileString;
		string plotDirectory;
		TFile * infile;
		TFile * outfile;
		int overwrite;
		
		void PrepareFiles();
		
		string GetName(int code);
		int GetColor(int code);	
		string GetParameterName(int code);
		
		enum mcCodes
		{
			NC =0,
			CC =1,
			Sig=2,
			Nue=2,
			Tau=3,
			Beam=4,
			Data=10
		};

};


void Stage2::PrintConfigs()
{
	printf("configuration:\n");
	printf("input file %s\n",infileString.c_str());
	printf("output file %s\n",outfileString.c_str());
	
}

Stage2::Stage2()
{
	infileString="";
	outfileString="";
	overwrite=0;
	plotDirectory="";
}

void Stage2::PrepareFiles()
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

//////////////////////
///// DrawRecoEByRun

void Stage2::DrawRecoEByRun()
{
	DrawRecoEByRun(
		GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/normal/"),
		GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/normal/"),
		"ann11",
		"Near Horn On Standard Normal Reconstruction ANN11"
	);
	
	DrawRecoEByRun(
		GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/normal/"),
		GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/normal/"),
		"ann11",
		"Near Horn On MRCC Normal Reconstruction ANN11"
	);

	DrawRecoEByRun(
		GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/ParticlePID/"),
		GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/ParticlePID/"),
		"pidF",
		"Near Horn On Standard ParticlePID pidF"
	);
	
	DrawRecoEByRun(
		GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/ParticlePID/"),
		GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/ParticlePID/"),
		"pidF",
		"Near Horn On MRCC ParticlePID pidF"
	);


	//a test plot
	DrawRecoEByRun(
		GetDirectory(infile,"/stage0/near/data/MRCC/horn_on/normal/"),
		GetDirectory(infile,"/stage0/near/MC/MRCC/horn_on/normal/"),
		"pidF",
		"Near Horn On MRCC Normal Reconstruction ANN11"
	);



}

void Stage2::DrawRecoEByRun(TDirectory * data, TDirectory *mc, string cut, string title)
{

	if(!data || !mc)return;

	TDirectory *tmp;
	
	tmp=data->GetDirectory(("Run1/"+cut).c_str());
	TH1D * r_data_1 = (TH1D*)tmp->Get("recoE")->Clone();
	
	tmp=data->GetDirectory(("Run2/"+cut).c_str());
	TH1D * r_data_2 = (TH1D*)tmp->Get("recoE")->Clone();
	
	tmp=data->GetDirectory(("Run3/"+cut).c_str());
	TH1D * r_data_3 = (TH1D*)tmp->Get("recoE")->Clone();
	
	
	tmp=mc->GetDirectory(("Run1/"+cut).c_str());
	TH1D * r_mc_1 = (TH1D*)tmp->Get("recoE")->Clone();
	
	tmp=mc->GetDirectory(("Run2/"+cut).c_str());
	TH1D * r_mc_2 = (TH1D*)tmp->Get("recoE")->Clone();
	
	tmp=mc->GetDirectory(("Run3/"+cut).c_str());
	TH1D * r_mc_3 = (TH1D*)tmp->Get("recoE")->Clone();
	

	r_data_1->SetLineColor(1);	
	r_data_2->SetLineColor(2);
	r_data_3->SetLineColor(4);
	
	r_mc_1->SetLineColor(1);	
	r_mc_2->SetLineColor(2);
	r_mc_3->SetLineColor(4);	

	r_data_1->SetLineStyle(2);	
	r_data_2->SetLineStyle(2);
	r_data_3->SetLineStyle(2);
	
	
	gStyle->SetOptStat(0);
	
//	TCanvas *c1 = new TCanvas("c","",500,600);
	TCanvas *c1 = new TCanvas();
	c1->cd();
	
	r_mc_1->SetTitle(title.c_str());
	r_mc_1->SetXTitle("Reconstructed Energy (GeV)");
	r_mc_1->SetYTitle("Events / 1e19 POT");
	r_mc_1->GetYaxis()->SetTitleOffset(1.2);
	
	r_mc_1->Draw("h");
	r_mc_2->Draw("sameh");
	r_mc_3->Draw("sameh");
	
	r_data_1->Draw("sameh");
	r_data_2->Draw("sameh");
	r_data_3->Draw("sameh");


	TLegend * l = new TLegend(0.7,0.7,0.88,0.88,"","NDC");
	l->AddEntry(r_mc_1,"MC Run 1");
	l->AddEntry(r_mc_2,"MC Run 2");
	l->AddEntry(r_mc_3,"MC Run 3");
	
	l->AddEntry(r_data_1,"Data Run 1");
	l->AddEntry(r_data_2,"Data Run 2");
	l->AddEntry(r_data_3,"Data Run 3");
	
	l->Draw("same");
	

	vector<string>paths;
	Tokenize(mc->GetPath(),paths,"/");
	
	string localPlotDirectory=plotDirectory+"/"+paths[2]+"/compare/";
	for(int j=4;j<(int)paths.size();j++)
	{
		localPlotDirectory+=paths[j]+"/";
	}	
	
	if(plotDirectory!="")MakeTrueDirectory(localPlotDirectory+"/gif");		
	if(plotDirectory!="")MakeTrueDirectory(localPlotDirectory+"/eps");
	
	if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/eps/recoEByRun_"+cut+".eps").c_str());
	if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/gif/recoEByRun_"+cut+".gif").c_str());

	
}



///// DrawRecoEByRun
//////////////////////
 

//////////////////////
///// MRCCDoubleRatioPlot

void Stage2::MRCCDoubleRatioPlot()
{

	string run = "All";
	
	
	MRCCDoubleRatioPlot(
		GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/ParticlePID/"+run+"/presel"),
		GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/ParticlePID/"+run+"/presel"),
		GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/ParticlePID/"+run+"/presel"),
		GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/ParticlePID/"+run+"/presel"),
		"Near Standard ParticlePID Presel (Run1+Run2+Run3)/3"	
	);

	MRCCDoubleRatioPlot(
		GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/ParticlePID/"+run+"/pidF"),
		GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/ParticlePID/"+run+"/pidF"),
		GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/ParticlePID/"+run+"/pidF"),
		GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/ParticlePID/"+run+"/pidF"),
		"Near Standard ParticlePID pidF (Run1+Run2+Run3)/3"	
	);

	MRCCDoubleRatioPlot(
		GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/normal/"+run+"/presel"),
		GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/normal/"+run+"/presel"),
		GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/normal/"+run+"/presel"),
		GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/normal/"+run+"/presel"),
		"Near Standard Normal Presel (Run1+Run2+Run3)/3"	
	);

	MRCCDoubleRatioPlot(
		GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/normal/"+run+"/ann11"),
		GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/normal/"+run+"/ann11"),
		GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/normal/"+run+"/ann11"),
		GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/normal/"+run+"/ann11"),
		"Near Standard Normal ANN11 (Run1+Run2+Run3)/3"	
	);

}


void Stage2::MRCCDoubleRatioPlot(TDirectory * data, TDirectory *mc, TDirectory * mrccdata, TDirectory * mrccmc, string title)
{

	TCanvas *c1=new TCanvas();

	gStyle->SetOptStat(0);
	
	TH1D * d = (TH1D*)data->FindObjectAny("recoE")->Clone();
	TH1D * m = (TH1D*)mc->FindObjectAny("recoE")->Clone();
	TH1D * dmrcc = (TH1D*)mrccdata->FindObjectAny("recoE")->Clone();
	TH1D * mmrcc = (TH1D*)mrccmc->FindObjectAny("recoE")->Clone();
			

	dmrcc->Divide(mmrcc);
	
	d->Divide(m);
	
	dmrcc->Divide(d);
	
	dmrcc->SetTitle(title.c_str());
	dmrcc->SetXTitle("Reconstructed Energy (GeV)");
	dmrcc->SetYTitle(" (MRCC Data / MRCC MC) / (Data / MC)");
	dmrcc->GetYaxis()->SetTitleOffset(1.2);
	dmrcc->GetYaxis()->SetRangeUser(0.8,1.2);
	dmrcc->Draw();

		
	vector<string>paths;
	Tokenize(mrccdata->GetPath(),paths,"/");
	
	string localPlotDirectory=plotDirectory+"/"+paths[2]+"/prediction/";
	for(int j=4;j<(int)paths.size();j++)
	{
		localPlotDirectory+=paths[j]+"/";
	}	
	
	if(plotDirectory!="")MakeTrueDirectory(localPlotDirectory+"/gif");		
	if(plotDirectory!="")MakeTrueDirectory(localPlotDirectory+"/eps");
	
	printf("%s\n",localPlotDirectory.c_str());
	
	if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/eps/mrccDoubleRatio.eps").c_str());
	if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/gif/mrccDoubleRatio.gif").c_str());

	delete c1;

}


///// 
//////////////////////


//////////////////////
/////  

void Stage2::DoPurEffPlot()
{

	//quick for now...
	gStyle->SetOptStat(0);
	
	
	TDirectory *dirpur = infile->GetDirectory("/stage1/full/purities/far/MC/standard/ParticlePID/All/pidF/");
	TH1D * pur = (TH1D*)dirpur->FindObjectAny("recoE_sig__over__recoE");


	TDirectory *direff = infile->GetDirectory("/stage1/full/efficiencies/far/MC/standard/ParticlePID/All/pidF__over__fiducial/");
	TH1D * eff = (TH1D*)direff->FindObjectAny("recoE_sig");
	
	if(!pur || !eff)return;
	
	pur->SetLineColor(2);
	eff->SetLineColor(4);
	
	eff->SetXTitle("Reconstructed Energy (GeV)");
	eff->SetTitle("");
	
	TCanvas *c1= new TCanvas();
	c1->cd();
	eff->Draw("h");
	pur->Draw("hsame");

	
	TLegend * l = new TLegend(0.7,0.7,0.88,0.88,"","NDC");
	l->AddEntry(pur,"Purity");
	l->AddEntry(eff,"Efficiency");
	
	
	l->Draw();
	
	
	c1->SaveAs("pureff_pidF.gif");
	
	

}

/////  
//////////////////////




//////////////////////
/////  DoMRCCPredictionNear

void Stage2::DoMRCCPredictionNear()
{


	
	for(int i=0;i<4;i++)
	{
		string run="";
		switch(i)
		{
			case 0:run="Run1";break;
			case 1:run="Run2";break;
			case 2:run="Run3";break;
			case 3:run="All";break;
		}
/*
		DoMRCCPredictionNear(
			GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/ParticlePID/"+run+"/presel"),
			GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/ParticlePID/"+run+"/presel"),
			GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/ParticlePID/"+run+"/presel"),
			GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/ParticlePID/"+run+"/presel")
		);

		DoMRCCPredictionNear(
			GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/ParticlePID/"+run+"/pidF"),
			GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/ParticlePID/"+run+"/pidF"),
			GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/ParticlePID/"+run+"/pidF"),
			GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/ParticlePID/"+run+"/pidF")
		);
*/
		DoMRCCPredictionNear(
			GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/normal/"+run+"/ann11/mrcc_all"),
			GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/normal/"+run+"/ann11/mrcc_all"),
			GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/normal/"+run+"/ann11"),
			GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/normal/"+run+"/ann11")
		);

                DoMRCCPredictionNear(
                        GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/normal/"+run+"/ann14/mrcc_all"),
                        GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/normal/"+run+"/ann14/mrcc_all"),
                        GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/normal/"+run+"/ann14"),
                        GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/normal/"+run+"/ann14")
                );


                DoMRCCPredictionNear(
                        GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/normal/"+run+"/pidF/mrcc_all"),
                        GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/normal/"+run+"/pidF/mrcc_all"),
                        GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/normal/"+run+"/pidF"),
                        GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/normal/"+run+"/pidF")
                );

	}

}


void Stage2::DoMRCCPredictionNear(TDirectory *dataMRCC, TDirectory *mcMRCC, TDirectory *data, TDirectory *mc)
{

	TCanvas *c1 = new TCanvas();
	c1->cd();

		TH1D * m = (TH1D*)mcMRCC->FindObjectAny("recoE")->Clone();
	//	TH1D * m_nc = (TH1D*)mcMRCC->FindObjectAny("recoE_nc")->Clone();
	//	TH1D * m_cc = (TH1D*)mcMRCC->FindObjectAny("recoE_cc")->Clone();	
	//	TH1D * m_beam = (TH1D*)mcMRCC->FindObjectAny("recoE_beam")->Clone();							
		TH1D * m_data = (TH1D*)dataMRCC->FindObjectAny("recoE")->Clone();								

		m_data->Divide(m);
		
		TH1D * p = (TH1D*)mc->FindObjectAny("recoE")->Clone();
		TH1D * p_nc = (TH1D*)mc->FindObjectAny("recoE_nc")->Clone();
		TH1D * p_beam = (TH1D*)mc->FindObjectAny("recoE_beam")->Clone();							
		TH1D * p_data = (TH1D*)data->FindObjectAny("recoE")->Clone();			
		p_nc->Multiply(m_data);
		

		TH1D * p_cc = (TH1D*)p_data->Clone();
		p_cc->Add(p_nc,-1);
		p_cc->Add(p_beam,-1);
		
		
		gStyle->SetOptStat(0);
		
		p_data->SetTitle("");
		p_data->SetXTitle("Reconstructed Energy (GeV)");
		p_data->SetYTitle("Events / 1e19 POT");
		p_data->GetYaxis()->SetTitleOffset(1.2);
		
		p_data->Draw();
		
	//	p->SetLineColor(2);
	//	p->Draw("sameh");
		
		p_cc->SetLineColor(3);
		p_cc->Draw("sameh");
		
		p_beam->SetLineColor(4);
		p_beam->Draw("sameh");
		
		p_nc->SetLineColor(6);
		p_nc->Draw("sameh");
		
		
		TH1D * all = (TH1D*)p_nc->Clone();
		all->Add(p_cc);
		all->Add(p_beam);
		
		all->SetLineColor(2);
		all->Draw("same");
		
		
				
	 	double lx=0.7;
	 	double ly=0.7;
		TLegend *l=new TLegend(lx,ly,lx+0.18,ly+0.18,"","NDC");
		l->AddEntry(all,"All");
		l->AddEntry(p_data,"Data");
		l->AddEntry(p_nc,GetName(NC).c_str());
		l->AddEntry(p_cc,GetName(CC).c_str());
		l->AddEntry(p_beam,GetName(Beam).c_str());
		l->Draw("Same");
		
		
		printf("%s\n%s\n%s\n%s\n",dataMRCC->GetPath(),mcMRCC->GetPath(),data->GetPath(),mc->GetPath());
		printf("bg %f nc %f cc %f beam %f\n",all->Integral(),p_nc->Integral(),p_cc->Integral(),p_beam->Integral());
		
		
	vector<string>paths;
	Tokenize(dataMRCC->GetPath(),paths,"/");
	
	string localPlotDirectory=plotDirectory+"/"+paths[2]+"/prediction/";
	for(int j=4;j<(int)paths.size();j++)
	{
		localPlotDirectory+=paths[j]+"/";
	}	
	
	if(plotDirectory!="")MakeTrueDirectory(localPlotDirectory+"/gif");		
	if(plotDirectory!="")MakeTrueDirectory(localPlotDirectory+"/eps");
	
	printf("%s\n",localPlotDirectory.c_str());
	
	if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/eps/mrccPrediction.eps").c_str());
	if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/gif/mrccPrediction.gif").c_str());

		
}


/////  DoMRCCPredictionNear
//////////////////////


//////////////////////
/////


void Stage2::DrawParsCompareNear()
{

	for(int hoo=0;hoo<1;hoo++)  //0 2
	for(int rec=1;rec<2;rec++)  // 0 2
	{
		int maxj= (rec==0?7:4);
	for(int j=0;j<1;j++)  //0 maxj
	for(int i=0;i<4;i++) // 0 4
	{
	
		printf("%d %d %d %d\n",hoo,rec,j,i);
	
		string horn="";
		switch(hoo)
		{
			case 0:horn="horn_on";break;
			case 1:horn="horn_off";break;
		}		
	
		string reco="";
		switch(rec)
		{
			case 0:reco="ParticlePID";break;
			case 1:reco="normal";break;
		}
	
		string run="";	
		switch(i)
		{
			case 0:run="Run1";break;
			case 1:run="Run2";break;
			case 2:run="Run3";break;
			case 3:run="All";break;
		}

		string cut="";
		if(rec==0)
		switch(j)
		{
			case 0:cut="presel";break;
			case 1:cut="pidA";break;
			case 2:cut="pidB";break;
			case 3:cut="pidC";break;
			case 4:cut="pidD";break;
			case 5:cut="pidE";break;
			case 6:cut="pidF";break;
		}

		if(rec==1)
		switch(j)
		{
			case 0:cut="presel";break;
			case 1:cut="ann11";break;
			case 2:cut="ann11_firebird";break;
			case 3:cut="pidF";break;
		}

		DrawParsCompareNear(GetDirectory(infile,"/stage0/full/near/data/standard/"+horn+"/"+reco+"/"+run+"/"+cut+"/"),
			GetDirectory(infile,"/stage0/full/near/MC/standard/"+horn+"/"+reco+"/"+run+"/"+cut+"/"));
		
		DrawParsCompareNear(GetDirectory(infile,"/stage0/full/near/data/MRCC/"+horn+"/"+reco+"/"+run+"/"+cut+"/"),
			GetDirectory(infile,"/stage0/full/near/MC/MRCC/"+horn+"/"+reco+"/"+run+"/"+cut+"/"));		
	}	
	
	}
	
}

TCanvas * Stage2::DrawParsCompareNear(string localPlotDirectory, string newpath, TDirectory *data, TDirectory *mc, string var, double xlow, double xhigh, int rebin, double lx, double ly, char * tname, string desc)
{

	

		TCanvas *c1 = new TCanvas("c","",500,600);
		c1->cd();
		
		TPad *ptop = new TPad("pt","",0,0.3,1,1);
		TPad *pbottom = new TPad("pb","",0,0,1,0.3);
		
		ptop->Draw();
		pbottom->Draw();
		
		ptop->cd();
		
		char tmp[200];
		sprintf(tmp,"%s",var.c_str());
		TH1D * p = (TH1D*)mc->FindObjectAny(tmp);
		if(p)p=(TH1D*)p->Clone();

		sprintf(tmp,"%s_nc",var.c_str());
		TH1D * p_nc = (TH1D*)mc->FindObjectAny(tmp);
		if(p_nc)p_nc=(TH1D*)p_nc->Clone();

		sprintf(tmp,"%s_cc",var.c_str());
		TH1D * p_cc = (TH1D*)mc->FindObjectAny(tmp);
		if(p_cc)p_cc=(TH1D*)p_cc->Clone();
	
		sprintf(tmp,"%s_beam",var.c_str());
		TH1D * p_beam = (TH1D*)mc->FindObjectAny(tmp);
		if(p_beam)p_beam=(TH1D*)p_beam->Clone();
								
		
		sprintf(tmp,"%s",var.c_str());
		TH1D * p_data = data ? (TH1D*)data->FindObjectAny(tmp) : 0;
		if(p_data)p_data=(TH1D*)p_data->Clone();								
//		else p_data=(TH1D*)p->Clone();
		printf("BAD\n");
		
		if(!p_data)
		{
			//printf("missing data\n");
//			return c1;
		}
		
		if(!p)
		{
			//printf("missing mc\n");
			return c1;
		}
		
			
		
		p->SetDirectory(0);
		p_nc->SetDirectory(0);
		p_cc->SetDirectory(0);
		p_beam->SetDirectory(0);
		if(p_data)p_data->SetDirectory(0);
		
		

		p->Rebin(rebin);
		p_nc->Rebin(rebin);
		p_cc->Rebin(rebin);
		p_beam->Rebin(rebin);
		if(p_data)p_data->Rebin(rebin);
		
		p->GetXaxis()->SetRangeUser(xlow,xhigh);
		if(p_data)p_data->GetXaxis()->SetRangeUser(xlow,xhigh);

		p_nc->SetLineColor(GetColor(NC));
		p_cc->SetLineColor(GetColor(CC));
		p_beam->SetLineColor(GetColor(Beam));
		if(p_data)p_data->SetLineColor(GetColor(Data));
		
		p->SetTitle("");
		if(p_data)p_data->SetTitle("");	

		p->SetXTitle(tname);
		p->GetXaxis()->SetTitleSize(0.03);
		p->GetXaxis()->SetTitleOffset(1.3);

		TH1D * p_comp =  p_data ? (TH1D*)p_data->Clone() : 0;

		if(p_data)p_data->SetXTitle(tname);
		if(p_data)p_data->GetXaxis()->SetTitleSize(0.03);
		if(p_data)p_data->GetXaxis()->SetTitleOffset(1.3);
		
		if(!p_data || ( p->GetMaximum() > p_data->GetMaximum() ) )
		{
			p->Draw("h");
			if(p_data)p_data->Draw("ssameh");
		}else{
			p_data->Draw("h");
			p->Draw("ssameh");
		}
		
		p_nc->Draw("ssameh");
		p_cc->Draw("ssameh");
		p_beam->Draw("ssameh");

		
	
		TLegend *l=new TLegend(lx,ly,lx+0.18,ly+0.18,"","NDC");
		l->AddEntry(p,"All");
		if(p_data)l->AddEntry(p_data,"Data");
		l->AddEntry(p_nc,GetName(NC).c_str());
		l->AddEntry(p_cc,GetName(CC).c_str());
		l->AddEntry(p_beam,GetName(Beam).c_str());
		l->Draw("Same");
		
		
		pbottom->cd();
		

		
		
		if(p_comp)
		{
		p_comp->SetDirectory(0);
		p_comp->SetName("");
		p_comp->SetTitle("");
		
		
//////////////		
		p_comp->Divide(p);
		p_comp->SetYTitle("Data/MC");
		p_comp->GetYaxis()->SetRangeUser(0.8,1.2);


//////////////
	//	p_comp->SetYTitle("(Data-MC)/Data");
	//	p_comp->GetYaxis()->SetRangeUser(-0.2,1.2);

		p_comp->GetXaxis()->SetRangeUser(xlow,xhigh);
		p->GetYaxis()->SetTitleSize(0.1);
		p->GetYaxis()->SetTitleOffset(0.4);
		p->GetYaxis()->CenterTitle(true);

		p_comp->Draw();		
		}		
		
		c1->cd(0);

		
		//save the plot to the save directory...
		TDirectory *sd = GetDirectory(outfile,newpath);
		
		c1->SetTitle(desc.c_str());

		sd->cd();

		c1->Clone(var.c_str())->Write(var.c_str());
		sd->Write();
		
		//output the plot to an image file...
	
		if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/eps/"+var+".eps").c_str());
		if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/gif/"+var+".gif").c_str());

		return c1;

}


void Stage2::DrawParsCompareNear(TDirectory *data, TDirectory *mc)
{

	if(!data && !mc)return;
	
	//save the plot to the save directory...
	vector<string>paths;
	Tokenize(mc->GetPath(),paths,"/");
	
	string newpath = "/stage2/plots/near/compare/";
	string localPlotDirectory=plotDirectory+"/near/compare/";
	for(int j=4;j<(int)paths.size();j++)
	{
		newpath+=paths[j]+"/";
		localPlotDirectory+=paths[j]+"/";
	}
	GetDirectory(outfile,newpath,2);
	
	
	
	if(plotDirectory!="")MakeTrueDirectory(localPlotDirectory+"/gif");		
	if(plotDirectory!="")MakeTrueDirectory(localPlotDirectory+"/eps");

	gStyle->SetOptStat(0);
	
	
	
	
	char tmp[200];
	for(int i=0;i<14;i++)
	{
	


		
		int rebin=1;
		double xlow=-10000;
		double xhigh=10000;
		
		//legend position
		double lx=0.7;
		double ly=0.7;
		
		switch(i)
		{
			case 0:
				xlow=0;
				xhigh=1.5;
				rebin=10;
				break;
			case 1:
				rebin=10;
				break;
			case 2:
				rebin=4;
				lx=0.2;
				break;
			case 3:
				xhigh=15;
				break;
			case 4:
				rebin=5;
				break;
			case 5:
				rebin=4;
				lx=0.2;
				break;
			case 6:
				rebin=4;
				xhigh=1.449;
//				rebin=4;
				break;
			case 7:
				rebin=4;
				break;
			case 8:
				xhigh=200;
				rebin=4;
				break;
			case 9:
				rebin=4;
				break;
			case 10:
				rebin=4;
				break;
			case 11:
				xlow=1;
				xhigh=5.5;
				rebin=4;
				break;
			case 12:
				rebin=2;
				break;	
			case 13:
				rebin=4;
				break;			
		};
		
		sprintf(tmp,"pars_%d",i);
		char tname[200];
		sprintf(tname,"Parameter %d - %s",i,GetParameterName(i).c_str());
		string desc = GetParameterName(i).c_str();
		
		TCanvas *c1 = DrawParsCompareNear(localPlotDirectory, newpath, data, mc, string(tmp), xlow, xhigh,  rebin,  lx,  ly,  tname, desc);
		c1->Close();
	}
	
	
	DrawParsCompareNear(localPlotDirectory, newpath, data, mc, "recoE", 0, 20,  1,  0.7,  0.7,  "Reconstructed Energy (Gev)", "Reconstructed Energy (Gev)")->Close();	
		
}


/////
//////////////////////


//////////////////////
///// SuperFOMs
void Stage2::DrawSuperFOMs(double sys)
{

	for(int i=0;i<4;i++)
	{
		string run="";	
		switch(i)
		{
			case 0:run="Run1";break;
			case 1:run="Run2";break;
			case 2:run="Run3";break;
			case 3:run="All";break;
		}


		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/fiducial"),"pidA",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/fiducial"),"pidB",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/fiducial"),"pidC",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/fiducial"),"pidD",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/fiducial"),"pidE",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/fiducial"),"pidF",sys);
	
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/presel"),"pidA",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/presel"),"pidB",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/presel"),"pidC",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/presel"),"pidD",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/presel"),"pidE",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/presel"),"pidF",sys);
	
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidA"),"pidA",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidB"),"pidB",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidC"),"pidC",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidD"),"pidD",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidE"),"pidE",sys);
		DrawSuperFOMs(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidF"),"pidF",sys);	

	}
}

void Stage2::DrawSuperFOMs(TDirectory *d, string pid, double sys)
{

	if(!d)return;
	TH1D * sig = (TH1D*)d->FindObjectAny((pid+"_sig").c_str())->Clone();
	TH1D * bg = (TH1D*)d->FindObjectAny((pid).c_str())->Clone();
	
	sig->SetDirectory(0);
	bg->SetDirectory(0);
	
	if(!sig || !bg)return;
	
	bg->Add(sig,-1);
	
	int bins = sig->GetNbinsX();
	TH1D * sf = new TH1D("sf","sf",bins,sig->GetXaxis()->GetXmin(),sig->GetXaxis()->GetXmax());
	sf->SetDirectory(0);
	
	for(int i=0;i<bins;i++)
	{
		double b=bg->Integral(i+1,bins+1);
		double s=sig->Integral(i+1,bins+1);
		double t=0;
		if(b)t=s/sqrt(b+b*b*sys*sys);
		sf->SetBinContent(i+1,t);
	}
	
	TCanvas *c1 = new TCanvas();
	char tmp[200];
	sprintf(tmp,"SuperFOM %s (%.2f%% sys)",pid.c_str(),sys*100.);
	sf->SetTitle(tmp);
	sf->SetXTitle("PID Cut");
	sf->Draw();
	
	
	//save the plot to the save directory...
	vector<string>paths;
	Tokenize(d->GetPath(),paths,"/");
	
	string newpath = "/stage2/plots/";
	string localPlotDirectory=plotDirectory+"/";	
	for(int j=2;j<(int)paths.size();j++)
	{
		newpath+=paths[j]+"/";
		localPlotDirectory+=paths[j]+"/";
	}
	GetDirectory(outfile,newpath,1);
	
	//save the plot to the save directory...
	TDirectory *sd = GetDirectory(outfile,newpath);
		
	c1->SetTitle(((string)"SuperFOM "+pid).c_str());

	sd->cd();
	c1->Clone(((string)"SuperFOM_"+pid).c_str())->Write(((string)"SuperFOM_"+pid).c_str());
	sd->Write();
		
	//output the plot to an image file...
	if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/eps/SuperFOM_"+pid+".eps").c_str());
	if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/gif/SuperFOM_"+pid+".gif").c_str());
	
}

///// SuperFOMs
//////////////////////



//////////////////////
///// DrawParsFar


TCanvas * Stage2::DrawParsFar(double sigscale, string localPlotDirectory, string newpath, TDirectory * /*data*/, TDirectory *mc, string var, double xlow, double xhigh, int rebin, double lx, double ly, char * tname, string desc)
{


		TCanvas *c1 = new TCanvas("c","",500,600);
		c1->cd();
		
		char tmp[200];
		sprintf(tmp,"%s",var.c_str());
		TH1D * p = (TH1D*)mc->FindObjectAny(tmp)->Clone();

		sprintf(tmp,"%s_nc",var.c_str());
		TH1D * p_nc = (TH1D*)mc->FindObjectAny(tmp)->Clone();

		sprintf(tmp,"%s_cc",var.c_str());
		TH1D * p_cc = (TH1D*)mc->FindObjectAny(tmp)->Clone();
	
		sprintf(tmp,"%s_beam",var.c_str());
		TH1D * p_beam = (TH1D*)mc->FindObjectAny(tmp)->Clone();
		
		sprintf(tmp,"%s_tau",var.c_str());
		TH1D * p_tau = (TH1D*)mc->FindObjectAny(tmp)->Clone();
		
		sprintf(tmp,"%s_sig",var.c_str());
		TH1D * p_sig = (TH1D*)mc->FindObjectAny(tmp)->Clone();								
		

		
		p->SetDirectory(0);
		p_sig->SetDirectory(0);
		p_nc->SetDirectory(0);
		p_cc->SetDirectory(0);
		p_beam->SetDirectory(0);
		p_tau->SetDirectory(0);
		
		
		p->Rebin(rebin);
		p_nc->Rebin(rebin);
		p_cc->Rebin(rebin);
		p_beam->Rebin(rebin);
		p_tau->Rebin(rebin);
		p_sig->Rebin(rebin);
		
		p->GetXaxis()->SetRangeUser(xlow,xhigh);
		p_sig->GetXaxis()->SetRangeUser(xlow,xhigh);

		p_nc->SetLineColor(GetColor(NC));
		p_cc->SetLineColor(GetColor(CC));
		p_beam->SetLineColor(GetColor(Beam));
		p_tau->SetLineColor(GetColor(Tau));
		p_sig->SetLineColor(GetColor(Sig));
		
		p->SetTitle("");
		p_sig->SetTitle("");		


		p->SetXTitle(tname);
		p->GetXaxis()->SetTitleSize(0.03);
		p->GetXaxis()->SetTitleOffset(1.3);

		p_sig->SetXTitle(tname);
		p_sig->GetXaxis()->SetTitleSize(0.03);
		p_sig->GetXaxis()->SetTitleOffset(1.3);

		
		p->SetYTitle("Events / 7e20 POT");
			
		p_sig->Scale(sigscale);
		
		if(p->GetMaximum() > p_sig->GetMaximum())
		{
			p->Draw("h");
			p_sig->Draw("ssameh");
		}else{
			p_sig->Draw("h");
			p->Draw("ssameh");
		}
				
		p_nc->Draw("ssameh");
		p_cc->Draw("ssameh");
		p_beam->Draw("ssameh");
		p_tau->Draw("ssameh");

	
	
		TLegend *l=new TLegend(lx,ly,lx+0.18,ly+0.18,"","NDC");
		l->AddEntry(p,"All");
		char zzz[200];
		if(sigscale!=1)
		{
			sprintf(zzz,"%sx%.1f",GetName(Sig).c_str(),sigscale);
			l->AddEntry(p_sig,zzz);
		}else{
			l->AddEntry(p_sig,GetName(Sig).c_str());
		}		
		l->AddEntry(p_nc,GetName(NC).c_str());
		l->AddEntry(p_cc,GetName(CC).c_str());
		l->AddEntry(p_beam,GetName(Beam).c_str());
		l->AddEntry(p_tau,GetName(Tau).c_str());
		l->Draw("Same");
		
		

		//save the plot to the save directory...
		TDirectory *sd = GetDirectory(outfile,newpath);
		
		c1->SetTitle(desc.c_str());

		sd->cd();
		c1->Clone(var.c_str())->Write(var.c_str());
		sd->Write();
		
		//output the plot to an image file...
		if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/eps/"+var+".eps").c_str());
		if(plotDirectory!="")c1->SaveAs((localPlotDirectory+"/gif/"+var+".gif").c_str());	


	return c1;

}


void Stage2::DrawParsFar()
{
	
	for(int i=0;i<4;i++)
	{
		string run="";	
		switch(i)
		{
			case 0:run="Run1";break;
			case 1:run="Run2";break;
			case 2:run="Run3";break;
			case 3:run="All";break;
		}

		DrawParsFar(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/fiducial"),7);
		DrawParsFar(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/presel"),10);
		DrawParsFar(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidA"));	
		DrawParsFar(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidB"));	
		DrawParsFar(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidC"));	
		DrawParsFar(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidD"));	
		DrawParsFar(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidE"));	
		DrawParsFar(GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/pidF"));	

	}
}


void Stage2::DrawParsFar(TDirectory *d, double sigscale)
{
	if(!d)return;
	//save the plot to the save directory...
	vector<string>paths;
	Tokenize(d->GetPath(),paths,"/");
	
	string newpath = "/stage2/plots/";
	string localPlotDirectory=plotDirectory+"/";	
	for(int j=2;j<(int)paths.size();j++)
	{
		newpath+=paths[j]+"/";
		localPlotDirectory+=paths[j]+"/";
	}
	GetDirectory(outfile,newpath,2);
		

	if(plotDirectory!="")MakeTrueDirectory(localPlotDirectory+"/gif");			
	if(plotDirectory!="")MakeTrueDirectory(localPlotDirectory+"/eps");	
	
	
	gStyle->SetOptStat(0);
	
	
	char tmp[200];
	for(int i=0;i<14;i++)
	{
	

		
		
		int rebin=1;
		double xlow=-10000;
		double xhigh=10000;
		
		//legend position
		double lx=0.7;
		double ly=0.7;
		
		switch(i)
		{
			case 0:
				xlow=0;
				xhigh=1.5;
				rebin=10;
				break;
			case 1:
				rebin=10;
				break;
			case 2:
				rebin=2;
				lx=0.2;
				break;
			case 3:
				xhigh=15;
				break;
			case 4:
				rebin=5;
				break;
			case 5:
				rebin=2;
				lx=0.2;
				break;
			case 6:
				xhigh=1.449;
//				rebin=4;
				break;
			case 7:
				break;
			case 11:
				xlow=1;
				xhigh=5.5;
				break;
			case 12:
				rebin=2;
				break;	
			case 13:
				break;			
		};
		
	
		sprintf(tmp,"pars_%d",i);
		char tname[200];
		sprintf(tname,"Parameter %d - %s",i,GetParameterName(i).c_str());
		string desc = GetParameterName(i).c_str();
		
		TCanvas *c1 = DrawParsFar(sigscale, localPlotDirectory, newpath, 0, d, string(tmp), xlow, xhigh,  rebin,  lx,  ly,  tname, desc);
		c1->Close();
		
	}
	
	
	DrawParsFar(sigscale, localPlotDirectory, newpath, 0, d, "recoE", 0, 20,  1,  0.7,  0.7,  "Reconstructed Energy (Gev)", "Reconstructed Energy (Gev)")->Close();	


	DrawParsFar(sigscale, localPlotDirectory, newpath, 0, d, "pidA", -0.4, 1.4,  2,  0.7,  0.7,  "pidA", "pidA")->Close();	
	DrawParsFar(sigscale, localPlotDirectory, newpath, 0, d, "pidB", -0.4, 1.4,  2,  0.7,  0.7,  "pidB", "pidB")->Close();	
	DrawParsFar(sigscale, localPlotDirectory, newpath, 0, d, "pidC", -0.4, 1.4,  2,  0.7,  0.7,  "pidC", "pidC")->Close();	
	DrawParsFar(sigscale, localPlotDirectory, newpath, 0, d, "pidD", -0.4, 1.4,  2,  0.7,  0.7,  "pidD", "pidD")->Close();	
	DrawParsFar(sigscale, localPlotDirectory, newpath, 0, d, "pidE", -0.4, 1.4,  2,  0.7,  0.7,  "pidE", "pidE")->Close();	
	DrawParsFar(sigscale, localPlotDirectory, newpath, 0, d, "pidF", -0.4, 1.4,  2,  0.7,  0.7,  "pidF", "pidF")->Close();	
	
	
}



/////
//////////////////////



//////////////////////
/////


void Stage2::CountComponentsByCutLevelFar(double sys)
{
		
	for(int i=0;i<4;i++)
	{
		string run="";
		switch(i)
		{
			case 0:run="Run1";break;
			case 1:run="Run2";break;
			case 2:run="Run3";break;
			case 3:run="All";break;
		}


		CountComponentsByCutLevelFar(sys,GetDirectory(infile,"/stage0/full/far/MC/standard/ParticlePID/"+run+"/"));
	
		CountComponentsByCutLevelFar(sys,GetDirectory(infile,"/stage0/full/far/MC/standard/normal/"+run+"/"));
	}
}

void Stage2::CountComponentsByCutLevelFar(double sys, TDirectory * base)
{
	if(!base)return;
	//for each stage0 far directory
	
	ostringstream raw;
	ostringstream latex;
	
	latex << "\\begin{tabular}{ l || r || r | r || r | r | r | r | } \n\\hline\n";
	
	raw << setw(20)<<"Cut"<<setw(10)<<"SFOM"<<setw(10)<<"Sig"<<setw(10)<<"BG"<<setw(10)<<"NC"<<setw(10)<<"CC"<<setw(10)<<"Tau"<<setw(10)<<"Beam\n";
	
	latex << "Cut&SFOM&Sig&BG&NC&CC&Tau&Beam\\\\\n\\hline\n";

	CountComponentsByCutLevelFar(sys, raw, latex, GetDirectory(base,"fiducial"));
	CountComponentsByCutLevelFar(sys, raw, latex, GetDirectory(base,"presel"));	
	CountComponentsByCutLevelFar(sys, raw, latex, GetDirectory(base,"pidA"));
	CountComponentsByCutLevelFar(sys, raw, latex, GetDirectory(base,"pidB"));
	CountComponentsByCutLevelFar(sys, raw, latex, GetDirectory(base,"pidC"));
	CountComponentsByCutLevelFar(sys, raw, latex, GetDirectory(base,"pidD"));
	CountComponentsByCutLevelFar(sys, raw, latex, GetDirectory(base,"pidE"));
	CountComponentsByCutLevelFar(sys, raw, latex, GetDirectory(base,"pidF"));
	CountComponentsByCutLevelFar(sys, raw, latex, GetDirectory(base,"ann11"));
	CountComponentsByCutLevelFar(sys, raw, latex, GetDirectory(base,"ann11_firebird"));
	

	latex <<"\\end{tabular}\n";
	
	
	cout<<base->GetPath()<<"\n";
	
	cout<<raw.str()<<"\n";
	//cout<<latex.str()<<"\n";
	
	//TLatex *ll=new TLatex(0.1,0.9,"#nu");//latex.str().c_str());
	//ll->Draw();

}


void Stage2::CountComponentsByCutLevelFar(double sys, ostringstream &sraw, ostringstream &slatex, TDirectory *d)
{
	if(!d)return;
	
	TH1D * recoE_sig  = (TH1D*)d->FindObjectAny("recoE_sig");
	TH1D * recoE_nc  = (TH1D*)d->FindObjectAny("recoE_nc");
	TH1D * recoE_cc  = (TH1D*)d->FindObjectAny("recoE_cc");
	TH1D * recoE_beam  = (TH1D*)d->FindObjectAny("recoE_beam");
	TH1D * recoE_tau  = (TH1D*)d->FindObjectAny("recoE_tau");
	

	double t[5];
	for(int i=0;i<5;i++)t[i]=0;
	
	t[0]=recoE_nc->Integral();
	t[1]=recoE_cc->Integral();
	t[2]=recoE_sig->Integral();
	t[3]=recoE_tau->Integral();
	t[4]=recoE_beam->Integral();
	
	
	
	
	double bg=t[0]+t[1]+t[3]+t[4];
	double sf = t[2]/sqrt(bg+bg*bg*sys*sys);
		
	string dirname="";
	string aname = d->GetName();
	if(aname=="fiducial")dirname="Fiducial";	
	else if(aname=="presel")dirname="Preselection";
	else if(aname=="pidA")dirname="Pid A";
	else if(aname=="pidB")dirname="Pid B";
	else if(aname=="pidC")dirname="Pid C";
	else if(aname=="pidD")dirname="Pid D";
	else if(aname=="pidE")dirname="Pid E";
	else if(aname=="pidF")dirname="Pid F";
	else if(aname=="ann11")dirname="ANN 11";
	else if(aname=="ann11_firebird")dirname="ANN 11 (Firebird)";
	else dirname="?";
	
		
	sraw << fixed << setprecision(2)<<setw(20)<<dirname<<setw(10)<< sf<<setw(10)<<t[2]<<setw(10)<<bg<<setw(10)<<t[0]<<setw(10)<<t[1]<<setw(10)<<t[3]<<setw(10)<<t[4]<<"\n";
	slatex << dirname<<"&"<< sf<<" | "<<t[2]<<"&"<<bg<<"&"<<t[0]<<"&"<<t[1]<<"&"<<t[3]<<"&"<<t[4]<<"\\\\\n";
	
	//reset the format
	std::cout.setf( std::ios::fmtflags(), std::ios::floatfield ) ;
	
}

/////
//////////////////////




//////////////////////
/////


void Stage2::CountComponentsByCutLevelNear(double sys)
{
		
	for(int i=0;i<4;i++)
	{
		string run="";
		switch(i)
		{
			case 0:run="Run1";break;
			case 1:run="Run2";break;
			case 2:run="Run3";break;
			case 3:run="All";break;
		}

		
		CountComponentsByCutLevelNear(sys,GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/ParticlePID/"+run+"/"),
			GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/ParticlePID/"+run+"/"));
	
		CountComponentsByCutLevelNear(sys,GetDirectory(infile,"/stage0/full/near/data/standard/horn_on/normal/"+run+"/"),
			GetDirectory(infile,"/stage0/full/near/MC/standard/horn_on/normal/"+run+"/"));
	
		CountComponentsByCutLevelNear(sys,GetDirectory(infile,"/stage0/full/near/data/standard/horn_off/ParticlePID/"+run+"/"),
			GetDirectory(infile,"/stage0/full/near/MC/standard/horn_off/ParticlePID/"+run+"/"));
	
		CountComponentsByCutLevelNear(sys,GetDirectory(infile,"/stage0/full/near/data/standard/horn_off/normal/"+run+"/"),
			GetDirectory(infile,"/stage0/full/near/MC/standard/horn_off/normal/"+run+"/"));
	
		CountComponentsByCutLevelNear(sys,GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/ParticlePID/"+run+"/"),
			GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/ParticlePID/"+run+"/"));
	
		CountComponentsByCutLevelNear(sys,GetDirectory(infile,"/stage0/full/near/data/MRCC/horn_on/normal/"+run+"/"),
			GetDirectory(infile,"/stage0/full/near/MC/MRCC/horn_on/normal/"+run+"/"));		
	}
}

void Stage2::CountComponentsByCutLevelNear(double sys, TDirectory * data, TDirectory * mc)
{
	if(!data || !mc)
	{
		printf("error loading data %d mc %d\n",data?1:0,mc?1:0);
		return;
	}
	//for each stage0 far directory
	
	ostringstream raw;
	ostringstream latex;

	cout<<data->GetPath()<<"\n";
	cout<<mc->GetPath()<<"\n";
	
	//reset buffers
	raw.str("");
	latex.str("");
	
	latex << "\\begin{tabular}{ l ||  r | r | r || r | r  | } \n\\hline\n";
	
	raw << setw(20)<<"Cut"<<setw(10)<<"Data"<<setw(10)<<"BG"<<setw(10)<<"NC"<<setw(10)<<"CC"<<setw(10)<<"Beam\n";
	
	latex << "Cut&Data&BG&NC&CC&Beam\\\\\n\\hline\n";

	CountComponentsByCutLevelNear(sys, raw, latex, GetDirectory(data,"fiducial"), GetDirectory(mc,"fiducial"));
	CountComponentsByCutLevelNear(sys, raw, latex, GetDirectory(data,"presel"), GetDirectory(mc,"presel"));	
	CountComponentsByCutLevelNear(sys, raw, latex, GetDirectory(data,"pidA"), GetDirectory(mc,"pidA"));
	CountComponentsByCutLevelNear(sys, raw, latex, GetDirectory(data,"pidB"), GetDirectory(mc,"pidB"));
	CountComponentsByCutLevelNear(sys, raw, latex, GetDirectory(data,"pidC"), GetDirectory(mc,"pidC"));
	CountComponentsByCutLevelNear(sys, raw, latex, GetDirectory(data,"pidD"), GetDirectory(mc,"pidD"));
	CountComponentsByCutLevelNear(sys, raw, latex, GetDirectory(data,"pidE"), GetDirectory(mc,"pidE"));
	CountComponentsByCutLevelNear(sys, raw, latex, GetDirectory(data,"pidF"), GetDirectory(mc,"pidF"));
	CountComponentsByCutLevelNear(sys, raw, latex, GetDirectory(data,"ann11"), GetDirectory(mc,"ann11"));
	CountComponentsByCutLevelNear(sys, raw, latex, GetDirectory(data,"ann11_firebird"),
		 GetDirectory(mc,"ann11_firebird"));
	

	latex <<"\\end{tabular}\n";
	
	

	
	cout<<raw.str()<<"\n";
	//cout<<latex.str()<<"\n";
	
	//TLatex *ll=new TLatex(0.1,0.9,"#nu");//latex.str().c_str());
	//ll->Draw();

}


void Stage2::CountComponentsByCutLevelNear(double /*sys*/, ostringstream &sraw, ostringstream &slatex, TDirectory *data, TDirectory *mc)
{
	if(!data || !mc)return;
	
//	TH1D * recoE_sig  = (TH1D*)mc->FindObjectAny("recoE_sig");
	TH1D * recoE_nc  = (TH1D*)mc->FindObjectAny("recoE_nc");
	TH1D * recoE_cc  = (TH1D*)mc->FindObjectAny("recoE_cc");
	TH1D * recoE_beam  = (TH1D*)mc->FindObjectAny("recoE_beam");
//	TH1D * recoE_tau  = (TH1D*)mc->FindObjectAny("recoE_tau");
	
	TH1D * recoE_data  = (TH1D*)data->FindObjectAny("recoE");
	
	double t[5];
	for(int i=0;i<5;i++)t[i]=0;
	
	t[0]=recoE_nc->Integral();
	t[1]=recoE_cc->Integral();
//	t[2]=recoE_sig->Integral();
//	t[3]=recoE_tau->Integral();
	t[4]=recoE_beam->Integral();
	
	double dd=recoE_data->Integral();
	
	
	double bg=t[0]+t[1]+t[3]+t[4];
//	double sf = t[2]/sqrt(bg+bg*bg*sys*sys);
		
	string dirname="";
	string aname = data->GetName();
	if(aname=="fiducial")dirname="Fiducial";	
	else if(aname=="presel")dirname="Preselection";
	else if(aname=="pidA")dirname="Pid A";
	else if(aname=="pidB")dirname="Pid B";
	else if(aname=="pidC")dirname="Pid C";
	else if(aname=="pidD")dirname="Pid D";
	else if(aname=="pidE")dirname="Pid E";
	else if(aname=="pidF")dirname="Pid F";
	else if(aname=="ann11")dirname="ANN 11";
	else if(aname=="ann11_firebird")dirname="ANN 11 (Firebird)";
	else dirname="?";
	
		
	sraw << fixed << setprecision(2)<<setw(20)<<dirname<< setw(10)<<dd<<setw(10)<<bg<<setw(10)<<t[0]<<setw(10)<<t[1]<<setw(10)<<t[4]<<"\n";
	slatex << dirname<<" | "<<dd<<"&"<<bg<<"&"<<t[0]<<"&"<<t[1]<<"&"<<t[4]<<"\\\\\n";
	
	//reset the format
	std::cout.setf( std::ios::fmtflags(), std::ios::floatfield ) ;
	
}

/////
/////////////////////








void Stage2::Run()
{

	PrepareFiles();
	

	

	//SaveAndClose();
}

void Stage2::SaveAndClose()
{

	outfile->Write();
	outfile->Close();

}

string Stage2::GetName(int code)
{
	switch(code)
	{
		case 0: return "NC";
		case 1: return "CC";
		case 2: return "#nu_{e}";
		case 3: return "#nu_{#tau}";
		case 4: return "Beam #nu_{e}";
		case Data: return "Data";
	}
	return "?";
}



int Stage2::GetColor(int code)
{
	switch(code)
	{
		case 0: return 3;
		case 1: return 4;
		case 2: return 2;
		case 3: return 5;
		case 4: return 8;
		case Data: return 2;  //will need a different color in the far...
	}
	return 0;
}



string Stage2::GetParameterName(int code)
{
	switch(code)
	{
		case 0: return "Longest Particle Path Length (m)";
		case 1: return "Molier Radius R (m)";
		case 2: return "Reconstructed EM Shower Fraction";
		case 3: return "Total number of Reconstructed Particles";
		case 4: return "Weighted Angle Off of Z";
		case 5: return "Energy of Largest Reconstructed Particle / Total Reconstructed Energy";
		case 6: return "EM Shower Fit Par b of Largest Reconstructed Particle";
		case 7: return "EM Shower Fit Par e0 of Largest Reconstructed Particle (MEU)";
		case 8: return "EM Shower Fit #Chi^{2} of Largest Reconstructed Particle";		
		case 9: return "Largest Particle Peak Difference (MEU)";
		case 10: return "Largest Particle EM Shower Profile Comparison #Chi^{2}/NDF";
		case 11: return "EM Shower Fit Par a of Largest Reconstructed Particle";
		case 12: return "Number of Reconstructed Clusters";
		case 13: return "EM Shower Fit of Largest Reconstructed Particle Parameters a/e0";
	}
	return "?";
}


