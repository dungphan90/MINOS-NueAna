#include <iostream>
#include "TFile.h"
#include "TChain.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TMultiGraph.h"
#include "TLine.h"
#include "TText.h"
#include "TF1.h"
#include "Conventions/DetectorType.h"
#include "NueAna/Reweight/NueRW.h" 
#include "NueAna/macros/chooz.C"

using namespace std;

void PlotESpect(const char* ndf, const char* fdf);
void FitNue(const char *challengefile, const char *rwfile, float ndscale, float fdscale, bool save=false);
void Sensitivity(const char *rwfile,Double_t ss2th=-1);
void StudyND(const char *rwfile,Bool_t MDC=false);
void ComputeChi2(TH1 *h1, TH1 *h2, double &chi2, int &ndof,int method=0);
void FindUpperBounds(TGraph *g,float &b68,float &b90,float &b99);
void FindLowerBounds(TGraph *g,float &b68,float &b90,float &b99);


double Chi2Prob(double *x, double *par)
{
   return TMath::Exp(-1.*x[0]/2.)*TMath::Power(x[0],par[0]/2.-1)/(TMath::Power(2,par[0]/2.)*TMath::Gamma(par[0]/2.));
}

void PlotESpect(const char* ndf, const char* fdf)
{

  TH1F *tnomnd=new TH1F("tnomnd","Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *munomnd=new TH1F("munomnd","#{mu} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ebnomnd=new TH1F("ebnomnd","e_{beam} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *taunomnd=new TH1F("taunomnd","#{tau} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ncnomnd=new TH1F("ncnomnd","NC Spectrum ND;GeV;# Events/GeV",20,0,20);

  TH1F *tnomfd=new TH1F("tnomfd","Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *munomfd=new TH1F("munomfd","#{mu} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ebnomfd=new TH1F("ebnomfd","e_{beam} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *taunomfd=new TH1F("taunomfd","#{tau} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ncnomfd=new TH1F("ncnomfd","NC Spectrum FD;GeV;# Events/GeV",20,0,20);

  TH1F *tqelupnd=new TH1F("tqelupnd","Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *muqelupnd=new TH1F("muqelupnd","#{mu} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ebqelupnd=new TH1F("ebqelupnd","e_{beam} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *tauqelupnd=new TH1F("tauqelupnd","#{tau} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ncqelupnd=new TH1F("ncqelupnd","NC Spectrum ND;GeV;# Events/GeV",20,0,20);

  TH1F *tqeldownnd=new TH1F("tqeldownnd","Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *muqeldownnd=new TH1F("muqeldownnd","#{mu} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ebqeldownnd=new TH1F("ebqeldownnd","e_{beam} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *tauqeldownnd=new TH1F("tauqeldownnd","#{tau} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ncqeldownnd=new TH1F("ncqeldownnd","NC Spectrum ND;GeV;# Events/GeV",20,0,20);


  TH1F *tqelupfd=new TH1F("tqelupfd","Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *muqelupfd=new TH1F("muqelupfd","#{mu} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ebqelupfd=new TH1F("ebqelupfd","e_{beam} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *tauqelupfd=new TH1F("tauqelupfd","#{tau} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ncqelupfd=new TH1F("ncqelupfd","NC Spectrum FD;GeV;# Events/GeV",20,0,20);

  TH1F *tqeldownfd=new TH1F("tqeldownfd","Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *muqeldownfd=new TH1F("muqeldownfd","#{mu} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ebqeldownfd=new TH1F("ebqeldownfd","e_{beam} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *tauqeldownfd=new TH1F("tauqeldownfd","#{tau} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ncqeldownfd=new TH1F("ncqeldownfd","NC Spectrum FD;GeV;# Events/GeV",20,0,20);

  TH1F *tresupnd=new TH1F("tresupnd","Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *muresupnd=new TH1F("muresupnd","#{mu} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ebresupnd=new TH1F("ebresupnd","e_{beam} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *tauresupnd=new TH1F("tauresupnd","#{tau} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ncresupnd=new TH1F("ncresupnd","NC Spectrum ND;GeV;# Events/GeV",20,0,20);

  TH1F *tresdownnd=new TH1F("tresdownnd","Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *muresdownnd=new TH1F("muresdownnd","#{mu} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ebresdownnd=new TH1F("ebresdownnd","e_{beam} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *tauresdownnd=new TH1F("tauresdownnd","#{tau} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ncresdownnd=new TH1F("ncresdownnd","NC Spectrum ND;GeV;# Events/GeV",20,0,20);


  TH1F *tresupfd=new TH1F("tresupfd","Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *muresupfd=new TH1F("muresupfd","#{mu} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ebresupfd=new TH1F("ebresupfd","e_{beam} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *tauresupfd=new TH1F("tauresupfd","#{tau} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ncresupfd=new TH1F("ncresupfd","NC Spectrum FD;GeV;# Events/GeV",20,0,20);

  TH1F *tresdownfd=new TH1F("tresdownfd","Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *muresdownfd=new TH1F("muresdownfd","#{mu} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ebresdownfd=new TH1F("ebresdownfd","e_{beam} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *tauresdownfd=new TH1F("tauresdownfd","#{tau} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ncresdownfd=new TH1F("ncresdownfd","NC Spectrum FD;GeV;# Events/GeV",20,0,20);


  TH1F *tknoupnd=new TH1F("tknoupnd","Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *muknoupnd=new TH1F("muknoupnd","#{mu} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ebknoupnd=new TH1F("ebknoupnd","e_{beam} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *tauknoupnd=new TH1F("tauknoupnd","#{tau} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ncknoupnd=new TH1F("ncknoupnd","NC Spectrum ND;GeV;# Events/GeV",20,0,20);

  TH1F *tknodownnd=new TH1F("tknodownnd","Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *muknodownnd=new TH1F("muknodownnd","#{mu} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ebknodownnd=new TH1F("ebknodownnd","e_{beam} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *tauknodownnd=new TH1F("tauknodownnd","#{tau} Spectrum ND;GeV;# Events/GeV",20,0,20);
  TH1F *ncknodownnd=new TH1F("ncknodownnd","NC Spectrum ND;GeV;# Events/GeV",20,0,20);

  TH1F *tknoupfd=new TH1F("tknoupfd","Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *muknoupfd=new TH1F("muknoupfd","#{mu} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ebknoupfd=new TH1F("ebknoupfd","e_{beam} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *tauknoupfd=new TH1F("tauknoupfd","#{tau} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ncknoupfd=new TH1F("ncknoupfd","NC Spectrum FD;GeV;# Events/GeV",20,0,20);

  TH1F *tknodownfd=new TH1F("tknodownfd","Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *muknodownfd=new TH1F("muknodownfd","#{mu} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ebknodownfd=new TH1F("ebknodownfd","e_{beam} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *tauknodownfd=new TH1F("tauknodownfd","#{tau} Spectrum FD;GeV;# Events/GeV",20,0,20);
  TH1F *ncknodownfd=new TH1F("ncknodownfd","NC Spectrum FD;GeV;# Events/GeV",20,0,20);

  TFile nd(ndf);
  TTree *ndtree=(TTree *)nd.Get("rwtree");

  NueRW *nrw = new NueRW();
  ndtree->SetBranchAddress("NueRW",&nrw);

  for(int z=0;z<34;z++){
    ndtree->GetEntry(z);
    if(z==0){
      for(int i=0;i<nrw->EBINS;i++){
	tnomnd->Fill(i,nrw->nbgE[i]);
	munomnd->Fill(i,nrw->nnumuE[i]);
	ebnomnd->Fill(i,nrw->nnuebE[i]);
	taunomnd->Fill(i,nrw->nnutauE[i]);
	ncnomnd->Fill(i,nrw->nncE[i]);
      }
    }

    if(nrw->res_ma==1&&nrw->kno_r112==1){
      if(nrw->qel_ma>1.47&&nrw->qel_ma<1.53){
	for(int i=0;i<nrw->EBINS;i++){
	  tqelupnd->Fill(i,nrw->nbgE[i]);
	  muqelupnd->Fill(i,nrw->nnumuE[i]);
	  ebqelupnd->Fill(i,nrw->nnuebE[i]);
	  tauqelupnd->Fill(i,nrw->nnutauE[i]);
	  ncqelupnd->Fill(i,nrw->nncE[i]);
	}
      }
      
      if(nrw->qel_ma>.47&&nrw->qel_ma<.53){
	for(int i=0;i<nrw->EBINS;i++){
	  tqeldownnd->Fill(i,nrw->nbgE[i]);
	  muqeldownnd->Fill(i,nrw->nnumuE[i]);
	  ebqeldownnd->Fill(i,nrw->nnuebE[i]);
	  tauqeldownnd->Fill(i,nrw->nnutauE[i]);
	  ncqeldownnd->Fill(i,nrw->nncE[i]);
	}
      }
    }

    if(nrw->qel_ma==1&&nrw->kno_r112==1){
      if(nrw->res_ma>1.07&&nrw->res_ma<1.13){
	for(int i=0;i<nrw->EBINS;i++){
	  tresupnd->Fill(i,nrw->nbgE[i]);
	  muresupnd->Fill(i,nrw->nnumuE[i]);
	  ebresupnd->Fill(i,nrw->nnuebE[i]);
	  tauresupnd->Fill(i,nrw->nnutauE[i]);
	  ncresupnd->Fill(i,nrw->nncE[i]);
	}
      }
      
      if(nrw->res_ma>.87&&nrw->res_ma<.93){
	for(int i=0;i<nrw->EBINS;i++){
	  tresdownnd->Fill(i,nrw->nbgE[i]);
	  muresdownnd->Fill(i,nrw->nnumuE[i]);
	  ebresdownnd->Fill(i,nrw->nnuebE[i]);
	  tauresdownnd->Fill(i,nrw->nnutauE[i]);
	  ncresdownnd->Fill(i,nrw->nncE[i]);
	}
      }
    }

    if(nrw->qel_ma==1&&nrw->res_ma==1){
      if(nrw->kno_r112>1.47&&nrw->kno_r112<1.53){
	for(int i=0;i<nrw->EBINS;i++){
	  tknoupnd->Fill(i,nrw->nbgE[i]);
	  muknoupnd->Fill(i,nrw->nnumuE[i]);
	  ebknoupnd->Fill(i,nrw->nnuebE[i]);
	  tauknoupnd->Fill(i,nrw->nnutauE[i]);
	  ncknoupnd->Fill(i,nrw->nncE[i]);
	}
      }
      //      cout<<nrw->kno_r112<<endl;
      if(nrw->kno_r112>.47&&nrw->kno_r112<.53){
	//	cout<<"filling knodown nd"<<endl;
	for(int i=0;i<nrw->EBINS;i++){
	  tknodownnd->Fill(i,nrw->nbgE[i]);
	  muknodownnd->Fill(i,nrw->nnumuE[i]);
	  ebknodownnd->Fill(i,nrw->nnuebE[i]);
	  tauknodownnd->Fill(i,nrw->nnutauE[i]);
	  ncknodownnd->Fill(i,nrw->nncE[i]);
	}
      }
    }

  }


  TFile fd(fdf);
  TTree *fdtree=(TTree *)fd.Get("rwtree");

  NueRW *frw = new NueRW();
  fdtree->SetBranchAddress("NueRW",&frw);

  for(int z=0;z<34;z++){
    fdtree->GetEntry(z);
    if(z==0){
      for(int i=0;i<frw->EBINS;i++){
	tnomfd->Fill(i,frw->nbgE[i]-frw->nnuebE[i]/2.);
	munomfd->Fill(i,frw->nnumuE[i]);
	ebnomfd->Fill(i,frw->nnuebE[i]/2.);
	taunomfd->Fill(i,frw->nnutauE[i]);
	ncnomfd->Fill(i,frw->nncE[i]);
      }
    }

    if(frw->res_ma==1&&frw->kno_r112==1){
      if(frw->qel_ma>1.47&&frw->qel_ma<1.53){
	for(int i=0;i<frw->EBINS;i++){
	  tqelupfd->Fill(i,frw->nbgE[i]-frw->nnuebE[i]/2.);
	  muqelupfd->Fill(i,frw->nnumuE[i]);
	  ebqelupfd->Fill(i,frw->nnuebE[i]/2.);
	  tauqelupfd->Fill(i,frw->nnutauE[i]);
	  ncqelupfd->Fill(i,frw->nncE[i]);
	}
      }
      
      if(frw->qel_ma>.47&&frw->qel_ma<.53){
	for(int i=0;i<nrw->EBINS;i++){
	  tqeldownfd->Fill(i,frw->nbgE[i]-frw->nnuebE[i]/2.);
	  muqeldownfd->Fill(i,frw->nnumuE[i]);
	  ebqeldownfd->Fill(i,frw->nnuebE[i]/2.);
	  tauqeldownfd->Fill(i,frw->nnutauE[i]);
	  ncqeldownfd->Fill(i,frw->nncE[i]);
	}
      }
    }

    if(frw->qel_ma==1&&frw->kno_r112==1){
      if(frw->res_ma>1.07&&frw->res_ma<1.13){
	for(int i=0;i<frw->EBINS;i++){
	  tresupfd->Fill(i,frw->nbgE[i]-frw->nnuebE[i]/2.);
	  muresupfd->Fill(i,frw->nnumuE[i]);
	  ebresupfd->Fill(i,frw->nnuebE[i]/2.);
	  tauresupfd->Fill(i,frw->nnutauE[i]);
	  ncresupfd->Fill(i,frw->nncE[i]);
	}
      }
      
      if(frw->res_ma>.87&&frw->res_ma<.93){
	for(int i=0;i<nrw->EBINS;i++){
	  tresdownfd->Fill(i,frw->nbgE[i]-frw->nnuebE[i]/2.);
	  muresdownfd->Fill(i,frw->nnumuE[i]);
	  ebresdownfd->Fill(i,frw->nnuebE[i]/2.);
	  tauresdownfd->Fill(i,frw->nnutauE[i]);
	  ncresdownfd->Fill(i,frw->nncE[i]);
	}
      }
    }

    if(frw->qel_ma==1&&frw->res_ma==1){
      if(frw->kno_r112>1.47&&frw->kno_r112<1.53){
	for(int i=0;i<frw->EBINS;i++){
	  tknoupfd->Fill(i,frw->nbgE[i]-frw->nnuebE[i]/2.);
	  muknoupfd->Fill(i,frw->nnumuE[i]);
	  ebknoupfd->Fill(i,frw->nnuebE[i]/2.);
	  tauknoupfd->Fill(i,frw->nnutauE[i]);
	  ncknoupfd->Fill(i,frw->nncE[i]);
	}
      }
      
      if(frw->kno_r112>.47&&frw->kno_r112<.53){
	for(int i=0;i<nrw->EBINS;i++){
	  tknodownfd->Fill(i,frw->nbgE[i]-frw->nnuebE[i]/2.);
	  muknodownfd->Fill(i,frw->nnumuE[i]);
	  ebknodownfd->Fill(i,frw->nnuebE[i]/2.);
	  tauknodownfd->Fill(i,frw->nnutauE[i]);
	  ncknodownfd->Fill(i,frw->nncE[i]);
	}
      }
    }

  }


tnomnd->SetLineColor(1);
munomnd->SetLineColor(1);
ebnomnd->SetLineColor(1);
taunomnd->SetLineColor(1);
ncnomnd->SetLineColor(1);

tnomfd->SetLineColor(1);
munomfd->SetLineColor(1);
ebnomfd->SetLineColor(1);
taunomfd->SetLineColor(1);
ncnomfd->SetLineColor(1);

tqelupnd->SetLineColor(4);
muqelupnd->SetLineColor(4);
ebqelupnd->SetLineColor(4);
tauqelupnd->SetLineColor(4);
ncqelupnd->SetLineColor(4);

tqeldownnd->SetLineColor(2);
muqeldownnd->SetLineColor(2);
ebqeldownnd->SetLineColor(2);
tauqeldownnd->SetLineColor(2);
ncqeldownnd->SetLineColor(2);

tqelupfd->SetLineColor(4);
muqelupfd->SetLineColor(4);
ebqelupfd->SetLineColor(4);
tauqelupfd->SetLineColor(4);
ncqelupfd->SetLineColor(4);

tqeldownfd->SetLineColor(2);
muqeldownfd->SetLineColor(2);
ebqeldownfd->SetLineColor(2);
tauqeldownfd->SetLineColor(2);
ncqeldownfd->SetLineColor(2);

tresupnd->SetLineColor(4);
muresupnd->SetLineColor(4);
ebresupnd->SetLineColor(4);
tauresupnd->SetLineColor(4);
ncresupnd->SetLineColor(4);

tresdownnd->SetLineColor(2);
muresdownnd->SetLineColor(2);
ebresdownnd->SetLineColor(2);
tauresdownnd->SetLineColor(2);
ncresdownnd->SetLineColor(2);

tresupfd->SetLineColor(4);
muresupfd->SetLineColor(4);
ebresupfd->SetLineColor(4);
tauresupfd->SetLineColor(4);
ncresupfd->SetLineColor(4);

tresdownfd->SetLineColor(2);
muresdownfd->SetLineColor(2);
ebresdownfd->SetLineColor(2);
tauresdownfd->SetLineColor(2);
ncresdownfd->SetLineColor(2);

tknoupnd->SetLineColor(4);
muknoupnd->SetLineColor(4);
ebknoupnd->SetLineColor(4);
tauknoupnd->SetLineColor(4);
ncknoupnd->SetLineColor(4);

tknodownnd->SetLineColor(2);
muknodownnd->SetLineColor(2);
ebknodownnd->SetLineColor(2);
tauknodownnd->SetLineColor(2);
ncknodownnd->SetLineColor(2);

tknoupfd->SetLineColor(4);
muknoupfd->SetLineColor(4);
ebknoupfd->SetLineColor(4);
tauknoupfd->SetLineColor(4);
ncknoupfd->SetLineColor(4);

tknodownfd->SetLineColor(2);
muknodownfd->SetLineColor(2);
ebknodownfd->SetLineColor(2);
tauknodownfd->SetLineColor(2);
ncknodownfd->SetLineColor(2);




 TCanvas *qelcan=new TCanvas("qelcan","qelcan",950,650);
 qelcan->Divide(2,2);
 qelcan->cd(1);
 tnomnd->Draw("e");
 tqelupnd->Draw("same");
 tqeldownnd->Draw("same");
 qelcan->cd(2);
 tnomfd->Draw("e");
 tqelupfd->Draw("same");
 tqeldownfd->Draw("same");
 qelcan->cd(3);
 gPad->Divide(2,2);
 gPad->cd(1);
 munomnd->Draw("e");
 muqelupnd->Draw("same");
 muqeldownnd->Draw("same");
 qelcan->cd(3);
 gPad->cd(2);
 ebnomnd->Draw("e");
 ebqelupnd->Draw("same");
 ebqeldownnd->Draw("same");
 qelcan->cd(3);
 gPad->cd(3);
 taunomnd->Draw("e");
 tauqelupnd->Draw("same");
 tauqeldownnd->Draw("same");
 qelcan->cd(3);
 gPad->cd(4);
 ncnomnd->Draw("e");
 ncqelupnd->Draw("same");
 ncqeldownnd->Draw("same");
 qelcan->cd(4);
 gPad->Divide(2,2);
 gPad->cd(1);
 munomfd->Draw("e");
 muqelupfd->Draw("same");
 muqeldownfd->Draw("same");
 qelcan->cd(4);
 gPad->cd(2);
 ebnomfd->Draw("e");
 ebqelupfd->Draw("same");
 ebqeldownfd->Draw("same");
 qelcan->cd(4);
 gPad->cd(3);
 taunomfd->Draw("e");
 tauqelupfd->Draw("same");
 tauqeldownfd->Draw("same");
 qelcan->cd(4);
 gPad->cd(4);
 ncnomfd->Draw("e");
 ncqelupfd->Draw("same");
 ncqeldownfd->Draw("same");



 TCanvas *rescan=new TCanvas("rescan","rescan",950,650);
 rescan->Divide(2,2);
 rescan->cd(1);
 tnomnd->Draw("e");
 tresupnd->Draw("same");
 tresdownnd->Draw("same");
 rescan->cd(2);
 tnomfd->Draw("e");
 tresupfd->Draw("same");
 tresdownfd->Draw("same");
 rescan->cd(3);
 gPad->Divide(2,2);
 gPad->cd(1);
 munomnd->Draw("e");
 muresupnd->Draw("same");
 muresdownnd->Draw("same");
 rescan->cd(3);
 gPad->cd(2);
 ebnomnd->Draw("e");
 ebresupnd->Draw("same");
 ebresdownnd->Draw("same");
 rescan->cd(3);
 gPad->cd(3);
 taunomnd->Draw("e");
 tauresupnd->Draw("same");
 tauresdownnd->Draw("same");
 rescan->cd(3);
 gPad->cd(4);
 ncnomnd->Draw("e");
 ncresupnd->Draw("same");
 ncresdownnd->Draw("same");
 rescan->cd(4);
 gPad->Divide(2,2);
 gPad->cd(1);
 munomfd->Draw("e");
 muresupfd->Draw("same");
 muresdownfd->Draw("same");
 rescan->cd(4);
 gPad->cd(2);
 ebnomfd->Draw("e");
 ebresupfd->Draw("same");
 ebresdownfd->Draw("same");
 rescan->cd(4);
 gPad->cd(3);
 taunomfd->Draw("e");
 tauresupfd->Draw("same");
 tauresdownfd->Draw("same");
 rescan->cd(4);
 gPad->cd(4);
 ncnomfd->Draw("e");
 ncresupfd->Draw("same");
 ncresdownfd->Draw("same");



 TCanvas *knocan=new TCanvas("knocan","knocan",950,650);
 knocan->Divide(2,2);
 knocan->cd(1);
 tnomnd->Draw("e");
 tknoupnd->Draw("same");
 tknodownnd->Draw("same");
 knocan->cd(2);
 tnomfd->Draw("e");
 tknoupfd->Draw("same");
 tknodownfd->Draw("same");
 knocan->cd(3);
 gPad->Divide(2,2);
 gPad->cd(1);
 munomnd->Draw("e");
 muknoupnd->Draw("same");
 muknodownnd->Draw("same");
 knocan->cd(3);
 gPad->cd(2);
 ebnomnd->Draw("e");
 ebknoupnd->Draw("same");
 ebknodownnd->Draw("same");
 knocan->cd(3);
 gPad->cd(3);
 taunomnd->Draw("e");
 tauknoupnd->Draw("same");
 tauknodownnd->Draw("same");
 knocan->cd(3);
 gPad->cd(4);
 ncnomnd->Draw("e");
 ncknoupnd->Draw("same");
 ncknodownnd->Draw("same");
 knocan->cd(4);
 gPad->Divide(2,2);
 gPad->cd(1);
 munomfd->Draw("e");
 muknoupfd->Draw("same");
 muknodownfd->Draw("same");
 knocan->cd(4);
 gPad->cd(2);
 ebnomfd->Draw("e");
 ebknoupfd->Draw("same");
 ebknodownfd->Draw("same");
 knocan->cd(4);
 gPad->cd(3);
 taunomfd->Draw("e");
 tauknoupfd->Draw("same");
 tauknodownfd->Draw("same");
 knocan->cd(4);
 gPad->cd(4);
 ncnomfd->Draw("e");
 ncknoupfd->Draw("same");
 ncknodownfd->Draw("same");



}

void FitNue(const char *challengefile, const char *rwfile, float ndscale, float fdscale, bool save)
{
   //when you add the extra row for nominal weights, remember to start at index 1
//decided not to do this, so don't worrry about it.

/*
   const int NQEL=3;
   const int NRES=3;
   const int NKNO=3;
   const int NDM2=3;
   const int NSST=3;
*/
   const int NUE32=21;

  TFile *f = new TFile(challengefile);
  TH1F *ndspec = (TH1F *)f->Get("ndspec");
  TH1F *fdspec = (TH1F *)f->Get("fdspec");
  ndspec->Sumw2();
  fdspec->Sumw2();

  TFile *g = new TFile(rwfile);
  TTree *fdt=(TTree *)g->Get("rwtree_FD_norm");
  TTree *ndt=(TTree *)g->Get("rwtree_ND_norm");

  if(fdt->GetEntries()!=ndt->GetEntries()){
    cout<<"TREE sizes don't match.  Won't continue"<<endl;
    return;
  }

  NueRW *frw = new NueRW();
  fdt->SetBranchAddress("NueRW",&frw);
  NueRW *nrw = new NueRW();
  ndt->SetBranchAddress("NueRW",&nrw);

  float minchi2=10000.;
  int minindex=0;
  float minchi2xsnom=10000.;
  int minindexxsnom=0;

  float minnchi2=10000.;
  int minnindex=0;
  float minfchi2=10000.;
  int minfindex=0;

  int NROWS=(int)(fdt->GetEntries());

  double qelarr[NROWS];
  double resarr[NROWS];
  double knoarr[NROWS];
  double dmarr[NROWS];
  double ssarr[NROWS];
  double ue32arr[NROWS];
  double chi2arr[NROWS];
  double ndofarr[NROWS];
  double nchi2arr[NROWS];
  double nndofarr[NROWS];
  double fchi2arr[NROWS];
  double fndofarr[NROWS];
  double chi2xsnom[NROWS];
  double ndofxsnom[NROWS];

  double index[NROWS];
  
  int indexminchi2ue320=0;
  int minndfue320=0;
  double minchi2ue320=10000.;
  double minchi2onlyue32[NUE32]={10000.};
  double ue32shortarr[NUE32];
  double farminchi2onlyue32[NUE32]={10000.};
  for(int i=0;i<NUE32;i++){
     minchi2onlyue32[i]=10000.;
     farminchi2onlyue32[i]=10000.;
     ue32shortarr[i]=0.;
  }

  int ndnom=0;
  int fdnom=0;
  for(int i=0;i<NROWS;i++){
/*
     int qelit=(int)(i/(NRES*NKNO*NDM2*NSST*NUE32));
     int resit=(int)((i-qelit*(NRES*NKNO*NDM2*NSST*NUE32))/(NKNO*NDM2*NSST*NUE32));
     int knoit=(int)((i-qelit*(NRES*NKNO*NDM2*NSST*NUE32)-resit*(NKNO*NDM2*NSST*NUE32))/(NDM2*NSST*NUE32));
     int dmit=(int)((i-qelit*(NRES*NKNO*NDM2*NSST*NUE32)
		     -resit*(NKNO*NDM2*NSST*NUE32)
		     -knoit*(NDM2*NSST*NUE32))/(NSST*NUE32));
     int ssit=(int)((i-qelit*(NRES*NKNO*NDM2*NSST*NUE32)
		      -resit*(NKNO*NDM2*NSST*NUE32)
		      -knoit*(NDM2*NSST*NUE32)
		      -dmit*(NSST*NUE32))/(NUE32));
*/
     int ueit=(i%NUE32);


//    cout<<"*********ON ENTRY "<<i<<"******************"<<endl;
//    cout<<qelit<<" "<<resit<<" "<<knoit<<" "<<dmit<<" "<<ssit<<" "<<ueit<<endl;


    index[i]=1.*i;
    fdt->GetEntry(i);
    ndt->GetEntry(i);

    TH1F ndh("ndh","ndh",frw->EBINS,0,frw->EBINS*frw->EBINW);
    TH1F fdh("fdh","fdh",nrw->EBINS,0,nrw->EBINS*nrw->EBINW);

    for(int j=0;j<frw->EBINS;j++){
      ndh.SetBinContent(j+1,nrw->nbgE[j]*ndscale);
      fdh.SetBinContent(j+1,(frw->nbgE[j]+frw->nsigE[j])*fdscale);
    }

    if(nrw->qel_ma==1&&nrw->res_ma==1&&nrw->kno_r112==1){
       ndnom=i;
    }

    int ndf=0, fndf=0, nndf=0;
    double chi2=0,nchi2=0,fchi2=0;
//    ComputeChi2(&ndh,ndspec,nchi2,nndf,0);
//    ComputeChi2(&fdh,fdspec,fchi2,fndf,0);
    ComputeChi2(&ndh,ndspec,nchi2,nndf,1);
    ComputeChi2(&fdh,fdspec,fchi2,fndf,1);
    
    chi2 = fchi2 + nchi2;
    ndf = fndf + nndf;

    if(frw->qel_ma==1&&frw->res_ma==1&&frw->kno_r112==1&&
       fabs(frw->dm2-2.175e-3)<1e-6&&fabs(frw->ss2th-0.925)<1e-6){
       chi2xsnom[i]=chi2;
       ndofxsnom[i]=ndf;
       if(minchi2xsnom>chi2){
	  minchi2xsnom=chi2;
	  minindexxsnom=i;
       }
       if(fabs(frw->UE32)<0.000001){
	  fdnom=i;
	  cout<<"fd nominal index "<<i<<endl;
       }
    }
    else{
       chi2xsnom[i]=-99.;
    }
    if(fabs(frw->UE32)<0.00000001){
       if(minchi2ue320>chi2){
	  minchi2ue320=chi2;
	  indexminchi2ue320=i;
	  minndfue320=ndf;
       }
    }

//    ue32arr[i]=frw->UE32;
    ue32arr[i]=4*frw->UE32*(1-frw->UE32);

    if(chi2<minchi2onlyue32[ueit]){
       minchi2onlyue32[ueit]=chi2;
//       ue32shortarr[ueit]=frw->UE32;
       ue32shortarr[ueit]=ue32arr[i];
    }
    if(fchi2<farminchi2onlyue32[ueit]){
       farminchi2onlyue32[ueit]=fchi2;
    }

    qelarr[i]=frw->qel_ma;
    resarr[i]=frw->res_ma;
    knoarr[i]=frw->kno_r112;
    dmarr[i]=frw->dm2;
    ssarr[i]=frw->ss2th;
    chi2arr[i]=chi2;
    ndofarr[i]=1.*ndf;
    nchi2arr[i]=nchi2;
    nndofarr[i]=1.*nndf;
    fchi2arr[i]=fchi2;
    fndofarr[i]=1.*fndf;

    if(chi2<minchi2){
//      cout<<"New min: "<<chi2<<" i "<<i<<endl;
      minchi2=chi2;
      minindex=i;
    }
    if(fchi2<minfchi2){
//      cout<<"New far min: "<<fchi2<<" i "<<i<<endl;
      minfchi2=fchi2;
      minfindex=i;
    }
    if(nchi2<minnchi2){
//      cout<<"New far min: "<<nchi2<<" i "<<i<<endl;
      minnchi2=nchi2;
      minnindex=i;
    }
  }    				     
  cout<<"MININDEX xs nom "<<minindexxsnom<<endl;

  cout<<"minimum chi2 for UE32==0 "<<minchi2ue320<<" index "<<indexminchi2ue320<<endl;
  TF1 *probfunc = new TF1("prob",Chi2Prob,0,5000,1);
  probfunc->SetParameter(0,minndfue320);
  cout<<"probability "<<1-probfunc->Integral(0,minchi2ue320)<<endl;


  TGraph *chi2vue32min = new TGraph(NUE32,ue32shortarr,minchi2onlyue32);
  TGraph *farchi2vue32min = new TGraph(NUE32,ue32shortarr,farminchi2onlyue32);
  TGraph *chi2vue32mincopy = new TGraph();
  TGraph *farchi2vue32mincopy = new TGraph();
  float tmpmin = 10000;
  float fartmpmin = 10000;
  float ue32answer=0.;
  for(int i=0;i<NUE32;i++){
     if(tmpmin>minchi2onlyue32[i]){
	tmpmin=minchi2onlyue32[i];
	ue32answer=ue32shortarr[i];
     }
     if(fartmpmin>farminchi2onlyue32[i]){
	fartmpmin=farminchi2onlyue32[i];
     }
  }
  cout<<"tmpmin "<<tmpmin<<" fartmpmin "<<fartmpmin<<endl;
  for(int i=0;i<NUE32;i++){
     chi2vue32mincopy->SetPoint(i,ue32shortarr[i],minchi2onlyue32[i]-tmpmin);
     farchi2vue32mincopy->SetPoint(i,ue32shortarr[i],farminchi2onlyue32[i]-fartmpmin);
  }


  TGraph *chi2vue32 = new TGraph(NROWS,ue32arr,chi2arr);
  TGraph *chi2vrow = new TGraph(NROWS,index,chi2arr);
  TGraph *ndfvrow = new TGraph(NROWS,index,ndofarr);

  TGraph *nchi2vue32 = new TGraph(NROWS,ue32arr,nchi2arr);
  TGraph *nchi2vrow = new TGraph(NROWS,index,nchi2arr);
  TGraph *nndfvrow = new TGraph(NROWS,index,nndofarr);

  TGraph *fchi2vue32 = new TGraph(NROWS,ue32arr,fchi2arr);
  TGraph *fchi2vrow = new TGraph(NROWS,index,fchi2arr);
  TGraph *fndfvrow = new TGraph(NROWS,index,fndofarr);

  TGraph *chi2vqelbest = new TGraph();
  TGraph *chi2vresbest = new TGraph();
  TGraph *chi2vknobest = new TGraph();
  TGraph *chi2vdm2best = new TGraph();
  TGraph *chi2vssbest = new TGraph();
  TGraph *chi2vue32best = new TGraph();
  TGraph *chi2vue32nomxs =  new TGraph();
  TGraph *chi2vue32nomall = new TGraph();
  TGraph *farchi2vue32nomall = new TGraph();
  TGraph *chi2vue32nomallcopy = new TGraph();
  TGraph *farchi2vue32nomallcopy = new TGraph();


  int ptqel=0;
  int ptres=0;
  int ptkno=0;
  int ptdm2=0;
  int ptss=0;
  int ptueb=0;
  int ptuenomall=0;
  int ptuenxs=0;

  float tmpminnom=10000.;
  float fartmpminnom=10000.;
  for(int i=0;i<NROWS;i++){
     //plot chi2 v qel for best fit of other params
     if(fabs(resarr[i]-resarr[minindex])<0.001&&fabs(knoarr[i]-knoarr[minindex])<0.001
	&&fabs(dmarr[i]-dmarr[minindex])<0.000001&&fabs(ssarr[i]-ssarr[minindex])<0.001
	&&fabs(ue32arr[i]-ue32arr[minindex])<0.001){
	chi2vqelbest->SetPoint(ptqel,qelarr[i],chi2arr[i]);
	ptqel++;
     }
     //plot chi2 v res for best fit of other params
     if(fabs(qelarr[i]-qelarr[minindex])<0.001&&fabs(knoarr[i]-knoarr[minindex])<0.001
	&&fabs(dmarr[i]-dmarr[minindex])<0.000001&&fabs(ssarr[i]-ssarr[minindex])<0.001
	&&fabs(ue32arr[i]-ue32arr[minindex])<0.001){
	chi2vresbest->SetPoint(ptres,resarr[i],chi2arr[i]);
	ptres++;
     }
     //plot chi2 v kno for best fit of other params
     if(fabs(qelarr[i]-qelarr[minindex])<0.001&&fabs(resarr[i]-resarr[minindex])<0.001
	&&fabs(dmarr[i]-dmarr[minindex])<0.000001&&fabs(ssarr[i]-ssarr[minindex])<0.001
	&&fabs(ue32arr[i]-ue32arr[minindex])<0.001){
	chi2vknobest->SetPoint(ptkno,knoarr[i],chi2arr[i]);
	ptkno++;
     }
     //plot chi2 v dm for best fit of other params
     if(fabs(qelarr[i]-qelarr[minindex])<0.001&&fabs(resarr[i]-resarr[minindex])<0.001
	&&fabs(knoarr[i]-knoarr[minindex])<0.001&&fabs(ssarr[i]-ssarr[minindex])<0.001
	&&fabs(ue32arr[i]-ue32arr[minindex])<0.001){
	chi2vdm2best->SetPoint(ptdm2,dmarr[i],chi2arr[i]);
	ptdm2++;
     }
     //plot chi2 v ss for best fit of other params
     if(fabs(qelarr[i]-qelarr[minindex])<0.001&&fabs(resarr[i]-resarr[minindex])<0.001
	&&fabs(knoarr[i]-knoarr[minindex])<0.001&&fabs(dmarr[i]-dmarr[minindex])<0.000001
	&&fabs(ue32arr[i]-ue32arr[minindex])<0.001){
	chi2vssbest->SetPoint(ptss,ssarr[i],chi2arr[i]);
	ptss++;
     }
     //plot chi2 v ue32 for best fit of other params
     if(fabs(qelarr[i]-qelarr[minindex])<0.001&&fabs(resarr[i]-resarr[minindex])<0.001
	&&fabs(knoarr[i]-knoarr[minindex])<0.001&&fabs(dmarr[i]-dmarr[minindex])<0.000001
	&&fabs(ssarr[i]-ssarr[minindex])<0.001){
	chi2vue32best->SetPoint(ptueb,ue32arr[i],chi2arr[i]);
	ptueb++;
     }
     //plot chi2 v ue32 for nominal values of other params
     if(fabs(qelarr[i]-qelarr[fdnom])<0.001&&fabs(resarr[i]-resarr[fdnom])<0.001
	&&fabs(knoarr[i]-knoarr[fdnom])<0.001&&fabs(dmarr[i]-dmarr[fdnom])<0.000001
	&&fabs(ssarr[i]-ssarr[fdnom])<0.001){
//	cout<<"filling chi2vue32nomall "<<ue32arr[i]<<" "<<chi2xsnom[i]<<endl;
	chi2vue32nomall->SetPoint(ptuenomall,ue32arr[i],chi2arr[i]);
	farchi2vue32nomall->SetPoint(ptuenomall,ue32arr[i],fchi2arr[i]);
	if(chi2arr[i]<tmpminnom){
	   tmpminnom=chi2arr[i];
	}
	if(fchi2arr[i]<fartmpminnom){
	   fartmpminnom=fchi2arr[i];
	}
	ptuenomall++;
     }
     //plot chi2 v ue32 for best fit when only osc params vary
     if(fabs(qelarr[i]-qelarr[minindexxsnom])<0.001&&fabs(resarr[i]-resarr[minindexxsnom])<0.001
	&&fabs(knoarr[i]-knoarr[minindexxsnom])<0.001){
        if(chi2xsnom[i]>0){
	   chi2vue32nomxs->SetPoint(ptuenxs,ue32arr[i],chi2xsnom[i]);
	   ptuenxs++;
	}
     }
  }
  cout<<"tmpminnom "<<tmpminnom<<" fartmpminnom "<<fartmpminnom<<endl;
  for(int i=0;i<chi2vue32nomall->GetN();i++){
     chi2vue32nomallcopy->SetPoint(i,chi2vue32nomall->GetX()[i],chi2vue32nomall->GetY()[i]-tmpminnom);
     farchi2vue32nomallcopy->SetPoint(i,farchi2vue32nomall->GetX()[i],farchi2vue32nomall->GetY()[i]-fartmpminnom);
  }

  
  chi2vue32min->SetNameTitle("chi2vss2thetamin",
			     "Best #chi^{2} for a given sin^{2}(2#theta_{13});sin^{2}(2#theta_{13});#chi^{2}");
  farchi2vue32min->SetNameTitle("chi2vss2thetamin",
				"Best #chi^{2} for a given sin^{2}(2#theta_{13}), Far Only;sin^{2}(2#theta_{13});#chi^{2}");
  chi2vue32mincopy->SetNameTitle("deltachi2vssthetamin",
				 "Min  #Delta #chi^{2} for a given sin^{2}(2#theta_{13});sin^{2}(2#theta_{13});#Delta #chi^{2}");
  farchi2vue32mincopy->SetNameTitle("deltafarchi2vss2thetamin",
				    "Min  #Delta #chi^{2} for a given sin^{2}(2#theta_{13}), Far Only;sin^{2}(2#theta_{13});#Delta #chi^{2}");

  chi2vue32->SetNameTitle("chi2vss2thetha","All #chi^{2} vs. sin^{2}(2#theta_{13});sin^{2}(2#theta_{13});#chi^{2}");
  chi2vrow->SetNameTitle("chi2vrow","#chi^{2} vs. Row;Row index;#chi^{2}");
  ndfvrow->SetNameTitle("ndfvrow","NDF vs. Row;Row index;NDF");
  
  nchi2vue32->SetNameTitle("nchi2vss2theta","Near #chi^{2} vs. sin^{2}(2#theta_{13});sin^{2}(2#theta_{13});#chi^{2}");
  nchi2vrow->SetNameTitle("nchi2vrow","Near #chi^{2} vs. Row;Row index;#chi^{2}");
  nndfvrow->SetNameTitle("nndfvrow","Near NDF vs. Row;Row index;NDF");
  
  fchi2vue32->SetNameTitle("fchi2vss2theta","Far #chi^{2} vs. sin^{2}(2#theta_{13});sin^{2}(2#theta_{13});#chi^{2}");
  fchi2vrow->SetNameTitle("fchi2vrow","Far #chi^{2} vs. Row;Row index;#chi^{2}");
  fndfvrow->SetNameTitle("fndfvrow","Far NDF vs. Row;Row index;NDF");
  
  chi2vqelbest->SetNameTitle("chi2vqelbest","#chi^{2} vs. QEL;QEL;#chi^{2}");
  chi2vresbest->SetNameTitle("chi2vresbest","#chi^{2} vs. RES;RES;#chi^{2}");
  chi2vknobest->SetNameTitle("chi2vknobest","#chi^{2} vs. KNO;KNO;#chi^{2}");
  chi2vdm2best->SetNameTitle("chi2vdm2best","#chi^{2} vs. #Delta m^{2};#Delta m^{2} (eV^{2});#chi^{2}");
  chi2vssbest->SetNameTitle("chi2vssbest","#chi^{2} vs. sin^{2}(2#theta_{23});sin^{2}(2#theta_{23});#chi^{2}");
  chi2vue32best->SetNameTitle("chi2vss2thetabest","#chi^{2} vs. sin^{2}(2#theta_{13});sin^{2}(2#theta_{13});#chi^{2}");
  chi2vue32nomall->SetNameTitle("chi2vss2thetanomall","#chi^{2} vs. sin^{2}(2#theta_{13}), All parameters nominal;sin^{2}(2#theta_{13});#chi^{2}");
  farchi2vue32nomall->SetNameTitle("chi2vss2thetanomall",
				   "#chi^{2} vs. sin^{2}(2#theta_{13}), All parameters nominal, Far Only;sin^{2}(2#theta_{13});#chi^{2}");
  chi2vue32nomallcopy->SetNameTitle("delatchi2vss2thetanomall",
				    "#Delta #chi^{2} vs. sin^{2}(2#theta_{13}), All parameters nominal;sin^{2}(2#theta_{13});#Delta #chi^{2}");
  farchi2vue32nomallcopy->SetNameTitle("fardeltachi2vss2thetanomall","#Delta #chi^{2} vs. sin^{2}(2#theta_{13}), All parameters nominal, Far Only;sin^{2}(2#theta_{13});#Delta #chi^{2}");


  chi2vue32nomxs->SetNameTitle("chi2vss2thetanomxs",
			       "#chi^{2} vs. sin^{2}(2#theta_{13}), X-section parameters nominal;sin^{2}(2#theta_{13});#chi^{2}");
  

  fdt->GetEntry(minnindex);
/*
  cout<<"min NEAR chi2 "<<minnchi2<<" nndof "<<nndofarr[minindex]<<" index "<<minnindex<<endl;
  cout<<"dm2 "<<frw->dm2<<" thetha "<<frw->ss2th
      <<"qelma "<<frw->qel_ma<<" res "<<frw->res_ma
      <<"kno "<<frw->kno_r112<<" UE32 "<<frw->UE32<<endl;
*/
  fdt->GetEntry(minfindex);
/*
  cout<<"min FAR chi2 "<<minfchi2<<" ndof "<<fndofarr[minindex]<<" index "<<minfindex<<endl;
  cout<<"dm2 "<<frw->dm2<<" thetha "<<frw->ss2th
      <<"qelma "<<frw->qel_ma<<" res "<<frw->res_ma
      <<"kno "<<frw->kno_r112<<" UE32 "<<frw->UE32<<endl;
*/
  fdt->GetEntry(minindex);

  cout<<"TOTAL minchi2 "<<minchi2<<" ndof "<<ndofarr[minindex]<<" index "<<minindex<<endl;
  cout<<"dm2 "<<frw->dm2<<" thetha "<<frw->ss2th
      <<"qelma "<<frw->qel_ma<<" res "<<frw->res_ma
      <<"kno "<<frw->kno_r112<<" UE32 "<<frw->UE32<<endl
      <<" tot sig "<<frw->nsig<<" "<<frw->nbg<<endl;
  
  ndt->GetEntry(minindex);
  TH1F *ndh=new TH1F("ndhbest","Reco E_{#nu}, ND, best fit",frw->EBINS,0,frw->EBINS*frw->EBINW);
  TH1F *fdh=new TH1F("fdhbest","Reco E_{#nu}, FD, best fit",nrw->EBINS,0,nrw->EBINS*nrw->EBINW);
  TH1F *fdsh=new TH1F("fdshbest","Reco E_{#nu}, Signal, best fit",nrw->EBINS,0,nrw->EBINS*nrw->EBINW);
  
  for(int j=0;j<frw->EBINS;j++){
     ndh->SetBinContent(j+1,nrw->nbgE[j]*ndscale);
     fdh->SetBinContent(j+1,(frw->nbgE[j]+frw->nsigE[j])*fdscale);
     fdsh->SetBinContent(j+1,(frw->nsigE[j])*fdscale);
  }
  TH1F *ndnh=new TH1F("ndhnom","Reco E_{#nu}, ND, nominal",frw->EBINS,0,frw->EBINS*frw->EBINW);
  ndt->GetEntry(ndnom);
  for(int j=0;j<frw->EBINS;j++){
     ndnh->SetBinContent(j+1,nrw->nbgE[j]*ndscale);
  }
  TH1F *fdnh=new TH1F("fdhnom","Reco E_{#nu}, FD, nominal",frw->EBINS,0,frw->EBINS*frw->EBINW);
  fdt->GetEntry(fdnom);
  for(int j=0;j<frw->EBINS;j++){
     fdnh->SetBinContent(j+1,(frw->nbgE[j]+frw->nsigE[j])*fdscale);
  }

  ndspec->SetTitle("Near Detector, Reco E_{#nu};E (GeV);Events/GeV");
  fdspec->SetTitle("Far Detector, Reco E_{#nu};E (GeV);Events/GeV");
  
  ndspec->SetLineColor(1);
  ndh->SetLineColor(4);
  ndnh->SetLineColor(2);
  fdspec->SetLineColor(1);
  fdh->SetLineColor(4);
  fdnh->SetLineColor(2);
  fdsh->SetLineColor(6);
  ndspec->SetStats(0);
  ndnh->SetStats(0);
  fdspec->SetStats(0);
  fdnh->SetStats(0);
     

  TCanvas *canfit = new TCanvas("canfit","canfit");
  canfit->Divide(1,2);
  canfit->cd(1);
  ndspec->Draw("e");
  ndspec->GetXaxis()->CenterTitle();
  ndspec->GetYaxis()->CenterTitle();
  ndh->Draw("same");
  ndnh->Draw("same");
  canfit->cd(2);
  fdspec->Draw("e");
  fdh->Draw("same");
  fdsh->Draw("same");
  fdnh->Draw("same");
  fdspec->GetXaxis()->CenterTitle();
  fdspec->GetYaxis()->CenterTitle();


  TCanvas *chi2can = new TCanvas("chi2can","chi2can",900,600);
  chi2can->Divide(3,2);
  chi2can->cd(1);
  chi2vqelbest->Draw("AP");
  chi2vqelbest->GetXaxis()->CenterTitle();
  chi2vqelbest->GetYaxis()->CenterTitle();
  chi2can->cd(2);
  chi2vresbest->Draw("AP");
  chi2vresbest->GetXaxis()->CenterTitle();
  chi2vresbest->GetYaxis()->CenterTitle();
  chi2can->cd(3);
  chi2vknobest->Draw("AP");
  chi2vknobest->GetXaxis()->CenterTitle();
  chi2vknobest->GetYaxis()->CenterTitle();
  chi2can->cd(4);
  chi2vdm2best->Draw("AP");
  chi2vdm2best->GetXaxis()->CenterTitle();
  chi2vdm2best->GetYaxis()->CenterTitle();
  chi2can->cd(5);
  chi2vssbest->Draw("AP");
  chi2vssbest->GetXaxis()->CenterTitle();
  chi2vssbest->GetYaxis()->CenterTitle();
  chi2can->cd(6);
  chi2vue32best->Draw("AP");
  chi2vue32best->GetXaxis()->CenterTitle();
  chi2vue32best->GetYaxis()->CenterTitle();

  farchi2vue32nomall->SetMarkerColor(6);
  farchi2vue32min->SetMarkerColor(2);
  chi2vue32nomall->SetMarkerColor(4);
  chi2vue32min->SetMarkerColor(1);
  farchi2vue32nomallcopy->SetMarkerColor(6);
  farchi2vue32mincopy->SetMarkerColor(2);
  chi2vue32nomallcopy->SetMarkerColor(4);
  chi2vue32mincopy->SetMarkerColor(1);

  float ub68,ub90, ub99;
  float lb68,lb90, lb99;

/*
  FindUpperBounds(farchi2vue32nomallcopy,b68,b90,b99);
  cout<<"FAR only, NOM pars: 68%: "<<b68<<" 90%: "<<b90<<" 99% "<<b99<<endl;
  FindUpperBounds(farchi2vue32mincopy,b68,b90,b99);
  cout<<"Far only, All Pars 68%: "<<b68<<" 90%: "<<b90<<" 99% "<<b99<<endl;
*/

  FindUpperBounds(chi2vue32mincopy,ub68,ub90,ub99);
  cout<<"THE ANSWER: 68% "<<ub68<<" 90%: "<<ub90<<" 99% "<<ub99<<endl;
  FindLowerBounds(chi2vue32mincopy,lb68,lb90,lb99);
  cout<<"THE ANSWER: 68% "<<lb68<<" 90%: "<<lb90<<" 99% "<<lb99<<endl;

  TCanvas *uecan = new TCanvas("uecan","uecan");
  uecan->Divide(2,1);
  uecan->cd(1);
  farchi2vue32nomall->Draw("AP");
  farchi2vue32nomall->GetXaxis()->CenterTitle();
  farchi2vue32nomall->GetYaxis()->CenterTitle();
  farchi2vue32min->Draw("p");
  uecan->cd(2);
  chi2vue32nomall->Draw("AP");
  chi2vue32nomall->GetXaxis()->CenterTitle();
  chi2vue32nomall->GetYaxis()->CenterTitle();
  chi2vue32min->Draw("P");

  TCanvas *uecan2 = new TCanvas("uecan2","uecan2");
  uecan2->cd();
  farchi2vue32nomallcopy->Draw("AP");
  farchi2vue32nomallcopy->GetXaxis()->CenterTitle();
  farchi2vue32nomallcopy->GetYaxis()->CenterTitle();
  farchi2vue32mincopy->Draw("p");
  chi2vue32nomallcopy->Draw("P");
  chi2vue32mincopy->Draw("P");

  cout<<"*************************************************"<<endl;
  cout<<"MIN chi2 of min chi2 over ue32: "<<tmpmin<<" at "<<ue32answer<<endl;
  cout<<" 68% ubound is: "<<ub68<<endl;
  cout<<" 68% lbound is: "<<lb68<<endl;
  cout<<" 90% ubound is: "<<ub90<<endl;
  cout<<" 90% lbound is: "<<lb90<<endl;
  cout<<" 99% ubound is: "<<ub99<<endl;
  cout<<" 99% lbound is: "<<lb99<<endl;


  TCanvas *uecan3 = new TCanvas("uecan3","uecan3");
  uecan3->cd();
  chi2vue32mincopy->Draw("AP");
  chi2vue32mincopy->GetXaxis()->CenterTitle();
  chi2vue32mincopy->GetYaxis()->CenterTitle();
  chi2vue32mincopy->SetMinimum(-0.1);
  TLine *ul1 = new TLine(ub90,-0.1,ub90,2.71);
  ul1->SetLineColor(1);
  ul1->SetLineWidth(2);
  TLine *ul2 = new TLine(ub99,-0.1,ub99,6.63);
  ul2->SetLineColor(2);
  ul2->SetLineWidth(2);
  TLine *ll1 = new TLine(lb90,-0.1,lb90,2.71);
  ll1->SetLineColor(1);
  ll1->SetLineWidth(2);
  TLine *ll2 = new TLine(lb99,-0.1,lb99,6.63);
  ll2->SetLineColor(2);
  ll2->SetLineWidth(2);

  TLine *l90 = new TLine(0,2.71,ub90,2.71);
  l90->SetLineColor(1);
  l90->SetLineWidth(2);
  TLine *l99 = new TLine(0,6.63,ub99,6.63);
  l99->SetLineColor(2);
  l99->SetLineWidth(2);

  float maxrange=chi2vue32mincopy->GetXaxis()->GetXmax();

  if(ub90<maxrange&&ub90>0){
     ul1->Draw();
  }
  if(ub99<maxrange&&ub90>0){     
     ul2->Draw();
  }
  if(lb90>0){
     ll1->Draw();
  }
  if(lb99>0){
     ll2->Draw();
  }
  float maxyrange=chi2vue32mincopy->GetYaxis()->GetXmax();
  if(2.71<maxyrange){
     l90->Draw();
  }
  if(6.63<maxyrange){
     l99->Draw();
  }

  if(save){
     canfit->SaveAs("spectra.ps");
     canfit->SaveAs("spectra.eps");
     canfit->SaveAs("spectra.gif");
     chi2can->SaveAs("chi2vpar.ps");
     chi2can->SaveAs("chi2vpar.eps");
     chi2can->SaveAs("chi2vpar.gif");
     uecan->SaveAs("ndfdsys.ps");
     uecan->SaveAs("ndfdsys.eps");
     uecan->SaveAs("ndfdsys.gif");
     uecan2->SaveAs("deltachi2.ps");
     uecan2->SaveAs("deltachi2.eps");
     uecan2->SaveAs("deltachi2.gif");
     uecan3->SaveAs("theanswer.ps");
     uecan3->SaveAs("theanswer.eps");
     uecan3->SaveAs("theanswer.gif");
  }


  TFile *out=new TFile("MDC-finalfithists.root","RECREATE");
  chi2vqelbest->Write();
  chi2vresbest->Write();
  chi2vknobest->Write();
  chi2vdm2best->Write();
  chi2vssbest->Write();
  chi2vue32best->Write();
  farchi2vue32nomall->Write();
  farchi2vue32min->Write();
  chi2vue32nomall->Write();
  chi2vue32min->Write();

  farchi2vue32nomallcopy->Write();
  farchi2vue32mincopy->Write();
  chi2vue32nomallcopy->Write();
  chi2vue32mincopy->Write();

  chi2vue32mincopy->Write();
  ul1->Write();
  ul2->Write();
  ll1->Write();
  ll2->Write();

  out->Write();

     
}

void DeltaPlot(const char *rwfile, float scale)
{

  TFile *g = new TFile(rwfile);
  TTree *fdt=(TTree *)g->Get("rwtree_FD_norm");

  NueRW *frw = new NueRW();
  fdt->SetBranchAddress("NueRW",&frw);

  const int NROWS=(int)(fdt->GetEntries());

  const int nomindex=0;
  fdt->GetEntry(nomindex);
  TH1F *nomhist=new TH1F("nomhist","nomhist",frw->EBINS,0,frw->EBINS*frw->EBINW);
  for(int j=0;j<frw->EBINS;j++){
     nomhist->SetBinContent(j+1,(frw->nbgE[j]+frw->nsigE[j])*scale);
  }

  vector<double> ue32norm;
  vector<double> ue32invert;
  vector<double> deltanorm;
  vector<double> deltainvert;
  vector<double> chi2norm;
  vector<double> chi2invert;

  double maxue32=0.;
  double maxdelta=0.;
     
  for(int i=0;i<NROWS;i++){
     fdt->GetEntry(i);
     if(i%1000==0){
	cout<<"on entry "<<i<<endl;
     }
     
     TH1F fdh("fdh","fdh",frw->EBINS,0,frw->EBINS*frw->EBINW);
     for(int j=0;j<frw->EBINS;j++){
	fdh.SetBinContent(j+1,(frw->nbgE[j]+frw->nsigE[j])*scale);
     }
     double chi2=0.;
     int ndf=0;
     ComputeChi2(nomhist,&fdh,chi2,ndf,1);
     if(frw->heirarchy==1){
	ue32norm.push_back(frw->UE32);
	deltanorm.push_back(frw->delta);
	chi2norm.push_back(chi2);
     }
     else if(frw->heirarchy==-1){
	ue32invert.push_back(frw->UE32);
	deltainvert.push_back(frw->delta);
	chi2invert.push_back(chi2);
     }
//     cout<<"on row "<<i<<" heirarchy "<<frw->heirarchy
//	 <<" delta "<<frw->delta<<" ue32 "<<frw->UE32
//	 <<" chi2 "<<chi2<<" ndof "<<ndf<<endl;
     if(frw->UE32>maxue32){
	maxue32=frw->UE32;
     }
     if(frw->delta>maxdelta){
	maxdelta=frw->delta;
     }
  }

  unsigned int nue32steps=(int)(sqrt(1.*ue32norm.size()));
  if(sqrt(1.*ue32invert.size())>nue32steps){
     nue32steps=(int)(sqrt(1.*ue32invert.size()));
  }
  unsigned int ndeltasteps=(int)(sqrt(1.*deltanorm.size()));
  if(deltainvert.size()>ndeltasteps){
     ndeltasteps=(int)(sqrt(1.*deltainvert.size()));
  }


  TH2F *deltacontnorm = new TH2F("deltacontnorm","90% Sensitivity, 7.4e20 POT,  #nu_{#mu} Oscillation Parameters, FD Only;U_{e3}^{2};#delta",
				 nue32steps,0,maxue32,ndeltasteps,0,maxdelta);
  TH2F *deltacontinvert = new TH2F("deltacontinvert","deltacontinvert",
				 nue32steps,0,maxue32,ndeltasteps,0,maxdelta);

  cout<<"Filling norm histos "<<ue32norm.size()<<endl;
  for(unsigned int i=0;i<ue32norm.size();i++){
     deltacontnorm->Fill(ue32norm[i],deltanorm[i],chi2norm[i]);
  }
  cout<<"Filling invert histos "<<ue32invert.size()<<endl;
  for(unsigned int i=0;i<ue32invert.size();i++){
     deltacontinvert->Fill(ue32invert[i],deltainvert[i],chi2invert[i]);
  }
     

  double levnorm[1] = {2.71};
  levnorm[0]+=deltacontnorm->GetMinimum();
  double levinvert[1] = {2.71};
  levinvert[0]+=deltacontinvert->GetMinimum();

  deltacontnorm->SetContour(1,levnorm);
  deltacontinvert->SetContour(1,levinvert);
  deltacontnorm->SetLineColor(2);
  deltacontinvert->SetLineColor(4);

  deltacontnorm->SetStats(0);
  deltacontinvert->SetStats(0);
  TCanvas *deltacan = new TCanvas("deltacan","deltacan");
//  deltacan->Divide(1,2);
  deltacan->cd(1);
  deltacontnorm->Draw("cont2");
  deltacontnorm->GetXaxis()->CenterTitle();
  deltacontnorm->GetXaxis()->SetRangeUser(0,.096);
  deltacontnorm->GetYaxis()->SetRangeUser(0,6.1);
  deltacontnorm->GetYaxis()->CenterTitle();
  
////  deltacontnorm->Draw("colz");
//  deltacan->cd(2);
  deltacontinvert->Draw("cont2 same");
//  deltacontinvert->Draw("colz");
}

void Sensitivity(const char *rwfile,Double_t input_ss2th)
{

  Double_t POT = 7.4;

  TFile *g = new TFile(rwfile);
  TTree *fdt=(TTree *)g->Get("rwtree_FD_norm");
  TTree *ndt=(TTree *)g->Get("rwtree_ND_norm");

  if(fdt->GetEntries()!=ndt->GetEntries()){
    cout<<"TREE sizes don't match.  Won't continue"<<endl;
    return;
  }

  NueRW *frw = new NueRW();
  fdt->SetBranchAddress("NueRW",&frw);
  NueRW *nrw = new NueRW();
  ndt->SetBranchAddress("NueRW",&nrw);


  //for standard fitting trees:

  Int_t ndelta = 3;
  float deltam[3] = {0.00202275,0.002175,0.00232725};
  float xbins[4]  = {};
  for(int i=0;i<=ndelta;i++){
    xbins[i]= deltam[0] + (deltam[1]-deltam[0])*(float(i)-0.5);
  }
  Int_t ntheta13 = 21;
  float ybins[22] = {-0.0025,0.0025,0.0075,0.0125,0.0175,0.0225,0.0275,
		     0.0325,0.0375,0.0425,0.0475,0.0525,0.0575,0.0625,
		     0.0675,0.0725,0.0775,0.0825,0.0875,0.0925,0.0975,
		     0.1025};  
  for(int i=1;i<22;i++) ybins[i] = 4.*(ybins[i])*(1-(ybins[i]));
  ybins[0] = -ybins[1];

  //for special sensitivity plot trees:
  /*  
  Int_t ndelta = 17;
  float deltam[17] = {0};
  for(int i=0;i<ndelta;i++){
    deltam[i] = 0.0009570 + float(i)*(1.5224375e-4);
  }
  float xbins[18]  = {0};
  for(int i=0;i<=ndelta;i++){
    xbins[i]= deltam[0] + (deltam[1]-deltam[0])*(float(i)-0.5);
  }
  Int_t ntheta13 = 17;
  float ybins[18] = {-0.0025,0.0025,0.0075,0.0125,0.0175,0.0225,0.0275,
		     0.0325,0.0375,0.0425,0.0475,0.0525,0.0575,0.0625,
		     0.0675,0.0725,0.0775,0.0825};
  for(int i=1;i<=ntheta13;i++) ybins[i] = 4.*(ybins[i])*(1-(ybins[i]));
  ybins[0] = -ybins[1];
  */

  Double_t nom_UE32 = 0;
  Double_t nom_ss2th = 0.925;
  if(input_ss2th>0) nom_ss2th = input_ss2th;
  Double_t nom_dm2 = 0.002175;
  Double_t nom_dis = 1.0;
  Double_t nom_qel = 1.0;
  Double_t nom_res = 1.0;

  //set up "no oscillation" histograms:
  TH1F *ndh_noosc[ndelta];
  TH1F *fdh_noosc[ndelta];
  for(int i=0;i<ndelta;i++){
    ndh_noosc[i] = NULL;
    fdh_noosc[i] = NULL;
  }

  //find histograms:
  for(int i=0;i<fdt->GetEntries();i++){
    fdt->GetEntry(i);
    ndt->GetEntry(i);
    int nnbins=(int)(nrw->EBINS/nrw->EBINW);
    int nfbins=(int)(frw->EBINS/frw->EBINW);
    if(fabs(frw->UE32 - nom_UE32)<1e-6    && 
       fabs(frw->kno_r112 - nom_dis)<1e-6 && 
       fabs(frw->qel_ma - nom_qel)<1e-6   &&
       fabs(frw->res_ma - nom_res)<1e-6   && 
       fabs(frw->ss2th - nom_ss2th)<1e-6 ){
      Int_t key = -1;
      for(int j=0;j<ndelta;j++) {
	//cout << frw->dm2 << " " << deltam[j] << endl;
	if(fabs(frw->dm2-deltam[j])<1e-6) key = j;
      }
      //cout << key << endl;
      if(key>=0 && key<ndelta){
	char nom[256];
	sprintf(nom,"ndh_noosc%i",key);
	ndh_noosc[key] = new TH1F(nom,"ndh noosc",nnbins,0,frw->EBINS);
	sprintf(nom,"fdh_noosc%i",key);
	fdh_noosc[key] = new TH1F(nom,"fdh noosc",nfbins,0,nrw->EBINS);
   	for(int j=0;j<frw->EBINS;j++){
	  ndh_noosc[key]->Fill(j,(nrw->nsigE[j]+nrw->nbgE[j])*POT);
	  fdh_noosc[key]->Fill(j,(frw->nsigE[j]+frw->nbgE[j])*POT);
	}
      }
    }
  }

  //loop through again to build sensitivity plot:
  float minchi2=10000.;
  int minindex=0;
  TFile savefile("SensitivityAF.root","RECREATE");
  TH2F *th13vsdm2_FD = new TH2F("ue3vsdm2_FD","Sensitivity Plot",
				ntheta13,ybins,ndelta,xbins);
  TH2F *th13vsdm2_FDSys = new TH2F("ue3vsdm2_FDSys","Sensitivity Plot",
				   ntheta13,ybins,ndelta,xbins);
  TH2F *th13vsdm2_FDND = new TH2F("ue3vsdm2_FDND","Sensitivity Plot",
				  ntheta13,ybins,ndelta,xbins);
  TH1F *th13_FD = new TH1F("ue3_FD","Sensitivity Plot",ntheta13,ybins);
  TH1F *th13_FDSys = new TH1F("ue3_FDSys","Sensitivity Plot",ntheta13,ybins);
  TH1F *th13_FDOscSys = new TH1F("ue3_FDOscSys","Sensitivity Plot",
				 ntheta13,ybins);
  TH1F *th13_FDND = new TH1F("ue3_FDND","Sensitivity Plot",ntheta13,ybins);
  TH1F *th13_ND = new TH1F("ue3_ND","Sensitivity Plot",ntheta13,ybins);

  for(int i=1;i<=ntheta13;i++){
    th13_FD->SetBinContent(i,9999);
    th13_FDSys->SetBinContent(i,9999);
    th13_FDOscSys->SetBinContent(i,9999);
    th13_FDND->SetBinContent(i,9999);
    th13_ND->SetBinContent(i,9999);
    for(int j=1;j<=ndelta;j++){
      th13vsdm2_FD->SetBinContent(i,j,9999);
      th13vsdm2_FDSys->SetBinContent(i,j,9999);
      th13vsdm2_FDND->SetBinContent(i,j,9999);
    }
  }
  
  for(int i=0;i<fdt->GetEntries();i++){
    fdt->GetEntry(i);
    ndt->GetEntry(i);
    
    int nnbins=(int)(nrw->EBINS/nrw->EBINW);
    int nfbins=(int)(frw->EBINS/frw->EBINW);
    
    TH1F ndh("ndh","ndh",nnbins,0,frw->EBINS);
    TH1F fdh("fdh","fdh",nfbins,0,nrw->EBINS);
    for(int j=0;j<frw->EBINS;j++){
      ndh.Fill(j,(nrw->nsigE[j]+nrw->nbgE[j])*POT);
      fdh.Fill(j,(frw->nsigE[j]+frw->nbgE[j])*POT);
    }
    
    TH1F *fdspec = NULL;
    TH1F *ndspec = NULL;
    Int_t key = -1;
    for(int j=0;j<ndelta;j++) {
      if(fabs(frw->dm2-deltam[j])<1e-6) key = j;
    }
    if(key>=0 && key<ndelta) {
      fdspec = fdh_noosc[key];
      ndspec = ndh_noosc[key];
    }
    
    double chi2 = 0;
    int ndf = 0;
    double theta13 = 4.*(frw->UE32)*(1-(frw->UE32));
    Int_t the1Dbin = th13_FD->FindBin(theta13);
    Int_t the2Dbin = th13vsdm2_FD->FindBin(theta13,frw->dm2);

    //Get FD chi2:
    ComputeChi2(fdspec,&fdh,chi2,ndf,1);
    
    if(fabs(frw->kno_r112 - nom_dis) <1e-6  && 
       fabs(frw->qel_ma - nom_qel)   <1e-6  &&
       fabs(frw->res_ma - nom_res)   <1e-6  && 
       fabs(frw->ss2th - nom_ss2th)  <1e-6 ){
      //everything nominal:
      if(chi2<th13vsdm2_FD->GetBinContent(the2Dbin)) {
	th13vsdm2_FD->SetBinContent(the2Dbin,chi2);
      }
      if(fabs(frw->dm2 - nom_dm2)<1e-6) {
	if(chi2<th13_FD->GetBinContent(the1Dbin)){
	  th13_FD->SetBinContent(the1Dbin,chi2);
	}
      }
    }

    //allow oscillation systematics
    if(fabs(frw->kno_r112 - nom_dis) <1e-6  && 
       fabs(frw->qel_ma - nom_qel)   <1e-6  &&
       fabs(frw->res_ma - nom_res)   <1e-6) {
      if(chi2<th13_FDOscSys->GetBinContent(the1Dbin)){
	th13_FDOscSys->SetBinContent(the1Dbin,chi2);
      }
    }

    //allow systematics:
    if(chi2<th13vsdm2_FDSys->GetBinContent(the2Dbin)) {
      th13vsdm2_FDSys->SetBinContent(the2Dbin,chi2);
    }
    if(chi2<th13_FDSys->GetBinContent(the1Dbin)){
      th13_FDSys->SetBinContent(the1Dbin,chi2);
    }

    //Add on ND chi2:
    ComputeChi2(ndspec,&ndh,chi2,ndf,1);

    //see how much ND constrains things:
    if(chi2<th13vsdm2_FDND->GetBinContent(the2Dbin)) {
      th13vsdm2_FDND->SetBinContent(the2Dbin,chi2);
    }
    if(chi2<th13_FDND->GetBinContent(the1Dbin)){
      th13_FDND->SetBinContent(the1Dbin,chi2);
    }

    //Just ND chi2:
    chi2 = 0; ndf = 0;
    ComputeChi2(ndspec,&ndh,chi2,ndf,1);
    if(chi2<th13_ND->GetBinContent(the1Dbin)) {
      th13_ND->SetBinContent(the1Dbin,chi2);
    }

    if(chi2<minchi2){
      cout<<"New min: "<<chi2<<" i "<<i<<endl;
      minchi2=chi2;
      minindex=i;
    }
  }

  fdt->GetEntry(minindex);
  cout<<"minchi2 "<<minchi2<<" index "<<minindex<<endl;
  cout<<" dm2 "<<frw->dm2<<" thetha "<< frw->ss2th
      <<" qelma "<<frw->qel_ma<<" res "<<frw->res_ma
      <<" kno "<<frw->kno_r112<<" UE32 "<<frw->UE32<<endl;

  gROOT->LoadMacro("NueAna/macros/chooz.C");
  TCanvas *can1 = new TCanvas();
  TCanvas *can2 = new TCanvas();
  TCanvas *can3 = new TCanvas();
  TCanvas *can4 = new TCanvas();
  TCanvas *can5 = new TCanvas();
  TCanvas *can6 = new TCanvas();
  //can->Divide(2,3);
  can1->cd();
  th13vsdm2_FD->Draw("lego2");
  can2->cd();
  TH2F *th13vsdm2_FD_copy = new TH2F(*th13vsdm2_FD);
  double lev[1] = {2.71};
  lev[0]+=th13vsdm2_FD->GetMinimum();
  th13vsdm2_FD_copy->SetContour(1,lev);
  th13vsdm2_FD_copy->SetTitle("90% Sensitivity Contours for 7.4e20POT, #nu_{#mu} Oscillation Parameters, FD Only");
  th13vsdm2_FD_copy->SetXTitle("sin^{2}(2#theta_{13})");
  th13vsdm2_FD_copy->SetYTitle("#Deltam^{2}_{23} (eV^{2})");
  th13vsdm2_FD_copy->Draw("cont2");
  ChoozDraw();
  TGraph *chooz = (TGraph*) gROOT->FindObject("chooz");

  can3->cd();
  th13vsdm2_FDSys->Draw("lego2");
  can4->cd();
  TH2F *th13vsdm2_FDSys_copy = new TH2F(*th13vsdm2_FDSys);
  lev[0]=4.61+th13vsdm2_FDSys->GetMinimum();
  th13vsdm2_FDSys_copy->SetContour(1,lev);
  th13vsdm2_FDSys_copy->SetTitle("90% Sensitivity Contours for 7.4e20POT, 1-sigma systematics on Neugen Parameters, #nu_{#mu} Oscillation Parameters, FD Only");
  th13vsdm2_FDSys_copy->SetXTitle("sin^{2}(2#theta_{13})");
  th13vsdm2_FDSys_copy->SetYTitle("#Deltam^{2}_{23} (eV^{2})");
  th13vsdm2_FDSys_copy->Draw("cont2");
  chooz->Draw("L");

  can5->cd();
  th13vsdm2_FDND->Draw("lego2");
  can6->cd();
  TH2F *th13vsdm2_FDND_copy = new TH2F(*th13vsdm2_FDND);
  lev[0]=4.61+th13vsdm2_FDND->GetMinimum();
  th13vsdm2_FDND_copy->SetContour(1,lev);
  th13vsdm2_FDND_copy->SetTitle("90% Sensitivity Contours for 7.4e20POT, 1-sigma systematics on Neugen Parameters, #nu_{#mu} Oscillation Parameters, ND+FD");
  th13vsdm2_FDND_copy->SetXTitle("sin^{2}(2#theta_{13})");
  th13vsdm2_FDND_copy->SetYTitle("#Deltam^{2}_{23} (eV^{2})");
  th13vsdm2_FDND_copy->Draw("cont2");
  chooz->Draw("L");
  
  TLegend *leg = new TLegend(0.2,0.2,0.5,0.5);
  leg->AddEntry(th13vsdm2_FD_copy,"TJ Analysis","l");
  leg->AddEntry(chooz,"CHOOZ","l");
  leg->Draw();

  TCanvas *can7 = new TCanvas();
  can7->cd();
  th13_FD->SetLineColor(2); th13_FD->Draw();
  th13_FDSys->SetLineColor(3); th13_FDSys->Draw("sames");
  th13_FDOscSys->SetLineColor(7); th13_FDOscSys->Draw("sames");
  th13_FDND->SetLineColor(4); th13_FDND->Draw("sames");
  th13_ND->SetLineColor(6); th13_ND->Draw("sames");
  savefile.Write();
  savefile.Close();

}


void StudyND(const char *rwfile,Bool_t MDC)
{
  Int_t nominal_entry = 2456;  //SpecialND
  //nominal_entry = 2541;      //normal
  Int_t nfiles = 1;

  TFile *g = new TFile(rwfile);
  TTree *ndt = NULL;
  Double_t scaleFact = 1.0;
  if(MDC){
    ndt=(TTree *)g->Get("rwtree_ND_norm");
    scaleFact =7.4;
  }
  else {
    ndt=(TTree *)g->Get("rwtree_ND");    
    scaleFact = 1.0;
  }

  if(ndt->GetEntries()==0){
    cout<<"ND TREE size = 0.  Won't continue"<<endl;
    return;
  }


  Int_t NQELBINS = 17;
  Float_t FIRSTQEL = 0.956;
  Float_t QELWIDTH = 0.0055;
  Int_t nomQELkey = Int_t(Float_t(NQELBINS)/2. - 0.5);
  Int_t NRESBINS = 17;
  Float_t FIRSTRES = 0.956;
  Float_t RESWIDTH = 0.0055;
  Int_t nomRESkey = Int_t(Float_t(NRESBINS)/2. - 0.5);
  Int_t NDISBINS = 17;
  Float_t FIRSTDIS = 0.994;
  Float_t DISWIDTH = 0.00075;
  Int_t nomDISkey = Int_t(Float_t(NDISBINS)/2. - 0.5);


  /*   
  Int_t NQELBINS = 17;
  Float_t FIRSTQEL = 0.912;
  Float_t QELWIDTH = 0.011;
  Int_t nomQELkey = Int_t(Float_t(NQELBINS)/2. - 0.5);
  Int_t NRESBINS = 17;
  Float_t FIRSTRES = 0.912;
  Float_t RESWIDTH = 0.011;
  Int_t nomRESkey = Int_t(Float_t(NRESBINS)/2. - 0.5);
  Int_t NDISBINS = 17;
  Float_t FIRSTDIS = 0.88;
  Float_t DISWIDTH = 0.015;
  Int_t nomDISkey = Int_t(Float_t(NDISBINS)/2. - 0.5);
  */

  /*  
  Int_t NQELBINS = 3;
  Float_t FIRSTQEL = 0.97;
  Float_t QELWIDTH = 0.03;
  Int_t nomQELkey = Int_t(Float_t(NQELBINS)/2. - 0.5);
  Int_t NRESBINS = 3;
  Float_t FIRSTRES = 0.97;
  Float_t RESWIDTH = 0.03;
  Int_t nomRESkey = Int_t(Float_t(NRESBINS)/2. - 0.5);
  Int_t NDISBINS = 3;
  Float_t FIRSTDIS = 0.96;
  Float_t DISWIDTH = 0.04;
  Int_t nomDISkey = Int_t(Float_t(NDISBINS)/2. - 0.5);
  */

  //Plan for many TGraphs:
  float *qel_ma = new float[NQELBINS];
  float *res_ma = new float[NRESBINS];
  float *dis_fac = new float[NDISBINS];
  float *chi2_qel_ma = new float[NQELBINS];
  float *chi2_res_ma = new float[NRESBINS];
  float *chi2_dis_fac = new float[NDISBINS];
  
  for(int i=0;i<NQELBINS;i++){
    qel_ma[i] = FIRSTQEL + Float_t(i)*QELWIDTH;
    chi2_qel_ma[i] = 0;
  }
  for(int i=0;i<NRESBINS;i++){
    res_ma[i] = FIRSTRES + Float_t(i)*RESWIDTH;
    chi2_res_ma[i] = 0;
  } 
  for(int i=0;i<NDISBINS;i++){
    dis_fac[i] = FIRSTDIS + Float_t(i)*DISWIDTH;  
    chi2_dis_fac[i] = 0;
  }
  
  TH2F *qelVsres = new TH2F("qelVsres",
			    "#Delta#chi^{2} vs qel_ma and res_ma at nominal dis_fac",
			    NRESBINS,FIRSTRES - 0.5*RESWIDTH,
			    FIRSTRES + Float_t(NRESBINS-0.5)*RESWIDTH,
			    NQELBINS,FIRSTQEL - 0.5*QELWIDTH,
			    FIRSTQEL + Float_t(NQELBINS-0.5)*QELWIDTH);
  qelVsres->SetXTitle("res_ma/res_ma_{nominal}");
  qelVsres->SetYTitle("qel_ma/qel_ma_{nominal}");
  TH2F *qelVsdis = new TH2F("qelVsdis",
			    "#Delta#chi^{2} vs qel_ma and dis_fac at nominal res_ma",
			    NDISBINS,FIRSTDIS - 0.5*DISWIDTH,
			    FIRSTDIS + Float_t(NDISBINS-0.5)*DISWIDTH,
			    NQELBINS,FIRSTQEL - 0.5*QELWIDTH,
			    FIRSTQEL + Float_t(NQELBINS-0.5)*QELWIDTH);
  qelVsdis->SetXTitle("dis_ma/dis_ma_{nominal}");
  qelVsdis->SetYTitle("qel_ma/qel_ma_{nominal}");
  TH2F *disVsres = new TH2F("disVsres",
			    "#Delta#chi^{2} vs dis_fac and res_ma at nominal qel_ma",
			    NRESBINS,FIRSTRES - 0.5*RESWIDTH,
			    FIRSTRES + Float_t(NRESBINS-0.5)*RESWIDTH,
			    NDISBINS,FIRSTDIS - 0.5*DISWIDTH,
			    FIRSTDIS + Float_t(NDISBINS-0.5)*DISWIDTH);
  disVsres->SetXTitle("res_ma/res_ma_{nominal}");
  disVsres->SetYTitle("dis_fac/dis_fac_{nominal}");

  NueRW *nrw = new NueRW();
  ndt->SetBranchAddress("NueRW",&nrw);

  ndt->GetEntry(nominal_entry);
  int nnbins=(int)(nrw->EBINS/nrw->EBINW);
  TH1F *nominal = new TH1F("nominal","nominal",nnbins,0,nrw->EBINS);
  double nom_ue32 = nrw->UE32;
  double nom_ss2th = nrw->ss2th;
  double nom_dm2 = nrw->dm2;
  for(int j=0;j<nrw->EBINS;j++){
    nominal->Fill(j,nrw->nsigE[j]+nrw->nbgE[j]);
  }
  nominal->Scale(scaleFact);

  TH1F *worstfit[3];
  worstfit[0] = NULL; worstfit[1] = NULL; worstfit[2] = NULL;
  Int_t worstline = 0;
  for(int i=0;i<ndt->GetEntries();i++){
    if(i==nominal_entry) continue;
    ndt->GetEntry(i);
    if(nrw->UE32!=nom_ue32 || nrw->ss2th!=nom_ss2th || nrw->dm2!=nom_dm2) continue;

    nfiles = nrw->nfiles;

    Int_t QELkey = 0;
    Int_t RESkey = 0;
    Int_t DISkey = 0;
    for(int j=0;j<NQELBINS;j++) {
      if(fabs(nrw->qel_ma-qel_ma[j])<1e-6) QELkey = j;
    }
    for(int j=0;j<NRESBINS;j++) {
      if(fabs(nrw->res_ma-res_ma[j])<1e-6) RESkey = j;
    }
    for(int j=0;j<NDISBINS;j++) {
      if(fabs(nrw->kno_r112-dis_fac[j])<1e-6) DISkey = j;
    }
    
    TH1F ndh("ndh","ndh",nnbins,0,nrw->EBINS);
    for(int j=0;j<nrw->EBINS;j++){
      ndh.Fill(j,nrw->nsigE[j]+nrw->nbgE[j]);
    }
    ndh.Scale(scaleFact);
    
    double chi2 = 0;
    int ndf = 0;
    ComputeChi2(nominal,&ndh,chi2,ndf,1);
    if(DISkey==nomDISkey) {
      qelVsres->Fill(res_ma[RESkey],qel_ma[QELkey],chi2);
      cout << i<< " " << DISkey << " " << RESkey << " " << QELkey << " "
	   << dis_fac[DISkey] << " " << res_ma[RESkey] << " "
	   << qel_ma[QELkey] << endl;
    }
    if(RESkey==nomRESkey) qelVsdis->Fill(dis_fac[DISkey],qel_ma[QELkey],chi2);
    if(QELkey==nomQELkey) disVsres->Fill(res_ma[RESkey],dis_fac[DISkey],chi2);

    if(chi2>chi2_qel_ma[QELkey]){
      chi2_qel_ma[QELkey] = chi2;
      if(worstfit[0]) delete worstfit[0];
      worstfit[0] = new TH1F(ndh);
      worstfit[0]->SetName("qel");
    }
    if(chi2>chi2_res_ma[RESkey]){
      chi2_res_ma[RESkey] = chi2;
      if(worstfit[1]) delete worstfit[1];
      worstfit[1] = new TH1F(ndh);
      worstfit[1]->SetName("res");
    }
    if(chi2>chi2_dis_fac[DISkey]){
      chi2_dis_fac[DISkey] = chi2;
      if(worstfit[2]) delete worstfit[2];
      worstfit[2] = new TH1F(ndh);      
      worstfit[2]->SetName("dis");
      worstline = i;
    }
  }
  
  TGraph *graph_qel_ma = new TGraph(NQELBINS,qel_ma,chi2_qel_ma);
  graph_qel_ma->SetTitle("qel_ma");
  graph_qel_ma->GetXaxis()->SetTitle("qel_ma/qel_ma_{nominal}");
  graph_qel_ma->GetYaxis()->SetTitle("#chi^{2}");
  TGraph *graph_res_ma = new TGraph(NQELBINS,res_ma,chi2_res_ma);
  graph_res_ma->SetTitle("res_ma");
  graph_res_ma->GetXaxis()->SetTitle("res_ma/res_ma_{nominal}");
  graph_res_ma->GetYaxis()->SetTitle("#chi^{2}");
  TGraph *graph_dis_fac = new TGraph(NQELBINS,dis_fac,chi2_dis_fac);
  graph_dis_fac->SetTitle("dis_fac");
  graph_dis_fac->GetXaxis()->SetTitle("dis_fac/dis_fac_{nominal}");
  graph_dis_fac->GetYaxis()->SetTitle("#chi^{2}");
  TCanvas *can = new TCanvas();
  can->Divide(1,3);
  can->cd(1);
  graph_qel_ma->Draw("ALP");
  can->cd(2);
  graph_res_ma->Draw("ALP");
  can->cd(3);
  graph_dis_fac->Draw("ALP");

  TCanvas *can2a = new TCanvas();
  TCanvas *can2b = new TCanvas();
  TCanvas *can2c = new TCanvas();
  //can2->Divide(1,3);
  can2a->cd();
  qelVsres->Draw("COLZ");
  can2b->cd();
  qelVsdis->Draw("COLZ");
  can2c->cd();
  disVsres->Draw("COLZ");

  TCanvas *can3a = new TCanvas();
  TCanvas *can3b = new TCanvas();
  TCanvas *can3c = new TCanvas();
  //can3->Divide(1,3);
  TH1F *dummy90 = new TH1F("dummy90","dummy90",1,0,1);
  dummy90->SetLineColor(3);
  TH1F *dummy95 = new TH1F("dummy95","dummy95",1,0,1);
  dummy95->SetLineColor(2);
  TLegend *leg = new TLegend(0.8,0.7,0.9,0.85,"");
  leg->AddEntry(dummy90,"90% CL","L");
  leg->AddEntry(dummy95,"95% CL","L");
  
  double lev[2] = {4.61,6.18};
  can3a->cd();
  TH2F *qelVsres_copy = new TH2F(*qelVsres);
  qelVsres_copy->SetContour(2,lev);
  qelVsres_copy->SetLineColor(4);
  qelVsres_copy->Draw("CONT1");
  leg->Draw();
  can3b->cd();
  TH2F *qelVsdis_copy = new TH2F(*qelVsdis);
  qelVsdis_copy->SetContour(2,lev);
  qelVsdis_copy->SetLineColor(4);
  qelVsdis_copy->Draw("CONT1");
  leg->Draw();
  can3c->cd();
  TH2F *disVsres_copy = new TH2F(*disVsres);
  disVsres_copy->SetContour(2,lev);
  disVsres_copy->SetLineColor(4);
  disVsres_copy->Draw("CONT1");
  leg->Draw();

  TCanvas *can4 = new TCanvas();
  //can4->Divide(1,3);
  gStyle->SetOptStat(1110);
  nominal->SetLineColor(4);
  can4->cd(1);
  nominal->SetTitle("Near Detector #nu_{e} Background Energy Spectrum");
  nominal->SetXTitle("Energy (GeV)");
  nominal->Draw("e1");
  worstfit[2]->Draw("sames");
  TLegend *leg2 = new TLegend(0.5,0.5,0.7,0.7,"");
  leg2->AddEntry(nominal,"Nominal Spectrum","lp");
  leg2->AddEntry(worstfit[2],"'Worst Fit', 1-#sigma Neugen Sys.","l");
  leg2->SetBorderSize(0);
  leg2->Draw();
  char text[256];
  sprintf(text,"#Delta#chi^{2}=%f",chi2_dis_fac[0]);
  TLatex *tex = new TLatex(0.3,0.7,text);
  tex->SetNDC();
  tex->Draw();

  double thePOT = 550*2.4e13*nfiles;
  if(scaleFact>1.1) thePOT = 7.4e20;
  sprintf(text,"POT = %fe20",thePOT/1e20);
  TLatex *tex2 = new TLatex(0.3,0.5,text);
  tex2->SetNDC();
  tex2->Draw();

  //can4->cd(2);
  //nominal->Draw("e1");
  //worstfit[1]->Draw("sames");
  //can4->cd(3);
  //nominal->Draw("e1");
  //worstfit[2]->Draw("sames");

  TTree *fdt = NULL;
  fdt=(TTree *)g->Get("rwtree_FD_norm");
  NueRW *frw = new NueRW();
  fdt->SetBranchAddress("NueRW",&frw);  

  fdt->GetEntry(nominal_entry);
  TH1F *nom_fdh = new TH1F("nom_fdh","Nominal FD",nnbins,0,frw->EBINS);
  for(int j=0;j<frw->EBINS;j++){
    nom_fdh->Fill(j,frw->nsigE[j]+frw->nbgE[j]);
  }
  nom_fdh->Scale(scaleFact);
  
  fdt->GetEntry(worstline);
  TH1F *fdh = new TH1F("fdh","FD",nnbins,0,frw->EBINS);
  for(int j=0;j<frw->EBINS;j++){
    fdh->Fill(j,frw->nsigE[j]+frw->nbgE[j]);
  }
  fdh->Scale(scaleFact);
  
  new TCanvas();
  nom_fdh->SetLineColor(4);
  nom_fdh->Draw("e1");
  fdh->Draw("sames");
  leg2->Draw();
  Double_t fdchi2 = 0;
  Int_t fdndf = 0;
  ComputeChi2(nom_fdh,fdh,fdchi2,fdndf,1);
  sprintf(text,"#Delta#chi^{2}=%f",fdchi2);
  TLatex *tex3 = new TLatex(0.3,0.7,text);
  tex3->SetNDC();
  tex3->Draw();

  sprintf(text,"POT = 22.2e20");
  TLatex *tex4 = new TLatex(0.3,0.5,text);
  tex4->SetNDC();
  tex4->Draw();
  return;
}


void ComputeChi2(TH1 *h1, TH1 *h2, double &chi2, int &ndof,int method)
{

  if(!h1 || !h2) {
    //cout<<"One or other histograms missing, won't do chi2"<<endl;
    chi2=-1.;
    ndof=0;
    return;
  }
  if(h1->GetNbinsX() != h2->GetNbinsX()) {
    //cout<<"Number of bins don't match, won't do chi2"<<endl;
    chi2=-1.;
    ndof=0;
    return;
  }

  for(int i=1;i<=h1->GetNbinsX();i++){
    
    double h1bc=h1->GetBinContent(i);
    double h2bc=h2->GetBinContent(i);
    
    if(method==0){
      double num = pow(h1bc-h2bc,2);
      double denom = h2bc;
      if(denom!=0){
	chi2+=num/denom;
	ndof++;
      }
    }
    else if(method==1){
      if(h1bc==0 || h2bc==0) chi2 += 2.*(h1bc - h2bc);
      else chi2 += 2.*(h1bc - h2bc) + 2.*h2bc*TMath::Log(h2bc/h1bc);
      ndof++;
    }
  }
  return;
}

void FindUpperBounds(TGraph *g, float &b68, float &b1, float &b2)
{
   b68=-1.;
   b1=-1.;
   b2=-1.;
   for(int i=1;i<g->GetN();i++){
      if((g->GetY()[i]>1&&g->GetY()[i-1]<1)){
	 float m=(g->GetY()[i]-g->GetY()[i-1])/(g->GetX()[i]-g->GetX()[i-1]);
	 b68=g->GetX()[i-1]+(1.-g->GetY()[i-1])/m;
      }
      if((g->GetY()[i]>2.71&&g->GetY()[i-1]<2.71)){
	 float m=(g->GetY()[i]-g->GetY()[i-1])/(g->GetX()[i]-g->GetX()[i-1]);
	 b1=g->GetX()[i-1]+(2.71-g->GetY()[i-1])/m;
      }
      if((g->GetY()[i]>6.63&&g->GetY()[i-1]<6.63)){
	 float m=(g->GetY()[i]-g->GetY()[i-1])/(g->GetX()[i]-g->GetX()[i-1]);
	 b2=g->GetX()[i-1]+(6.63-g->GetY()[i-1])/m;       
      }
   }

}


void FindLowerBounds(TGraph *g, float &b68, float &b1, float &b2)
{
   b68=-1.;
   b1=-1.;
   b2=-1.;
   for(int i=1;i<g->GetN();i++){
      if(g->GetY()[i]<1&&g->GetY()[i-1]>1){
	 float m=(g->GetY()[i]-g->GetY()[i-1])/(g->GetX()[i]-g->GetX()[i-1]);
	 b68=g->GetX()[i-1]+(1.-g->GetY()[i-1])/m;
      }
      if(g->GetY()[i]<2.71&&g->GetY()[i-1]>2.71){
	 float m=(g->GetY()[i]-g->GetY()[i-1])/(g->GetX()[i]-g->GetX()[i-1]);
	 b1=g->GetX()[i-1]+(2.71-g->GetY()[i-1])/m;
      }
      if(g->GetY()[i]<6.63&&g->GetY()[i-1]>6.63){
	 float m=(g->GetY()[i]-g->GetY()[i-1])/(g->GetX()[i]-g->GetX()[i-1]);
	 b2=g->GetX()[i-1]+(6.63-g->GetY()[i-1])/m;       
      }
   }

}




void PlotSummary()
{
   float tj[4]={0,0.0396,0.0199,0.0396};
   float tj_l90[4]={-1,-1,-1,-1};
   float tj_u90[4]={0.170783,0.136346,0.152249,0.121706};
   float tj_l99[4]={-1,-1,-1,-1};
   float tj_u99[4]={0.284566,0.20688,0.246099,0.180711};

   float jbms[4]={0.1164,0.1164,0.1351,0.1164};
   float jbms_l90[4]={-1,0.00376866,-1,0.0470831};
   float jbms_u90[4]={0.317089,0.233849,0.309662,0.236962};
   float jbms_l99[4]={-1,-1,-1,0.00539893};
   float jbms_u99[4]={0.450445,0.303079,0.416635,0.295629};

   float as[4]={0.0784,0.0784,0.0784,0.0975};
   float as_l90[4]={-1,0.0233034,-1,0.0330758};
   float as_u90[4]={0.246546,0.207458,0.218101,0.19038};
   float as_l99[4]={-1,-1,-1,0.0114224};
   float as_u99[4]={0.357107,0.270169,0.310813,0.241444};


   TGraphAsymmErrors *ls_nm_ans90=new TGraphAsymmErrors();
   ls_nm_ans90->SetPoint(0,tj[0],1);
   ls_nm_ans90->SetPointError(0,tj[0]-tj_l90[0],tj_u90[0]-tj[0],0,0);
   ls_nm_ans90->SetPoint(1,jbms[0],2);
   ls_nm_ans90->SetPointError(1,jbms[0]-jbms_l90[0],jbms_u90[0]-jbms[0],0,0);
   ls_nm_ans90->SetPoint(2,as[0],3);
   ls_nm_ans90->SetPointError(2,as[0]-as_l90[0],as_u90[0]-as[0],0,0);
   ls_nm_ans90->SetNameTitle("ls_nm_ans90",
			     "No Matter Effect, 7.4x10^{20} POT;sin^{2}(2#theta_{13})");
   
   TGraphAsymmErrors *ls_nm_ans99=new TGraphAsymmErrors();
   ls_nm_ans99->SetPoint(0,tj[0],1);
   ls_nm_ans99->SetPointError(0,tj[0]-tj_l99[0],tj_u99[0]-tj[0],0,0);
   ls_nm_ans99->SetPoint(1,jbms[0],2);
   ls_nm_ans99->SetPointError(1,jbms[0]-jbms_l99[0],jbms_u99[0]-jbms[0],0,0);
   ls_nm_ans99->SetPoint(2,as[0],3);
   ls_nm_ans99->SetPointError(2,as[0]-as_l99[0],as_u99[0]-as[0],0,0);
   ls_nm_ans99->SetNameTitle("ls_nm_ans99",
			     "No Matter Effect, 7.4x10^{20} POT;sin^{2}(2#theta_{13})");
   
   TGraphAsymmErrors *hs_nm_ans90=new TGraphAsymmErrors();
   hs_nm_ans90->SetPoint(0,tj[1],1);
   hs_nm_ans90->SetPointError(0,tj[1]-tj_l90[1],tj_u90[1]-tj[1],0,0);
   hs_nm_ans90->SetPoint(1,jbms[1],2);
   hs_nm_ans90->SetPointError(1,jbms[1]-jbms_l90[1],jbms_u90[1]-jbms[1],0,0);
   hs_nm_ans90->SetPoint(2,as[1],3);
   hs_nm_ans90->SetPointError(2,as[1]-as_l90[1],as_u90[1]-as[1],0,0);
   hs_nm_ans90->SetNameTitle("hs_nm_ans90",
			     "No Matter Effect, 22.2x10^{20} POT;sin^{2}(2#theta_{13})");
   
   TGraphAsymmErrors *hs_nm_ans99=new TGraphAsymmErrors();
   hs_nm_ans99->SetPoint(0,tj[1],1);
   hs_nm_ans99->SetPointError(0,tj[1]-tj_l99[1],tj_u99[1]-tj[1],0,0);
   hs_nm_ans99->SetPoint(1,jbms[1],2);
   hs_nm_ans99->SetPointError(1,jbms[1]-jbms_l99[1],jbms_u99[1]-jbms[1],0,0);
   hs_nm_ans99->SetPoint(2,as[1],3);
   hs_nm_ans99->SetPointError(2,as[1]-as_l99[1],as_u99[1]-as[1],0,0);
   hs_nm_ans99->SetNameTitle("hs_nm_ans99",
			     "No Matter Effect, 22.2x10^{20} POT;sin^{2}(2#theta_{13})");
   
   
   TGraphAsymmErrors *ls_wm_ans90=new TGraphAsymmErrors();
   ls_wm_ans90->SetPoint(0,tj[2],1);
   ls_wm_ans90->SetPointError(0,tj[2]-tj_l90[2],tj_u90[2]-tj[2],0,0);
   ls_wm_ans90->SetPoint(1,jbms[2],2);
   ls_wm_ans90->SetPointError(1,jbms[2]-jbms_l90[2],jbms_u90[2]-jbms[2],0,0);
   ls_wm_ans90->SetPoint(2,as[2],3);
   ls_wm_ans90->SetPointError(2,as[2]-as_l90[2],as_u90[2]-as[2],0,0);
   ls_wm_ans90->SetNameTitle("ls_wm_ans90",
			     "With Matter Effect, 7.4x10^{20} POT;sin^{2}(2#theta_{13})");
   
   TGraphAsymmErrors *ls_wm_ans99=new TGraphAsymmErrors();
   ls_wm_ans99->SetPoint(0,tj[2],1);
   ls_wm_ans99->SetPointError(0,tj[2]-tj_l99[2],tj_u99[2]-tj[2],0,0);
   ls_wm_ans99->SetPoint(1,jbms[2],2);
   ls_wm_ans99->SetPointError(1,jbms[2]-jbms_l99[2],jbms_u99[2]-jbms[2],0,0);
   ls_wm_ans99->SetPoint(2,as[2],3);
   ls_wm_ans99->SetPointError(2,as[2]-as_l99[2],as_u99[2]-as[2],0,0);
   ls_wm_ans99->SetNameTitle("ls_wm_ans99",
			     "With Matter Effect, 7.4x10^{20} POT;sin^{2}(2#theta_{13})");
   
   TGraphAsymmErrors *hs_wm_ans90=new TGraphAsymmErrors();
   hs_wm_ans90->SetPoint(0,tj[3],1);
   hs_wm_ans90->SetPointError(0,tj[3]-tj_l90[3],tj_u90[3]-tj[3],0,0);
   hs_wm_ans90->SetPoint(1,jbms[3],2);
   hs_wm_ans90->SetPointError(1,jbms[3]-jbms_l90[3],jbms_u90[3]-jbms[3],0,0);
   hs_wm_ans90->SetPoint(2,as[3],3);
   hs_wm_ans90->SetPointError(2,as[3]-as_l90[3],as_u90[3]-as[3],0,0);
   hs_wm_ans90->SetNameTitle("hs_wm_ans90",
			     "With Matter Effect, 22.2x10^{20} POT;sin^{2}(2#theta_{13})");
   
   TGraphAsymmErrors *hs_wm_ans99=new TGraphAsymmErrors();
   hs_wm_ans99->SetPoint(0,tj[3],1);
   hs_wm_ans99->SetPointError(0,tj[3]-tj_l99[3],tj_u99[3]-tj[3],0,0);
   hs_wm_ans99->SetPoint(1,jbms[3],2);
   hs_wm_ans99->SetPointError(1,jbms[3]-jbms_l99[3],jbms_u99[3]-jbms[3],0,0);
   hs_wm_ans99->SetPoint(2,as[3],3);
   hs_wm_ans99->SetPointError(2,as[3]-as_l99[3],as_u99[3]-as[3],0,0);
   hs_wm_ans99->SetNameTitle("hs_wm_ans99",
			     "With Matter Effect, 22.2x10^{20} POT;sin^{2}(2#theta_{13})");
   
      
   TGraph *stupid=new TGraph();
   stupid->SetPoint(0,0.00001,0);
   stupid->SetPoint(1,0.25,4);

   stupid->SetMarkerColor(0);
   stupid->SetMarkerSize(0.001);
   stupid->SetMarkerStyle(1);
   ls_nm_ans90->SetLineColor(1);
   ls_nm_ans99->SetLineColor(2);
   hs_nm_ans90->SetLineColor(1);
   hs_nm_ans99->SetLineColor(2);
   ls_wm_ans90->SetLineColor(1);
   ls_wm_ans99->SetLineColor(2);
   hs_wm_ans90->SetLineColor(1);
   hs_wm_ans99->SetLineColor(2);

   TMultiGraph *mlsnm = new TMultiGraph("mlsnm",
					"No Matter Effect, 7.4x10^{20} POT;sin^{2}(2#theta_{13})");
   mlsnm->Add(ls_nm_ans99);
   mlsnm->Add(ls_nm_ans90);
   mlsnm->Add(stupid);

   TMultiGraph *mhsnm = new TMultiGraph("mhsnm",
					"No Matter Effect, 22.2x10^{20} POT;sin^{2}(2#theta_{13})");
   mhsnm->Add(hs_nm_ans99);
   mhsnm->Add(hs_nm_ans90);
   mhsnm->Add(stupid);

   TMultiGraph *mlswm = new TMultiGraph("mlswm",
					"With Matter Effect, 7.4x10^{20} POT;sin^{2}(2#theta_{13})");
   mlswm->Add(ls_wm_ans99);
   mlswm->Add(ls_wm_ans90);
   mlswm->Add(stupid);

   TMultiGraph *mhswm = new TMultiGraph("mhswm",
					"With Matter Effect, 22.2x10^{20} POT;sin^{2}(2#theta_{13})");
   mhswm->Add(hs_wm_ans99);
   mhswm->Add(hs_wm_ans90);
   mhswm->Add(stupid);

//   float choozelimit=.21982;
   float choozelimit=.185;
   TLine *cl = new TLine(choozelimit,0,choozelimit,4);
   cl->SetLineColor(4);
   cl->SetLineWidth(2);
   cl->SetLineStyle(2);
   
   char *names[3]={"TJ","JBMS","AS"};
   TText *t1=new TText(0.0396,1.1,names[0]);
   TText *t2=new TText(0.0396,2.1,names[1]);
   TText *t3=new TText(0.0396,3.1,names[2]);

   TCanvas *fa = new TCanvas("fa","fa");
   fa->Divide(2,2);
   fa->cd(1);
   mlsnm->Draw("AP");
   mlsnm->GetXaxis()->CenterTitle();
   mlsnm->GetXaxis()->SetRangeUser(0.0001,0.3);
   cl->Draw();
   t1->Draw();
   t2->Draw();
   t3->Draw();
   fa->cd(2);
   mhsnm->Draw("AP");
   mhsnm->GetXaxis()->CenterTitle();
   mhsnm->GetXaxis()->SetRangeUser(0.0001,0.3);
   cl->Draw();
   t1->Draw();
   t2->Draw();
   t3->Draw();
   fa->cd(3);
   mlswm->Draw("AP");
   mlswm->GetXaxis()->CenterTitle();
   mlswm->GetXaxis()->SetRangeUser(0.0001,0.3);
   cl->Draw();
   t1->Draw();
   t2->Draw();
   t3->Draw();
   fa->cd(4);
   mhswm->Draw("AP");
   mhswm->GetXaxis()->CenterTitle();
   mhswm->GetXaxis()->SetRangeUser(0.0001,0.3);
   cl->Draw();
   t1->Draw();
   t2->Draw();
   t3->Draw();

}

