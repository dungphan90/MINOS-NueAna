#include "TMLPAnalyzer.h"
#include "TMultiLayerPerceptron.h"
#include "TFile.h"
#include "TCanvas.h"

void TMLP()
{
	TFile * fin = TFile::Open("MLP.root");
	TMultiLayerPerceptron * mlp=0;
	mlp=(TMultiLayerPerceptron*)fin->Get("MLP");
	
	if(!mlp)return;
	
	  TCanvas* mlpa_canvas = new TCanvas("mlpa_canvas","Network analysis");
   mlpa_canvas->Divide(2,2);
   TMLPAnalyzer ana(mlp);
   // Initialisation
   ana.GatherInformations();
   // output to the console
   ana.CheckNetwork();
   mlpa_canvas->cd(1);
   // shows how each variable influences the network
   ana.DrawDInputs();
   mlpa_canvas->cd(2);
   // shows the network structure
   mlp->Draw();
   mlpa_canvas->cd(3);
   // draws the resulting network
   ana.DrawNetwork(0,"type==1","type==0");

	
}