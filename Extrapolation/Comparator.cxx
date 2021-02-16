#include "NueAna/Extrapolation/Comparator.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TMath.h"

using namespace std;

Comparator::Comparator(string whichHist) :
  fNDDataPOT(0),fFDDataPOT(0),fHistType(whichHist),fDoPrint(0)
{  
  sprintf(fSelection,"Unknown");
  fOscSysString = "Delta23";

  Int_t colourArray[20] = {2,3,4,6,7,8,9,11,13,28,29,30,34,36,38,42,46,48,49,50};
  for(int i=0;i<20;i++) fColourArray[i] = colourArray[i];
  
  fCanvasForPredictions = new TCanvas("CanvasForPredictions","Canvas For Predictions",
				      200,50,550,400);
  fCanvasForPredictions->SetGridx(); fCanvasForPredictions->SetGridy();
  fCanvasForRatios = new TCanvas("CanvasForRatios","Canvas For Ratios",
				 200,500,550,400);
  fCanvasForRatios->SetGridx(); fCanvasForRatios->SetGridy();
  fCanvasForIntegrals = new TCanvas("CanvasForIntegrals","Canvas For Integrals",
				    800,50,550,400);
  fCanvasForIntegrals->SetGridx(); fCanvasForIntegrals->SetGridy();
  fCanvasForSummary = new TCanvas("CanvasForSummary","Canvas For Summary",
				  800,500,550,400);
  fCanvasForSummary->SetGridx(); fCanvasForSummary->SetGridy();
  
  fPredictionsForDrawing = new THStack("PredictionsForDrawing","Predictions");
  fRatiosForDrawing = new THStack("RatiosForDrawing","Ratios");
  fIntegralsForDrawing = new TMultiGraph("IntegralsForDrawing","Integrals");
  fSummaryForDrawing = new TMultiGraph("SummaryForDrawing","Summary");
}

Comparator::~Comparator()
{
  //delete predictions:
  std::map<Background::Background_t,NueBackground*>::iterator dataBeg = 
    fDataHists[Detector::kNear].begin();
  std::map<Background::Background_t,NueBackground*>::iterator dataEnd = 
    fDataHists[Detector::kNear].end();
  while(dataBeg!=dataEnd){
    std::map<Systematic::Systematic_t,std::string>::iterator fileBeg = fFileMap.begin();
    std::map<Systematic::Systematic_t,std::string>::iterator fileEnd = fFileMap.end();
    while(fileBeg!=fileEnd){
      std::map<NueSystematic*,TH1D*>::iterator histBeg = 
	fPredictionMap[dataBeg->first][fileBeg->first].begin();
      std::map<NueSystematic*,TH1D*>::iterator histEnd = 
	fPredictionMap[dataBeg->first][fileBeg->first].end();
      while(histBeg!=histEnd){
	delete histBeg->first;
	delete histBeg->second;
	histBeg++;
      }
      fileBeg++;
    }
    dataBeg++;
  }

  //delete data hists:
  dataBeg = fDataHists[Detector::kNear].begin();
  while(dataBeg!=dataEnd) {
    delete dataBeg->second;
    dataBeg++;
  }
  dataBeg = fDataHists[Detector::kFar].begin();
  dataEnd = fDataHists[Detector::kFar].end(); 
  while(dataBeg!=dataEnd) {
    delete dataBeg->second;
    dataBeg++;
  }

  //delete histos/graphs held for plotting:
  this->DeletePredictionsForDrawing();
  this->DeleteRatiosForDrawing();
  this->DeleteIntegralsForDrawing();
  this->DeleteSummaryForDrawing();
  delete fPredictionsForDrawing;
  delete fRatiosForDrawing;
  delete fIntegralsForDrawing;
  delete fSummaryForDrawing;
  delete fCanvasForPredictions;
  delete fCanvasForRatios;
  delete fCanvasForIntegrals;
  delete fCanvasForSummary;
}

void Comparator::AddSysFile(Systematic::Systematic_t sys,string fname)
{  
  fFileMap[sys] = fname;
  if(fDataHists.begin()==fDataHists.end()) {
    TFile *f = new TFile(fname.c_str(),"READ");
    if(!this->ExtractDataHists(f)) 
      std::cout << "Error: Could not extract data histograms from file: " 
		<< fname << endl;
    delete f;
  }
}

Bool_t Comparator::ExtractDataHists(TFile *file)
{
  std::cout << "Extracting Data Hists" << std::endl;
  if(!file) return false;
 
  TTree *tree = (TTree*) file->Get("energytree");
  if(!tree) return false;
  
  tree->SetBranchAddress("nearPOT",&fNDDataPOT);
  tree->SetBranchAddress("farPOT",&fFDDataPOT);
  tree->SetBranchAddress("Selection",fSelection);
  tree->GetEvent(0);
  tree->ResetBranchAddresses();

  std::vector<Background::Background_t>::iterator bgBeg = fBgVec.begin();
  std::vector<Background::Background_t>::iterator bgEnd = fBgVec.end();
  while(bgBeg!=bgEnd){
    string histname = string(Background::AsString(*bgBeg)) + "_Standard/ND_" + fHistType;
    if((*bgBeg)==Background::kNuTauCC || (*bgBeg)==Background::kNueCC) {
      histname = string(Background::AsString(Background::kNuMuCC)) + 
	"_Standard/ND_" + fHistType;
    }
    TH1D *hist = (TH1D*) file->Get(histname.c_str());
    if(!hist) return false;
    string nuename = "ND_" + fHistType + "_" + string(fSelection) + "_" + 
      string(Background::AsString(*bgBeg)) + "_DATA";
    NueBackground *nbg = new NueBackground(nuename,hist,Detector::kNear,
					   (*bgBeg),Selection::StringToEnum(fSelection),
					   fNDDataPOT);
    fDataHists[Detector::kNear][*bgBeg] = nbg;

    histname = string(Background::AsString(*bgBeg)) + "_Standard/FD_" + fHistType;
    hist = (TH1D*) file->Get(histname.c_str());
    if(!hist) return false;
    nuename = "FD_" + fHistType + "_" + string(fSelection) + "_" + 
      string(Background::AsString(*bgBeg)) + "_DATA";
    NueBackground *fbg = new NueBackground(nuename,hist,Detector::kFar,
					   (*bgBeg),Selection::StringToEnum(fSelection),
					   fFDDataPOT);
    fDataHists[Detector::kFar][*bgBeg] = fbg;
    bgBeg++;
  }
  return true;
}

void Comparator::ComputeAll()
{
  std::map<Background::Background_t,NueBackground*>::iterator dataBeg = 
    fDataHists[Detector::kNear].begin();
  std::map<Background::Background_t,NueBackground*>::iterator dataEnd = 
    fDataHists[Detector::kNear].end();
  while(dataBeg!=dataEnd){
    std::map<Systematic::Systematic_t,std::string>::iterator fileBeg = fFileMap.begin();
    std::map<Systematic::Systematic_t,std::string>::iterator fileEnd = fFileMap.end();
    while(fileBeg!=fileEnd){
      fPredictionMap[dataBeg->first][fileBeg->first] = GetPredictions(fileBeg->first,
								      dataBeg->first);
      fileBeg++;
    }
    dataBeg++;
  }
}

std::map<NueSystematic*,TH1D*> Comparator::GetPredictions(Systematic::Systematic_t sysType,
							  Background::Background_t bgType)
{
  cout << "Getting Predictions for " << Systematic::AsString(sysType) << " and "
       << Background::AsString(bgType) << endl;
  std::map<NueSystematic*,TH1D*> histVec;
  NueFNExtrapolation *extra = new NueFNExtrapolation(fFileMap[sysType]);
  Int_t max_sys = extra->GetMaxTreeEntry();
  for(int i=0;i<max_sys;i++){
    extra->UseTreeEntry(i);
    NueSystematic *nueSys = extra->GetCurrentSystematic();
    TH1D *histo = (TH1D*) extra->GetSpectrum(Detector::kFar,fFDDataPOT,
					     fDataHists[Detector::kNear][bgType] );
    string nomus = string(histo->GetName()) + "_" + string(nueSys->GetName());
    histo->SetName(nomus.c_str());
    histo->SetDirectory(0);
    histVec[nueSys] = histo;
  }
  delete extra;
  return histVec;
}

void Comparator::DrawAll(Int_t whichSys) 
{
  this->DrawPrediction(whichSys);
  this->DrawRatio(whichSys);
  this->DrawIntegral(whichSys);
}

void Comparator::DrawAll(Int_t whichSys,Int_t whichBG) 
{
  this->DrawPrediction(whichSys,whichBG);
  this->DrawRatio(whichSys,whichBG);
  this->DrawIntegral(whichSys,whichBG);
}

void Comparator::DrawPrediction(Int_t whichSys)
{
  //delete previously made/drawn histograms
  this->DeletePredictionsForDrawing();

  //get/make TLegend
  fCanvasForPredictions->cd();
  TLegend *leg = 0;
  if((leg = (TLegend*) fCanvasForPredictions->FindObject("CFPLeg"))) leg->Clear();
  else {
    leg = new TLegend(0.2,0.2,0.4,0.4);
    leg->SetName("CFPLeg");
    leg->SetBorderSize(1);
  }

  Systematic::Systematic_t sys = Systematic::ESystematic(whichSys);
  if(fFileMap.find(sys)==fFileMap.end()) return;
  std::map<string,TH1D*> summedMap;
  TH1D *dataHist = NULL;

  std::map<Background::Background_t,NueBackground*>::iterator dataBeg = 
    fDataHists[Detector::kNear].begin();
  std::map<Background::Background_t,NueBackground*>::iterator dataEnd = 
    fDataHists[Detector::kNear].end();
  while(dataBeg!=dataEnd){
    Background::Background_t bg = dataBeg->first;
    std::map<NueSystematic*,TH1D*>::iterator predBeg = fPredictionMap[bg][sys].begin();
    std::map<NueSystematic*,TH1D*>::iterator predEnd = fPredictionMap[bg][sys].end();
    Int_t colCounter = 0;
    char label[256];
    while(predBeg!=predEnd){
      string sysNom = predBeg->first->GetName();
      Double_t sysVal = predBeg->first->GetSysValue(sys);
      if(sys==Systematic::kOscProb) {
	TString ts(sysNom);
	if(!ts.Contains(fOscSysString)) { predBeg++; continue;}
	Double_t tmpDub, theta23, delta23;
	Int_t tmpInt;
	predBeg->first->GetOscParams(tmpDub,theta23,tmpDub,tmpDub,delta23,tmpDub,tmpInt);	
	if(ts.Contains("Theta23")) sysVal = theta23;
	else if(ts.Contains("Delta23")) sysVal = delta23;
      }
      if(summedMap.find(sysNom)!=summedMap.end()) {
	summedMap[sysNom]->Add(predBeg->second);
      }
      else {
	string name = "FD_Total_Prediction_" + string(fSelection) + "_" + 
	  string(Systematic::AsString(sys));
	summedMap[sysNom] = (TH1D*) predBeg->second->Clone(name.c_str());
	summedMap[sysNom]->SetLineColor(fColourArray[colCounter]);
	summedMap[sysNom]->SetMarkerColor(fColourArray[colCounter]);
	if(sys!=Systematic::kOscProb) sprintf(label,"%s = %.2f",
					      Systematic::AsString(sys),sysVal);
	else sprintf(label,"%s = %.4f",fOscSysString.c_str(),sysVal);
	leg->AddEntry(summedMap[sysNom],label,"l");
	colCounter++;
      }
      predBeg++;
    }
    if(dataHist==0) {
      dataHist = (TH1D*) fDataHists[Detector::kFar][dataBeg->first]->
	GetHist()->Clone("DataHist");
      leg->AddEntry(dataHist,"Pseudo-Data","lp");
    }
    else dataHist->Add(fDataHists[Detector::kFar][dataBeg->first]->GetHist());
    dataBeg++;
  }

  std::map<string,TH1D*>::iterator sumBeg = summedMap.begin();
  std::map<string,TH1D*>::iterator sumEnd = summedMap.end();
  while(sumBeg!=sumEnd){
    sumBeg->second->SetAxisRange(0,10);
    fPredictionsForDrawing->Add(sumBeg->second,"hist");
    sumBeg++;
  }
  fPredictionsForDrawing->Add(dataHist);

  //draw:
  fPredictionsForDrawing->Draw("nostack");
  fPredictionsForDrawing->GetXaxis()->SetRangeUser(0,10);
  string title = "FD Total Prediction with " + 
    string(fSelection) + " and Tweaks in " + string(Systematic::AsString(sys));
  fPredictionsForDrawing->SetTitle(title.c_str());
  leg->Draw();
  gROOT->ProcessLine("SortOutStats(gPad,0.25,0.6,1,0.98,0.81);STS(0.04);");
  if(fDoPrint) {
    string printName = "FD_Total_Prediction_" + 
      string(fSelection) + "_" + string(Systematic::AsString(sys)) + ".gif";
    fCanvasForPredictions->Print(printName.c_str());
  }
}

void Comparator::DrawPrediction(Int_t whichSys,Int_t whichBG)
{
  //delete previously made/drawn histograms  
  this->DeletePredictionsForDrawing();
 
  //get/make TLegend
  fCanvasForPredictions->cd();  
  TLegend *leg = 0;
  if((leg = (TLegend*) fCanvasForPredictions->FindObject("CFPLeg"))) leg->Clear();
  else {
    leg = new TLegend(0.2,0.2,0.4,0.4);
    leg->SetName("CFPLeg");
    leg->SetBorderSize(1);
  }

  Int_t colCounter = 0;
  Systematic::Systematic_t sys = Systematic::ESystematic(whichSys);
  Background::Background_t bg = Background::EBackground(whichBG);
  if(!this->InBGVector(bg)) return;
  if(fFileMap.find(sys)==fFileMap.end()) return;
  std::map<NueSystematic*,TH1D*>::iterator predBeg = fPredictionMap[bg][sys].begin();
  std::map<NueSystematic*,TH1D*>::iterator predEnd = fPredictionMap[bg][sys].end();
  char label[256];
  while(predBeg!=predEnd) {
    NueSystematic *tmpSys = predBeg->first;
    if(sys==Systematic::kOscProb) {
      TString ts(tmpSys->GetName());
      if(!ts.Contains(fOscSysString)) { predBeg++; continue;}
      Double_t tmpDub, theta23, delta23;
      Int_t tmpInt;
      tmpSys->GetOscParams(tmpDub,theta23,tmpDub,tmpDub,delta23,tmpDub,tmpInt);
      if(ts.Contains("Theta23")) sprintf(label,"%s = %.4f",tmpSys->GetName(),theta23);
      else if(ts.Contains("Delta23")) sprintf(label,"%s = %.4f",tmpSys->GetName(),delta23);
    }
    else sprintf(label,"%s = %.2f",Systematic::AsString(sys),tmpSys->GetSysValue(sys));
    TH1D *tmpPred = (TH1D*) predBeg->second->
      Clone(string(string(predBeg->second->GetName())+string("_copy")).c_str());
    tmpPred->SetLineColor(fColourArray[colCounter]);
    tmpPred->SetMarkerColor(fColourArray[colCounter]);
    tmpPred->SetAxisRange(0,10);
    fPredictionsForDrawing->Add(tmpPred,"hist");
    leg->AddEntry(tmpPred,label,"lp");
    colCounter++;
    predBeg++;
  }

  //Draw:
  fPredictionsForDrawing->Draw("nostack");
  fPredictionsForDrawing->GetXaxis()->SetRangeUser(0,10);
  TH1D *dataHist = 
        (TH1D*) fDataHists[Detector::kFar][bg]->GetHist()->Clone("DataHist");
  dataHist->SetAxisRange(0,10);
  dataHist->Draw("same");
  leg->AddEntry(dataHist,"Pseudo-Data","lp");
  string title = "FD Prediction for " + string(Background::AsString(bg)) + " with " + 
    string(fSelection) + " and Tweaks in " + string(Systematic::AsString(sys));
  fPredictionsForDrawing->SetTitle(title.c_str());
  leg->Draw();
  gROOT->ProcessLine("SortOutStats(gPad,0.25,0.6,1,0.98,0.81);STS(0.04);");
  if(fDoPrint) {
    string printName = "FD_Prediction_" + string(Background::AsString(bg)) + "_" +
      string(fSelection) + "_" + string(Systematic::AsString(sys)) + ".gif";
    fCanvasForPredictions->Print(printName.c_str());
  }
}

void Comparator::DrawRatio(Int_t whichSys)
{
  //delete previously made/drawn histograms
  this->DeleteRatiosForDrawing();
  
  //get/make TLegend
  fCanvasForRatios->cd();  
  TLegend *leg = 0;
  if((leg = (TLegend*) fCanvasForRatios->FindObject("CFRLeg"))) leg->Clear();
  else {
    leg = new TLegend(0.2,0.2,0.4,0.4);
    leg->SetName("CFRLeg");
    leg->SetBorderSize(1);
  }

  Systematic::Systematic_t sys = Systematic::ESystematic(whichSys);
  if(fFileMap.find(sys)==fFileMap.end()) return;
  std::map<string,TH1D*> summedMap;
  TH1D *dataHist = NULL;

  std::map<Background::Background_t,NueBackground*>::iterator dataBeg = 
    fDataHists[Detector::kNear].begin();
  std::map<Background::Background_t,NueBackground*>::iterator dataEnd = 
    fDataHists[Detector::kNear].end();
  while(dataBeg!=dataEnd){
    Background::Background_t bg = dataBeg->first;
    std::map<NueSystematic*,TH1D*>::iterator predBeg = fPredictionMap[bg][sys].begin();
    std::map<NueSystematic*,TH1D*>::iterator predEnd = fPredictionMap[bg][sys].end();
    Int_t colCounter = 0;
    char label[256];
    while(predBeg!=predEnd){
      Double_t sysVal = predBeg->first->GetSysValue(sys);
      string sysNom = predBeg->first->GetName();
      if(sys==Systematic::kOscProb) {
	TString ts(sysNom);
	if(!ts.Contains(fOscSysString)) { predBeg++; continue;}
	Double_t tmpDub, theta23, delta23;
	Int_t tmpInt;
	predBeg->first->GetOscParams(tmpDub,theta23,tmpDub,tmpDub,delta23,tmpDub,tmpInt);
	if(ts.Contains("Theta23")) sysVal = theta23;
	else if(ts.Contains("Delta23")) sysVal = delta23;
      }
      if(summedMap.find(sysNom)!=summedMap.end()) 
	summedMap[sysNom]->Add(predBeg->second);
      else {
	string name = "FN_Total_Ratio_" + string(fSelection) + "_" + 
	  string(Systematic::AsString(sys));
	summedMap[sysNom] = (TH1D*) predBeg->second->Clone(name.c_str());
	summedMap[sysNom]->SetLineColor(fColourArray[colCounter]);
	summedMap[sysNom]->SetMarkerColor(fColourArray[colCounter]);
	if(sys!=Systematic::kOscProb) sprintf(label,"%s = %.2f",
					      Systematic::AsString(sys),sysVal);
	else sprintf(label,"%s = %.4f",fOscSysString.c_str(),sysVal);
	leg->AddEntry(summedMap[sysNom],label,"lp");
	colCounter++;
      }
      predBeg++;
    }
    if(dataHist==0) {
      dataHist = (TH1D*) fDataHists[Detector::kNear][dataBeg->first]
	->GetHist()->Clone("tmpHist");
      //leg->AddEntry(dataHist,"Pseudo-Data","lp");
    }
    else dataHist->Add(fDataHists[Detector::kNear][dataBeg->first]->GetHist());
    dataBeg++;
  }

  std::map<string,TH1D*>::iterator sumBeg = summedMap.begin();
  std::map<string,TH1D*>::iterator sumEnd = summedMap.end();
  while(sumBeg!=sumEnd){
    sumBeg->second->Divide(dataHist);
    sumBeg->second->Scale(fNDDataPOT/fFDDataPOT);
    sumBeg->second->SetAxisRange(0,10);
    fRatiosForDrawing->Add(sumBeg->second,"hist");
    sumBeg++;
  }
  delete dataHist;

  //draw:
  fRatiosForDrawing->Draw("nostack");
  fRatiosForDrawing->GetXaxis()->SetRangeUser(0,10);
  string title = "F/N Total Ratio with " + string(fSelection) + 
    " and Tweaks in " + string(Systematic::AsString(sys));
  fRatiosForDrawing->SetTitle(title.c_str());
  leg->Draw();
  gROOT->ProcessLine("SortOutStats(gPad,0.25,0.6,1,0.98,0.81);STS(0.04);");
  if(fDoPrint) {
    string printName = "FN_Total_Ratio_" + string(fSelection) + "_" + 
      string(Systematic::AsString(sys)) + ".gif";
    fCanvasForRatios->Print(printName.c_str());
  }
}

void Comparator::DrawRatio(Int_t whichSys,Int_t whichBG)
{
  //delete previously made/drawn histograms  
  this->DeleteRatiosForDrawing();
  
  //get/make TLegend
  fCanvasForRatios->cd();  
  TLegend *leg = 0;
  if((leg = (TLegend*) fCanvasForRatios->FindObject("CFRLeg"))) leg->Clear();
  else {
    leg = new TLegend(0.2,0.2,0.4,0.4);
    leg->SetName("CFRLeg");
    leg->SetBorderSize(1);
  }
 
  Systematic::Systematic_t sys = Systematic::ESystematic(whichSys);
  Background::Background_t bg = Background::EBackground(whichBG);
  if(!this->InBGVector(bg)) return;
  if(fFileMap.find(sys)==fFileMap.end()) return;
  std::map<NueSystematic*,TH1D*>::iterator predBeg = fPredictionMap[bg][sys].begin();
  std::map<NueSystematic*,TH1D*>::iterator predEnd = fPredictionMap[bg][sys].end();
  Int_t colCounter = 0;
  char label[256];
  while(predBeg!=predEnd) {
    NueSystematic *tmpSys = predBeg->first;
    if(sys==Systematic::kOscProb) {
      TString ts(tmpSys->GetName());
      if(!ts.Contains(fOscSysString)) { predBeg++; continue;}
      Double_t tmpDub, theta23, delta23;
      Int_t tmpInt;
      tmpSys->GetOscParams(tmpDub,theta23,tmpDub,tmpDub,delta23,tmpDub,tmpInt);
      if(ts.Contains("Theta23")) sprintf(label,"%s = %.4f",fOscSysString.c_str(),theta23);
      else if(ts.Contains("Delta23")) sprintf(label,"%s = %.4f",fOscSysString.c_str(),delta23);
    }
    else sprintf(label,"%s = %.2f",Systematic::AsString(sys),tmpSys->GetSysValue(sys));
    TH1D *tmpPred = (TH1D*) predBeg->second->Clone(string(predBeg->second->GetName() + 
							  string("_Copy")).c_str());
    tmpPred->Divide(fDataHists[Detector::kNear][bg]->GetHist());
    tmpPred->Scale(fNDDataPOT/fFDDataPOT);
    tmpPred->SetLineColor(fColourArray[colCounter]);
    tmpPred->SetMarkerColor(fColourArray[colCounter]);
    tmpPred->SetAxisRange(0,10);
    fRatiosForDrawing->Add(tmpPred,"hist");
    leg->AddEntry(tmpPred,label,"lp");
    colCounter++;
    predBeg++;
  }

  //draw:
  fRatiosForDrawing->Draw("nostack");
  fRatiosForDrawing->GetXaxis()->SetRangeUser(0,10);
  string title = "F/N Ratio for " + string(Background::AsString(bg)) + " with " + 
    string(fSelection) + " and Tweaks in " + string(Systematic::AsString(sys));
  fRatiosForDrawing->SetTitle(title.c_str());
  leg->Draw();
  gROOT->ProcessLine("SortOutStats(gPad,0.25,0.6,1,0.98,0.81);STS(0.04);");
  if(fDoPrint) {
    string printName = "FN_Ratio_" + string(Background::AsString(bg)) + "_" + 
      string(fSelection) + "_" + string(Systematic::AsString(sys)) + ".gif";
    fCanvasForRatios->Print(printName.c_str());
  }
}

void Comparator::DrawIntegral(Int_t whichSys)
{
  //delete previous TGraphs
  this->DeleteIntegralsForDrawing();
  //get/make TLegend
  fCanvasForIntegrals->cd();
  TLegend *leg = (TLegend*) fCanvasForIntegrals->FindObject("CFILeg");
  if(leg) leg->Clear();
  else {
    leg = new TLegend(0.2,0.2,0.4,0.4);
    leg->SetName("CFILeg");
    leg->SetBorderSize(0);
  }
  fCanvasForIntegrals->Clear();

  Int_t n = 0;
  Double_t sysValues[100] = {0}; //unfeasibly large...
  Double_t mcIntegValues[100] = {0};
  Double_t extrapIntegValues[100] = {0};
  Double_t dataTot = 0;
  Bool_t firstPass = true;
  Systematic::Systematic_t sys = Systematic::ESystematic(whichSys);
  if(fFileMap.find(sys)==fFileMap.end()) return;

  std::map<Background::Background_t,NueBackground*>::iterator dataBeg = 
    fDataHists[Detector::kNear].begin();
  std::map<Background::Background_t,NueBackground*>::iterator dataEnd = 
    fDataHists[Detector::kNear].end();
  while(dataBeg!=dataEnd){
    Background::Background_t bg = dataBeg->first;
    std::map<NueSystematic*,TH1D*>::iterator predBeg = fPredictionMap[bg][sys].begin();
    std::map<NueSystematic*,TH1D*>::iterator predEnd = fPredictionMap[bg][sys].end();    
    while(predBeg!=predEnd){
      Double_t theVal = 0;
      if(sys==Systematic::kOscProb) {
	TString ts(predBeg->first->GetName());
	if(!ts.Contains(fOscSysString)) { predBeg++; continue;}
	Double_t tmpDub, theta23, delta23;
	Int_t tmpInt;
	predBeg->first->GetOscParams(tmpDub,theta23,tmpDub,tmpDub,delta23,tmpDub,tmpInt);	
	if(ts.Contains("Theta23")) theVal = theta23;
	else if(ts.Contains("Delta23")) theVal = delta23;
      }
      else theVal = predBeg->first->GetSysValue(sys);
      if(firstPass) {
	sysValues[n] = theVal;
	extrapIntegValues[n] += predBeg->second->Integral();
	mcIntegValues[n] += this->GetFDSpectrumIntegral(bg,sys,sysValues[n]);
	n++;
      }
      else {
	for(int i=0;i<n;i++) {
	  if(TMath::Abs(sysValues[i]-theVal)<1e-6) {
	    extrapIntegValues[i] += predBeg->second->Integral();
	    mcIntegValues[i] += this->GetFDSpectrumIntegral(bg,sys,sysValues[i]);
	    break;
	  }
	}
      }
      predBeg++;
    }
    dataTot += fDataHists[Detector::kFar][bg]->GetHist()->Integral();
    firstPass = false;
    dataBeg++;
  }
  if(dataTot>0){
    for(int i=0;i<n;i++){
      mcIntegValues[i]/=dataTot;
      extrapIntegValues[i]/=dataTot;
    }
  }

  TGraph *extrapGr = new TGraph(n,sysValues,extrapIntegValues);
  TGraph *mcGr = new TGraph(n,sysValues,mcIntegValues);
  extrapGr->SetMarkerColor(4);
  mcGr->SetMarkerColor(2);
  fIntegralsForDrawing->Add(extrapGr,"P");
  fIntegralsForDrawing->Add(mcGr,"P");
  leg->AddEntry(extrapGr,"Data_{ND}#times(F/N)_{sys}/Data_{FD}","p");
  leg->AddEntry(mcGr,"MC^{sys}_{FD}/Data_{FD}","p");

  fIntegralsForDrawing->Draw("A");
  string title = "FD Total Spectrum Integral Ratios with " + string(fSelection) + 
    " and Tweaks in " + string(Systematic::AsString(sys));
  fIntegralsForDrawing->SetTitle(title.c_str());  
  leg->Draw();
  gROOT->ProcessLine("SortOutStats(gPad,0.3,0.075,1,0.7,0.09);STS(0.04);");
  if(fDoPrint) {
    string printName = "FD_Integral_" + string(fSelection) + "_" + 
      string(Systematic::AsString(sys)) + ".gif";
    fCanvasForIntegrals->Print(printName.c_str());
  }
}

void Comparator::DrawIntegral(Int_t whichSys,Int_t whichBG)
{
  //delete previous TGraphs
  this->DeleteIntegralsForDrawing();
  //get/make TLegend
  fCanvasForIntegrals->cd();  
  TLegend *leg = (TLegend*) fCanvasForIntegrals->FindObject("CFILeg");
  if(leg) leg->Clear();
  else {
    leg = new TLegend(0.2,0.2,0.4,0.4);
    leg->SetName("CFILeg");
    leg->SetBorderSize(0);
  }
  fCanvasForIntegrals->Clear();

  Int_t n = 0;
  Double_t sysValues[100] = {0}; //unfeasibly large...
  Double_t mcIntegValues[100] = {0};
  Double_t extrapIntegValues[100] = {0};
  Systematic::Systematic_t sys = Systematic::ESystematic(whichSys);
  Background::Background_t bg = Background::EBackground(whichBG);
  if(!this->InBGVector(bg)) return;
  if(fFileMap.find(sys)==fFileMap.end()) return;
  std::map<NueSystematic*,TH1D*>::iterator predBeg = fPredictionMap[bg][sys].begin();
  std::map<NueSystematic*,TH1D*>::iterator predEnd = fPredictionMap[bg][sys].end();
  while(predBeg!=predEnd) {
    if(sys==Systematic::kOscProb) {
      TString ts(predBeg->first->GetName());
      if(!ts.Contains(fOscSysString)) { predBeg++; continue;}
      Double_t tmpDub, theta23, delta23;
      Int_t tmpInt;
      predBeg->first->GetOscParams(tmpDub,theta23,tmpDub,tmpDub,delta23,tmpDub,tmpInt);	
      if(ts.Contains("Theta23")) sysValues[n] = theta23;
      else if(ts.Contains("Delta23")) sysValues[n] = delta23;
    }
    else sysValues[n] = predBeg->first->GetSysValue(sys);
    //if(sysValues[n]<-9998) sysValues[n] = 0;
    Double_t dataTot = fDataHists[Detector::kFar][bg]->GetHist()->Integral();
    extrapIntegValues[n] = predBeg->second->Integral();
    mcIntegValues[n] = this->GetFDSpectrumIntegral(bg,sys,sysValues[n]);
    if(dataTot>0){
      extrapIntegValues[n]/=dataTot;
      mcIntegValues[n]/=dataTot;
    }
    predBeg++;
    n++;
  }
  TGraph *extrapGr = new TGraph(n,sysValues,extrapIntegValues);
  TGraph *mcGr = new TGraph(n,sysValues,mcIntegValues);
  extrapGr->SetMarkerColor(4);
  mcGr->SetMarkerColor(2);
  fIntegralsForDrawing->Add(extrapGr,"P");
  fIntegralsForDrawing->Add(mcGr,"P");
  leg->AddEntry(extrapGr,"Data_{ND}#times(F/N)_{sys}/Data_{FD}","p");
  leg->AddEntry(mcGr,"MC^{sys}_{FD}/Data_{FD}","p");

  fIntegralsForDrawing->Draw("A");
  string title = "FD Spectrum Integral Ratios for " + string(Background::AsString(bg)) + 
    " with " + string(fSelection) + " and Tweaks in " + string(Systematic::AsString(sys));
  fIntegralsForDrawing->SetTitle(title.c_str());  
  leg->Draw();
  gROOT->ProcessLine("SortOutStats(gPad,0.3,0.075,1,0.7,0.09);STS(0.04);");
  if(fDoPrint) {
    string printName = "FD_Integral_" + string(Background::AsString(bg)) + "_" + 
      string(fSelection) + "_" + string(Systematic::AsString(sys)) + ".gif";
    fCanvasForIntegrals->Print(printName.c_str());
  }
}

Double_t Comparator::GetFDSpectrumIntegral(Background::Background_t bg,
					   Systematic::Systematic_t sys, Double_t inputVal)
{
  TFile *f = new TFile(fFileMap[sys].c_str(),"READ");
  TTree *tree = (TTree*) f->Get("energytree");  
  Double_t sysVal = 0;
  tree->SetBranchAddress(Systematic::AsString(sys),&sysVal);
  char sysname[256];
  tree->SetBranchAddress("SysName",sysname);
  Double_t farPOT = 0;
  tree->SetBranchAddress("farPOT",&farPOT);
  Double_t theta23,delta23;
  tree->SetBranchAddress("Theta23",&theta23); 
  tree->SetBranchAddress("DeltaMSq23",&delta23);
  for(int i=0;i<tree->GetEntries();i++){
    tree->GetEntry(i);
    if(sys==Systematic::kOscProb) {
      TString ts(sysname);
      if(!ts.Contains(fOscSysString)) continue;
      if(ts.Contains("Theta23")) {
	if(TMath::Abs(theta23-inputVal)<1e-6) break;
      }
      else if(ts.Contains("Delta23")) {
	if(TMath::Abs(delta23-inputVal)<1e-6) break;
      }
    }
   else if(TMath::Abs(sysVal-inputVal)<1e-6) break;
  }
  string hname = string(Background::AsString(bg)) + "_" + string(sysname) + "/FD_" + fHistType;
  TH1D *hist = (TH1D*) f->Get(hname.c_str());
  Double_t integ = hist->Integral()*(fFDDataPOT/farPOT);
  f->Close();
  delete f;
  return integ;
}

void Comparator::DrawSummary()
{
  this->DeleteSummaryForDrawing();
  fCanvasForSummary->cd();  
  TLegend *leg = (TLegend*) fCanvasForSummary->FindObject("CFSLeg");
  if(leg) leg->Clear();
  else {
    leg = new TLegend(0.2,0.2,0.4,0.4);
    leg->SetName("CFSLeg");
    leg->SetBorderSize(1);
  }
  fCanvasForSummary->Clear();
    
  TGraph *tot_min_sum = this->GetSummary();
  leg->AddEntry(tot_min_sum,"Total","p");
  TIter it(fSummaryForDrawing->GetListOfGraphs()->MakeIterator());
  TGraph *tot_max_sum = 0;
  while((tot_max_sum = (TGraph*) it())) {
    if(tot_min_sum != tot_max_sum) break;
  }
  
  Double_t min_sys = 0;
  Double_t max_sys = 0;
  for(int i=0;i<tot_min_sum->GetN();i++) {
    min_sys += TMath::Power(1-tot_min_sum->GetY()[i],2);
    max_sys += TMath::Power(tot_max_sum->GetY()[i]-1,2);
  }
  min_sys = 100*TMath::Sqrt(min_sys);
  max_sys = 100*TMath::Sqrt(max_sys);
  
  std::map<Background::Background_t,NueBackground*>::iterator dataBeg = 
    fDataHists[Detector::kNear].begin();
  std::map<Background::Background_t,NueBackground*>::iterator dataEnd = 
    fDataHists[Detector::kNear].end();
  while(dataBeg!=dataEnd){
    leg->AddEntry(this->GetSummary(Int_t(dataBeg->first)),
		  Background::AsString(dataBeg->first),"p");
    dataBeg++;
  }

  fCanvasForSummary->Clear();  
  fSummaryForDrawing->Draw("A");
  for(int i=0;i<tot_min_sum->GetN();i++){
    Systematic::Systematic_t sys = Systematic::ESystematic((int)tot_min_sum->GetX()[i]);
    fSummaryForDrawing->GetXaxis()->SetBinLabel(fSummaryForDrawing->GetXaxis()->
						FindBin(tot_min_sum->GetX()[i]),
						Systematic::AsString(sys));
  }
  leg->Draw();

  char tex_string[100]; sprintf(tex_string,"Total Systematic = +%.1f%% -%.1f%%",max_sys,min_sys);
  TLatex *tex = new TLatex(0.015,0.9,tex_string);
  tex->SetNDC();
  tex->SetTextSize(0.045);
  tex->Draw();

  gROOT->ProcessLine("SortOutStats(gPad,0.2,0.2,1,0.98,0.98);STS(0.04);");
  if(fDoPrint) {
    string printName = "SystematicSummary_" + string(fSelection) + ".gif";
    fCanvasForSummary->Print(printName.c_str());
  }
}

TGraph *Comparator::GetSummary(Int_t whichBG)
{

  Int_t n = 0;
  Double_t one[100] = {0};
  Double_t zero[100] = {0};
  Double_t sysNum[100] = {0};
  Double_t sysMin[100] = {0};
  Double_t sysMax[100] = {0};

  //hold current value of OscSysString. Will do Theta23 first then Delta23.
  string oscSt = fOscSysString;
  fOscSysString = "Theta23";

  std::map<Systematic::Systematic_t,std::string>::iterator fileBeg = 
    fFileMap.begin();
  std::map<Systematic::Systematic_t,std::string>::iterator fileEnd = 
    fFileMap.end();
  while(fileBeg!=fileEnd){
    if(whichBG<0) this->DrawIntegral(Int_t(fileBeg->first));
    else this->DrawIntegral(Int_t(fileBeg->first),whichBG);
    TGraph *gr = 0;
    if(fIntegralsForDrawing->GetListOfGraphs()){
      TIter it(fIntegralsForDrawing->GetListOfGraphs()->MakeIterator());
      while((gr = (TGraph*) it())) {
	if(gr->GetMarkerColor()==4) break;
      }
    }
    if(gr) {      
      sysNum[n] = Int_t(fileBeg->first);
      Double_t x=0,y=0;
      for(int i=0;i<gr->GetN();i++){
	gr->GetPoint(i,x,y);
	if(i==0 || y>sysMax[n]) sysMax[n] = y;
	if(i==0 || y<sysMin[n]) sysMin[n] = y;
      }
    }
    one[n] = 1;
    zero[n] = 0;
    n++;

    //if we are doing OscProb, and fOscSysString is Theta23, 
    //then do Delta23, otherwise move onto the next systematic
    if(fileBeg->first==Systematic::kOscProb && 
       fOscSysString!="Delta23") fOscSysString = "Delta23";
    else fileBeg++;
  }

  //revert back to original value of fOscSysString
  fOscSysString = oscSt;

  TGraph *minGr = new TGraph(n,sysNum,sysMin);
  TGraph *maxGr = new TGraph(n,sysNum,sysMax);
  if(whichBG<0) cout << "Summary info for Total background" << endl;
  else {
    Background::Background_t bg = Background::EBackground(whichBG);
    cout << "Summary info for " << Background::AsString(bg) << endl;;
  }
  for(int i=0;i<n;i++){
    Systematic::Systematic_t sys = Systematic::ESystematic((int)sysNum[i]);
    cout << Systematic::AsString(sys) << " " 
	 << sysMin[i] << " " << sysMax[i] << endl;
  }
  if(whichBG>0){
    minGr->SetMarkerColor(fColourArray[whichBG+1]);
    maxGr->SetMarkerColor(fColourArray[whichBG+1]);
    minGr->SetMarkerStyle(24); //open circles
    maxGr->SetMarkerStyle(24);
    minGr->SetMarkerSize(1); //open circles
    maxGr->SetMarkerSize(1);
  }
  fSummaryForDrawing->Add(minGr,"P");
  fSummaryForDrawing->Add(maxGr,"P");
  return minGr;
}


void Comparator::DeletePredictionsForDrawing()
{
  if(fPredictionsForDrawing->GetHists()) {
    TIter it(fPredictionsForDrawing->GetHists()->MakeIterator());
    TObject *ob = 0;
    while((ob = it())) delete ob;
  }
}

void Comparator::DeleteRatiosForDrawing()
{
  if(fRatiosForDrawing->GetHists()) {
    TIter it(fRatiosForDrawing->GetHists()->MakeIterator());
    TObject *ob = 0;
    while((ob = it())) delete ob;
  }
}

void Comparator::DeleteIntegralsForDrawing()
{
  if(fIntegralsForDrawing->GetListOfGraphs()){
    TIter it(fIntegralsForDrawing->GetListOfGraphs()->MakeIterator());
    TObject *ob = 0;
    while((ob = it())) delete ob;
  }
}

void Comparator::DeleteSummaryForDrawing()
{
  if(fSummaryForDrawing->GetListOfGraphs()){
    TIter it(fSummaryForDrawing->GetListOfGraphs()->MakeIterator());
    TObject *ob = 0;
    while((ob = it())) delete ob;
  }
}

Bool_t Comparator::InBGVector(Background::Background_t bg)
{
  std::vector<Background::Background_t>::iterator beg = fBgVec.begin();
  std::vector<Background::Background_t>::iterator end = fBgVec.end();
  while(beg!=end) {
    if((*beg)==bg) return true;
    beg++;
  }
  return false;
}

