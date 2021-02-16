#ifndef JBCOMPARATOR_H
#define JBCOMPARATOR_H
#include "NueAna/NueAnaTools/Selection.h"
#include "NueAna/Extrapolation/Background.h"
#include "NueAna/Extrapolation/NuePrediction.h"
#include "NueAna/Extrapolation/NueFNExtrapolation.h"
#include "NueAna/Extrapolation/NueBackground.h"
#include "Conventions/Detector.h"
#include "THStack.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TGraphAsymmErrors.h"
#include "TFile.h"
#include "TString.h"
#include "TH1D.h"
#include "TLegend.h"
#include "TStyle.h"
#include <string>
#include <vector>
#include <iostream>

class JBComparator
{

 public:
  
  JBComparator(std::string);  //contructor take hist type as string (RecoEnergy/TrueEnergy)
  ~JBComparator();

  void AddBackground(Background::Background_t bg) {fBgVec.push_back(bg);}  
  void AddSysFile(Systematic::Systematic_t,string);  //add file name and associated systematic
  Bool_t ExtractDataHists(TFile *);  //extract "Standard" histograms from 
                                     //first file passed to use as data
  void ComputeAll();

  //draw predictions, ratios, integrals for total prediction or individual backgrounds
  void DrawAll(Int_t);
  void DrawAll(Int_t,Int_t);

  void DrawPrediction(Int_t whichSys);
  void DrawPrediction(Int_t whichSys,Int_t whichBG);

  void DrawRatio(Int_t whichSys);
  void DrawRatio(Int_t whichSys,Int_t whichBG);

  void DrawIntegral(Int_t whichSys);
  void DrawIntegral(Int_t whichSys,Int_t whichBG);
  
  void DrawSummary();
  TGraph *GetSummary(Int_t bg=-1);

  void DoPrint(Bool_t doPrint) {fDoPrint = doPrint;}
  void SetOscSysString(string oss) {fOscSysString = oss;}

  void GetError(Int_t whichSys, double* error, double* minerr, TH1D* Base);
  void DetermineError();

 protected:

  //general information about the files/hists:
  Double_t fNDDataPOT;
  Double_t fFDDataPOT;
  std::string fHistType;
  char fSelection[256];
  Int_t fColourArray[20];
  Bool_t fDoPrint;
  std::string fOscSysString;
  
  //structure to hold systematic<->filename mapping
  std::map<Systematic::Systematic_t,string> fFileMap;

  //structure to hold which backgrounds to consider here:
  std::vector<Background::Background_t> fBgVec;

  //structure to hold data histograms based on detector and background
  //data histograms held as NueBackground objects since these have POT info, etc.
  std::map<Detector::Detector_t,std::map<Background::Background_t,NueBackground*> > fDataHists;

  //map to hold Predictions for all backgrounds and systematics
  std::map<Background::Background_t,std::map<Systematic::Systematic_t,
    std::map<NueSystematic*,TH1D*> > > fPredictionMap;
  
  //function to extract predictions for a specific background and specific systematic
  std::map<NueSystematic*,TH1D*> GetPredictions(Systematic::Systematic_t,
						Background::Background_t);

  Double_t GetFDSpectrumIntegral(Background::Background_t,Systematic::Systematic_t, Double_t);
  Bool_t InBGVector(Background::Background_t);

  TCanvas *fCanvasForPredictions;
  TCanvas *fCanvasForRatios;
  TCanvas *fCanvasForIntegrals;
  TCanvas *fCanvasForSummary;

  void DeletePredictionsForDrawing();
  void DeleteRatiosForDrawing();
  void DeleteIntegralsForDrawing();
  void DeleteSummaryForDrawing();

  THStack *fPredictionsForDrawing;
  THStack *fRatiosForDrawing;  
  TMultiGraph *fIntegralsForDrawing;
  TMultiGraph *fSummaryForDrawing;

  TFile *fRatioFile;

};
#endif //COMPARATOR_H
