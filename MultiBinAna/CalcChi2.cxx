#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include <sys/stat.h>
#include <math.h>
#include <assert.h>
#include <complex>

#include "TH1D.h"
#include "TH2D.h"
#include "TString.h"
#include "THStack.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDEigen.h"
#include "TFile.h"
#include "TDecompSVD.h"
#include "TDecompChol.h"
#include "TMarker.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TLeaf.h"
#include "TMarker.h"
#include "TGraph.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TString.h"

#include "NueAna/MultiBinAna/CalcChi2.h"

using namespace std;

const double k1267 = 1.26693276;
const double kKmUnits = 1000.;

//---------------------------------------------------------------------------------
void PrintParStatus(params my_pars)
{
  std::cout << "" << std::endl;
  std::cout << "/==========Parameter Status==========/" << std::endl;
  std::cout << "" << std::endl;
  std::cout << "Dm232                =      " << my_pars.Dm232 << std::endl;
  std::cout << "Dm221                =      " << my_pars.Dm221 << std::endl;
  std::cout << "Theta23              =      " << my_pars.th23 << std::endl;
  std::cout << "Theta12              =      " << my_pars.th12 << std::endl;
  std::cout << "Theta13              =      " << my_pars.th13 << std::endl;
  std::cout << "DeltaCP              =      " << my_pars.deltaCP << std::endl;
  std::cout << "Dm241                =      " << my_pars.Dm241 << std::endl;
  std::cout << "Theta24              =      " << my_pars.th24 << std::endl;
  std::cout << "Theta34              =      " << my_pars.th34 << std::endl;
  std::cout << "Theta14              =      " << my_pars.th14 << std::endl;
  std::cout << "Delta14              =      " << my_pars.delta14 << std::endl;
  std::cout << "Delta24              =      " << my_pars.delta24 << std::endl;
  std::cout << "" << std::endl;
  std::cout << "/====================================/" << std::endl;
  std::cout << "" << std::endl;
}
//---------------------------------------------------------------------------------
CalcChi2::CalcChi2()
{
  got_beamOptics  = false;
  got_inputHistos = false;
  got_covMx       = false;
  got_unOscHistos = false;
  doDataFit       = false;
  
  fakeDataCC = 0;
  LoadInputHistograms("$SRT_PRIVATE_CONTEXT/NueAna/MultiBinAna/dataRelease.root");
  //Extract covariance matrices
  if (!got_covMx)
    {
      CoVarCC_relative = (TMatrixD*)InFile->Get("TotalCCCovar")->Clone("TotalCCCovar"); assert(CoVarCC_relative);
      got_covMx = true;
    }

  //Construct MINOS/MINOS+ two detector data spectra
  TH1D* h2det_data_CC_minos = (TH1D*)GetTwoDetSpectrum(ND_dataCC_minos,FD_dataCC_minos);
  TH1D* h2det_data_CC_minosPlus = (TH1D*)GetTwoDetSpectrum(ND_dataCC_minosPlus,FD_dataCC_minosPlus);

  //Combine MINOS & MINOS+ data spectra
  dataCC = (TH1D*)h2det_data_CC_minos->Clone();
  dataCC->Add(h2det_data_CC_minosPlus);

  h2det_data_CC_minos->Delete();
  h2det_data_CC_minosPlus->Delete();
  
  if (fakeDataCC == 0)
    {
      double dm232NH = 0.002524;
      double dm232IH = -0.002514;
      double th23NHLO = 0.72626;
      double th23NHUO = 0.84454;
      double th23IHLO = 0.69796;
      double th23IHUO = 0.87284;

      params fake_pars;
      fake_pars.Dm232   = dm232NH;
      //fake_pars.Dm232   = -2.514e-3;
      fake_pars.Dm221   = 0.0000754;
      fake_pars.th23    = th23NHLO;
      fake_pars.th12    = 0.5540758073;
      fake_pars.th13    = 0.149116;
      fake_pars.deltaCP = 0.0;
      fake_pars.Dm241   = 0.0;
      fake_pars.th24    = 0.0;
      fake_pars.th34    = 0.0;
      fake_pars.th14    = 0.0;
      fake_pars.delta24 = 0.0;

      GenerateOscillatedSpectra(fake_pars);
      TH1D* h2det_MC_CC_minos     = (TH1D*)GetTwoDetSpectrum(NDOscCC_MC_minos,FDOscCC_MC_minos);
      TH1D* h2det_MC_CC_minosPlus = (TH1D*)GetTwoDetSpectrum(NDOscCC_MC_minosPlus,FDOscCC_MC_minosPlus);

      fakeDataCC = (TH1D*)h2det_MC_CC_minos->Clone();
      fakeDataCC->Add(h2det_MC_CC_minosPlus);

      h2det_MC_CC_minos->Delete();
      h2det_MC_CC_minosPlus->Delete();
    }

  //Generate predictions
  GenerateUnOscillatedSpectra();
}
//---------------------------------------------------------------------------------
double CalcChi2::Chi2(double Dm232, double Dm221, double th23,
		      double th12, double th13,double deltaCP,
		      double Dm241, double th24, double th34,
		      double th14, double delta14, double delta24)
{
  params my_pars;
  my_pars.Dm232   = Dm232;
  my_pars.Dm221   = Dm221;
  my_pars.th23    = th23;
  my_pars.th12    = th12;
  my_pars.th13    = th13;
  my_pars.deltaCP = deltaCP;
  my_pars.Dm241   = Dm241;
  my_pars.th24    = th24;
  my_pars.th34    = th34;
  my_pars.th14    = th14;
  my_pars.delta14 = delta14;
  my_pars.delta24 = delta24;
  
  GenerateOscillatedSpectra(my_pars);

  //Construct MINOS/MINOS+ two detector prediction spectra
  TH1D* h2det_MC_CC_minos = (TH1D*)GetTwoDetSpectrum(NDOscCC_MC_minos,FDOscCC_MC_minos);
  TH1D* h2det_MC_CC_minosPlus = (TH1D*)GetTwoDetSpectrum(NDOscCC_MC_minosPlus,FDOscCC_MC_minosPlus);

  //Combine MINOS & MINOS+ prediction spectra
  TH1D* predCC = (TH1D*)h2det_MC_CC_minos->Clone();
  predCC->Add(h2det_MC_CC_minosPlus);

  //Scale and invert covariance matrices for fitting
  CoVarCC_inverted = (TMatrixD*)ScaleCovarianceMatrix(predCC,CoVarCC_relative);

  Double_t chi2(0.0);

  if (doDataFit)
    {
      chi2 = ComparePredWithData(predCC,
				 dataCC,
				 CoVarCC_inverted,
				 my_pars.Dm232);
    }
  else
    {
      chi2 = ComparePredWithData(predCC,
				 fakeDataCC,
				 CoVarCC_inverted,
				 my_pars.Dm232);
    }
  predCC->Delete();
  h2det_MC_CC_minos->Delete();
  h2det_MC_CC_minosPlus->Delete();
  NDOscCC_MC_minos->Delete();
  FDOscCC_MC_minos->Delete();
  NDOscCC_MC_minosPlus->Delete();
  FDOscCC_MC_minosPlus->Delete();
  return chi2;
}
//---------------------------------------------------------------------------------
void CalcChi2::Zombie(TFile* f){

  if(f->IsZombie() && (!f->IsOpen())){
    std::cout << "File " << f->GetName() << " failed to open." << std::endl;
    assert(false);
  }
  else{
    std::cout << "File " << f->GetName() << " opened successfully" << std::endl;
  }
}
//---------------------------------------------------------------------------------
void CalcChi2::LoadInputHistograms(TString fileName)
{
  if(got_inputHistos) return;

  InFile  = new TFile(fileName);

  Zombie(InFile);
  
  //Extract RecoToTrue MC simulations for MINOS
  InFile->GetObject("hRecoToTrueNDCCSelectedTrueNC_minos",   NDCC_TrueNC_minos);   assert(NDCC_TrueNC_minos);
  InFile->GetObject("hRecoToTrueNDCCSelectedNuMu_minos",     NDCC_NuMu_minos);     assert(NDCC_NuMu_minos);
  InFile->GetObject("hRecoToTrueNDCCSelectedBeamNue_minos",  NDCC_BeamNue_minos);  assert(NDCC_BeamNue_minos);
  InFile->GetObject("hRecoToTrueNDCCSelectedAppNue_minos",   NDCC_AppNue_minos);   assert(NDCC_AppNue_minos);
  InFile->GetObject("hRecoToTrueNDCCSelectedAppNuTau_minos", NDCC_AppNuTau_minos); assert(NDCC_AppNuTau_minos);

  InFile->GetObject("hRecoToTrueFDCCSelectedTrueNC_minos",   FDCC_TrueNC_minos);   assert(FDCC_TrueNC_minos);
  InFile->GetObject("hRecoToTrueFDCCSelectedNuMu_minos",     FDCC_NuMu_minos);     assert(FDCC_NuMu_minos);
  InFile->GetObject("hRecoToTrueFDCCSelectedBeamNue_minos",  FDCC_BeamNue_minos);  assert(FDCC_BeamNue_minos);
  InFile->GetObject("hRecoToTrueFDCCSelectedAppNue_minos",   FDCC_AppNue_minos);   assert(FDCC_AppNue_minos);
  InFile->GetObject("hRecoToTrueFDCCSelectedAppNuTau_minos", FDCC_AppNuTau_minos); assert(FDCC_AppNuTau_minos);

  //Extract RecoToTrue MC simulations for MINOS+
  InFile->GetObject("hRecoToTrueNDCCSelectedTrueNC_minosPlus",   NDCC_TrueNC_minosPlus);   assert(NDCC_TrueNC_minosPlus);
  InFile->GetObject("hRecoToTrueNDCCSelectedNuMu_minosPlus",     NDCC_NuMu_minosPlus);     assert(NDCC_NuMu_minosPlus);
  InFile->GetObject("hRecoToTrueNDCCSelectedBeamNue_minosPlus",  NDCC_BeamNue_minosPlus);  assert(NDCC_BeamNue_minosPlus);
  InFile->GetObject("hRecoToTrueNDCCSelectedAppNue_minosPlus",   NDCC_AppNue_minosPlus);   assert(NDCC_AppNue_minosPlus);
  InFile->GetObject("hRecoToTrueNDCCSelectedAppNuTau_minosPlus", NDCC_AppNuTau_minosPlus); assert(NDCC_AppNuTau_minosPlus);

  InFile->GetObject("hRecoToTrueFDCCSelectedTrueNC_minosPlus",   FDCC_TrueNC_minosPlus);   assert(FDCC_TrueNC_minosPlus);
  InFile->GetObject("hRecoToTrueFDCCSelectedNuMu_minosPlus",     FDCC_NuMu_minosPlus);     assert(FDCC_NuMu_minosPlus);
  InFile->GetObject("hRecoToTrueFDCCSelectedBeamNue_minosPlus",  FDCC_BeamNue_minosPlus);  assert(FDCC_BeamNue_minosPlus);
  InFile->GetObject("hRecoToTrueFDCCSelectedAppNue_minosPlus",   FDCC_AppNue_minosPlus);   assert(FDCC_AppNue_minosPlus);
  InFile->GetObject("hRecoToTrueFDCCSelectedAppNuTau_minosPlus", FDCC_AppNuTau_minosPlus); assert(FDCC_AppNuTau_minosPlus);
  
  //Extract data histograms
  //MINOS
  InFile->GetObject("dataFDCC_minos", FD_dataCC_minos); assert(FD_dataCC_minos);
  InFile->GetObject("dataNDCC_minos", ND_dataCC_minos); assert(ND_dataCC_minos);
  
  //MINOS+
  InFile->GetObject("dataFDCC_minosPlus", FD_dataCC_minosPlus); assert(FD_dataCC_minosPlus);
  InFile->GetObject("dataNDCC_minosPlus", ND_dataCC_minosPlus); assert(ND_dataCC_minosPlus);  
  
  got_inputHistos = true;
}
//---------------------------------------------------------------------------------
void CalcChi2::GenerateUnOscillatedSpectra()
{
  if (got_unOscHistos) return;

  params my_pars0;
  my_pars0.Dm232 = 0.0;
  my_pars0.Dm221 = 0.0;
  my_pars0.th23  = 0.0;
  my_pars0.th13  = 0.0;
  my_pars0.th12  = 0.0;
  my_pars0.deltaCP = 0.0;
  my_pars0.Dm241   = 0.0;
  my_pars0.th24    = 0.0;
  my_pars0.th34    = 0.0;
  my_pars0.th14    = 0.0;
  my_pars0.delta14 = 0.0;
  my_pars0.delta24 = 0.0;

  //UnOscillated CC MC -- MINOS
  NDUnOscCC_MC_minos = (TH1D*)CreateTotalSpectrum(	my_pars0, 
							NDCC_TrueNC_minos, 
							NDCC_NuMu_minos, 
							NDCC_BeamNue_minos, 
							NDCC_AppNue_minos, 
							NDCC_AppNuTau_minos, 
							1.04*kKmUnits);

  FDUnOscCC_MC_minos = (TH1D*)CreateTotalSpectrum(	my_pars0, 
							FDCC_TrueNC_minos, 
							FDCC_NuMu_minos, 
							FDCC_BeamNue_minos, 
							FDCC_AppNue_minos, 
							FDCC_AppNuTau_minos, 
							735.0*kKmUnits);

  //UnOscillated CC MC -- MINOS+
  NDUnOscCC_MC_minosPlus = (TH1D*)CreateTotalSpectrum(	my_pars0, 
							NDCC_TrueNC_minosPlus, 
							NDCC_NuMu_minosPlus, 
							NDCC_BeamNue_minosPlus, 
							NDCC_AppNue_minosPlus, 
							NDCC_AppNuTau_minosPlus, 
							1.04*kKmUnits);
  FDUnOscCC_MC_minosPlus = (TH1D*)CreateTotalSpectrum(	my_pars0, 
							FDCC_TrueNC_minosPlus, 
							FDCC_NuMu_minosPlus, 
							FDCC_BeamNue_minosPlus, 
							FDCC_AppNue_minosPlus, 
							FDCC_AppNuTau_minosPlus, 
							735.0*kKmUnits);

  got_unOscHistos = true;
}
//---------------------------------------------------------------------------------
void CalcChi2::GenerateOscillatedSpectra(params my_pars)
{
  //Oscillated CC MC -- MINOS
  NDOscCC_MC_minos = (TH1D*)CreateTotalSpectrum(	my_pars, 
							NDCC_TrueNC_minos, 
							NDCC_NuMu_minos, 
							NDCC_BeamNue_minos, 
							NDCC_AppNue_minos, 
							NDCC_AppNuTau_minos, 
							1.04*kKmUnits);
  FDOscCC_MC_minos = (TH1D*)CreateTotalSpectrum(	my_pars, 
							FDCC_TrueNC_minos, 
							FDCC_NuMu_minos, 
							FDCC_BeamNue_minos, 
							FDCC_AppNue_minos, 
							FDCC_AppNuTau_minos, 
							735.0*kKmUnits);

  //Oscillated CC MC -- MINOS+
  NDOscCC_MC_minosPlus = (TH1D*)CreateTotalSpectrum(	my_pars, 
							NDCC_TrueNC_minosPlus, 
							NDCC_NuMu_minosPlus, 
							NDCC_BeamNue_minosPlus, 
							NDCC_AppNue_minosPlus, 
							NDCC_AppNuTau_minosPlus, 
							1.04*kKmUnits);
  FDOscCC_MC_minosPlus = (TH1D*)CreateTotalSpectrum(	my_pars, 
							FDCC_TrueNC_minosPlus, 
							FDCC_NuMu_minosPlus, 
							FDCC_BeamNue_minosPlus, 
							FDCC_AppNue_minosPlus, 
							FDCC_AppNuTau_minosPlus, 
							735.0*kKmUnits);
}
//---------------------------------------------------------------------------------
TH1D* CalcChi2::GetTwoDetSpectrum(TH1D* hND, TH1D* hFD){

  int NDbins = hND->GetNbinsX();
  int FDbins = hFD->GetNbinsX();
  int Nbins = NDbins + FDbins;
  const int Nedges = Nbins + 1; 
 
  vector<double> edges(Nedges);

  edges[0]=0;
  
  double shift = 40.0;	//shift bin edges of ND spectrum by maximum energy 
			//of first spectrum to force increasing bin edges
  
  for(int i=1;i<=FDbins;i++){
    edges[i] = hFD->GetXaxis()->GetBinUpEdge(i);
  }
  for(int i=1;i<=NDbins;i++){
    edges[i+FDbins] = hND->GetXaxis()->GetBinUpEdge(i) + shift;
  }

  TH1D* hSpec = new TH1D("","",Nbins,&edges[0]);
  hSpec->Sumw2();

  for(int i=1;i<=FDbins;i++){
    hSpec->SetBinContent(i,hFD->GetBinContent(i));
    hSpec->SetBinError(i,hFD->GetBinError(i));
  }
  for(int i=1;i<=NDbins;i++){
    hSpec->SetBinContent(i+FDbins,hND->GetBinContent(i));
    hSpec->SetBinError(i+FDbins,hND->GetBinError(i));
  }

  return hSpec;
}
//---------------------------------------------------------------------------------
TMatrixD* CalcChi2::ScaleCovarianceMatrix(TH1D* pred, TMatrixD* mtx)
{
  Double_t binNumber = pred->GetNbinsX();
  TMatrixD* scaled_mtx = (TMatrixD*)mtx->Clone();
  TMatrixD* inverted_mtx = (TMatrixD*)mtx->Clone();

  Double_t bci, bcj, stat, syst, sig2;
  for(Int_t i=1; i<=binNumber; ++i){
    for(Int_t j=1; j<=binNumber; ++j){
      bci = pred->GetBinContent(i);
      bcj = pred->GetBinContent(j);
      //Poisson statistical uncertainty
      stat = 0;
      if(i==j){
	stat = bci;
      }
      //Systematic uncertainty
      syst = 0;
      syst = mtx->operator()(i-1,j-1);
      syst = bci*bcj*syst;
      //Sum stat and syst uncertainty
      sig2 = stat + syst;
      scaled_mtx->operator()(i-1,j-1) = sig2;
    }
  }
  
  TDecompSVD* DeCom = new TDecompSVD(*scaled_mtx);
  *inverted_mtx = DeCom->Invert();
  
  delete DeCom;
  delete scaled_mtx;
    
  return inverted_mtx;
}
//---------------------------------------------------------------------------------
Double_t CalcChi2::PenaltyTermDm232(Double_t dm232)
{
  Double_t dm232_pen = 0.0;
  dm232_pen = TMath::Power( (TMath::Abs(dm232) - 0.0025) , 2); // numerator
  dm232_pen /= TMath::Power( 0.0005, 2);

  return dm232_pen;
}
//---------------------------------------------------------------------------------
Double_t CalcChi2::PenaltyTermNuisance(Double_t par, Double_t mean, Double_t sigma)
{
  Double_t pen = 0.0;
  pen = TMath::Power( (par - mean) , 2); // numerator
  pen /= TMath::Power( sigma, 2);

  return pen;
}
//---------------------------------------------------------------------------------
Double_t CalcChi2::ChiSqFunction(TH1D* rPred, TH1D* rData, TMatrixD* CoVarInvert)
{
  if(!(rPred->GetNbinsX() == rData->GetNbinsX())){ 
    std::cout << "Binning Error. Asserting" << std::endl;
    assert(false);
  }

  Int_t NumberOfBins = rPred->GetNbinsX();

  TVectorD Difference(NumberOfBins);

  for(Int_t i=1; i<=NumberOfBins; ++i){
    Difference(i-1) = (rData->GetBinContent(i) - rPred->GetBinContent(i));
  }

  TVectorD temp = Difference;
  temp *= (*CoVarInvert);

  Double_t TotalChiSq = temp*Difference;

  return TotalChiSq;

}
//----------------------------------------------------------------------
Double_t CalcChi2::ComparePredWithData(TH1D* predCC, 
				       TH1D* dataCC, 
				       TMatrixD* CoVarCC_inverted, 
				       Double_t Dm2
				       )
{
  totalChi2_CC = ChiSqFunction(predCC, dataCC, CoVarCC_inverted);
  
  //Atmospheric mass-splitting penalty term
  Penalty_dm232   = PenaltyTermDm232(Dm2);

  totalChi2 = totalChi2_CC + Penalty_dm232;
  
  return totalChi2;
}
//---------------------------------------------------------------------------------
TH1D* CalcChi2::CreateTotalSpectrum(params my_pars,
				    TH2D* TrueNC,
				    TH2D* NuMu,
				    TH2D* BeamNue,
				    TH2D* AppNue,
				    TH2D* AppNuTau,
				    double baseline
				    )
{
  TH1D* vtruenc   = (TH1D*)CreateSpectrumComponent(my_pars, "TrueNC",   TrueNC,   baseline);
  TH1D* vnumu     = (TH1D*)CreateSpectrumComponent(my_pars, "NuMu",     NuMu,     baseline);
  TH1D* vbeamnue  = (TH1D*)CreateSpectrumComponent(my_pars, "BeamNue",  BeamNue,  baseline);
  TH1D* vappnue   = (TH1D*)CreateSpectrumComponent(my_pars, "AppNue",   AppNue,   baseline);
  TH1D* vappnutau = (TH1D*)CreateSpectrumComponent(my_pars, "AppNuTau", AppNuTau, baseline);

  TString name = vtruenc->GetName();
  name += "_total";
  TH1D* hTotal = (TH1D*)vtruenc->Clone(name);
  
  //TH1D* hTotal = new TH1D(*vtruenc);
  hTotal->Add(vnumu);
  hTotal->Add(vbeamnue);
  hTotal->Add(vappnue);
  hTotal->Add(vappnutau);

  vtruenc->Delete();
  vnumu->Delete();
  vbeamnue->Delete();
  vappnue->Delete();
  vappnutau->Delete();
  
  return hTotal;
}
//---------------------------------------------------------------------------------
double CalcChi2::FourFlavourNuMuToNuSProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline)
{ 

  // Calculate other mass splittings
  const double dm231 = dm221 + dm232;
  const double dm241 = dm231 + dm243;

  const double c12 = cos(theta12); const double s12 = sin(theta12);
  const double c13 = cos(theta13); const double s13 = sin(theta13);
  const double c14 = cos(theta14); const double s14 = sin(theta14);
  const double c23 = cos(theta23); const double s23 = sin(theta23);
  const double c24 = cos(theta24); const double s24 = sin(theta24);
  const double c34 = cos(theta34); const double s34 = sin(theta34);

  complex<double> expNegCP13 = complex<double>(cos(delta1), -sin(delta1));
  complex<double> expNegCP14 = complex<double>(cos(delta2), -sin(delta2));
  complex<double> expNegCP24 = complex<double>(cos(delta3), -sin(delta3));

  complex<double> Umu2  =  c12 * c23 * c24
                        -  c24 * s12 * s13 * s23 * conj(expNegCP13)
                        -  c13 * s12 * s14 * s24 * expNegCP24 * conj(expNegCP14);

  complex<double> Umu3  =  c13 * c24 * s23
                        -  s13 * s14 * s24 * expNegCP13 * expNegCP24 * conj(expNegCP14);
  
  complex<double> Umu4  =  c14 * s24 * expNegCP24;
  

  complex<double> Us2   =  -c13 * c24 * c34 * s12 * s14 * conj(expNegCP14)
                           -c12 * c23 * c34 * s24 * conj(expNegCP24)
                           +c34 * s12 * s13 * s23 * s24 * conj(expNegCP13 * expNegCP24)
                           +c23 * s12 * s13 * s34 * conj(expNegCP13)
                           +c12 * s23 * s34;
  
  complex<double> Us3   =  -c24 * c34 * s13 * s14 * expNegCP13 * conj(expNegCP14)
                           -c13 * c34 * s23 * s24 * conj(expNegCP24)
                           -c13 * c23 * s34;
  
  complex<double> Us4   =  c14 * c24 * c34;
  
  
  complex<double> i(0.0, 1.0);
  
  double DeltaM21 = (k1267 * dm221 * baseline / kKmUnits) / energy;
  double DeltaM31 = (k1267 * dm231 * baseline / kKmUnits) / energy;
  double DeltaM41 = (k1267 * dm241 * baseline / kKmUnits) / energy;
  
  complex<double> expDeltaM21 = complex<double>(cos(DeltaM21), -sin(DeltaM21));
  complex<double> expDeltaM31 = complex<double>(cos(DeltaM31), -sin(DeltaM31));
  complex<double> expDeltaM41 = complex<double>(cos(DeltaM41), -sin(DeltaM41));
  
  double oscProb = norm(-2.0 * i * conj(Umu2) * Us2 * sin(DeltaM21) * expDeltaM21
                        -2.0 * i * conj(Umu3) * Us3 * sin(DeltaM31) * expDeltaM31
			-2.0 * i * conj(Umu4) * Us4 * sin(DeltaM41) * expDeltaM41);

  
  return oscProb;
}
//---------------------------------------------------------------------------------
double CalcChi2::FourFlavourDisappearanceWeight
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline)
{

  // Calculate other mass splittings
  const double dm231 = dm221 + dm232;
  const double dm241 = dm231 + dm243;

  const double c12 = cos(theta12); const double s12 = sin(theta12);
  const double c13 = cos(theta13); const double s13 = sin(theta13);
  const double c14 = cos(theta14); const double s14 = sin(theta14);
  const double c23 = cos(theta23); const double s23 = sin(theta23);
  const double c24 = cos(theta24); const double s24 = sin(theta24);

  complex<double> expNegCP13 = complex<double>(cos(delta1), -sin(delta1));
  complex<double> expNegCP14 = complex<double>(cos(delta2), -sin(delta2));
  complex<double> expNegCP24 = complex<double>(cos(delta3), -sin(delta3));

  complex<double> Umu2  =  c12 * c23 * c24
                        -  c24 * s12 * s13 * s23 * conj(expNegCP13)
                        -  c13 * s12 * s14 * s24 * expNegCP24 * conj(expNegCP14);

  complex<double> Umu3  =  c13 * c24 * s23
                        -  s13 * s14 * s24 * expNegCP13 * expNegCP24 * conj(expNegCP14);
  
  complex<double> Umu4  =  c14 * s24 * expNegCP24;
  
  complex<double> i(0.0, 1.0);
  
  double DeltaM21 = (k1267 * dm221 * baseline / kKmUnits) / energy;
  double DeltaM31 = (k1267 * dm231 * baseline / kKmUnits) / energy;
  double DeltaM41 = (k1267 * dm241 * baseline / kKmUnits) / energy;
  
  complex<double> expDeltaM21 = complex<double>(cos(DeltaM21), -sin(DeltaM21));
  complex<double> expDeltaM31 = complex<double>(cos(DeltaM31), -sin(DeltaM31));
  complex<double> expDeltaM41 = complex<double>(cos(DeltaM41), -sin(DeltaM41));
  
  double oscProb  = norm(1.0 
			 - 2.0 * i * conj(Umu2) * Umu2 * sin(DeltaM21) * expDeltaM21 
			 - 2.0 * i * conj(Umu3) * Umu3 * sin(DeltaM31) * expDeltaM31 
			 - 2.0 * i * conj(Umu4) * Umu4 * sin(DeltaM41) * expDeltaM41);

  return oscProb;
}
//---------------------------------------------------------------------------------
double CalcChi2::FourFlavourNuESurvivalProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline)
{

  // Calculate other mass splittings
  const double dm231 = dm221 + dm232;
  const double dm241 = dm231 + dm243;

  const double s12 = sin(theta12);
  const double c13 = cos(theta13); const double s13 = sin(theta13);
  const double c14 = cos(theta14); const double s14 = sin(theta14);

  complex<double> expNegCP13 = complex<double>(cos(delta1), -sin(delta1));
  complex<double> expNegCP14 = complex<double>(cos(delta2), -sin(delta2));

  complex<double> Ue2   =  c13 * c14 * s12;
  complex<double> Ue3   =  c14 * s13 * expNegCP13;
  complex<double> Ue4   =  s14 * expNegCP14;
  
  complex<double> i(0.0, 1.0);
  
  double DeltaM21 = (k1267 * dm221 * baseline / kKmUnits) / energy;
  double DeltaM31 = (k1267 * dm231 * baseline / kKmUnits) / energy;
  double DeltaM41 = (k1267 * dm241 * baseline / kKmUnits) / energy;
  
  complex<double> expDeltaM21 = complex<double>(cos(DeltaM21), -sin(DeltaM21));
  complex<double> expDeltaM31 = complex<double>(cos(DeltaM31), -sin(DeltaM31));
  complex<double> expDeltaM41 = complex<double>(cos(DeltaM41), -sin(DeltaM41));
  
  double oscProb = norm(1.0 
			- 2.0 * i * conj(Ue2) * Ue2 * sin(DeltaM21) * expDeltaM21 
			- 2.0 * i * conj(Ue3) * Ue3 * sin(DeltaM31) * expDeltaM31 
			- 2.0 * i * conj(Ue4) * Ue4 * sin(DeltaM41) * expDeltaM41);

  return oscProb;
}
//---------------------------------------------------------------------------------
double CalcChi2::FourFlavourNuMuToNuEProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline)
{

  // Calculate other mass splittings
  const double dm231 = dm221 + dm232;
  const double dm241 = dm231 + dm243;

  const double c12 = cos(theta12); const double s12 = sin(theta12);
  const double c13 = cos(theta13); const double s13 = sin(theta13);
  const double c14 = cos(theta14); const double s14 = sin(theta14);
  const double c23 = cos(theta23); const double s23 = sin(theta23);
  const double c24 = cos(theta24); const double s24 = sin(theta24);

  complex<double> expNegCP13 = complex<double>(cos(delta1), -sin(delta1));
  complex<double> expNegCP14 = complex<double>(cos(delta2), -sin(delta2));
  complex<double> expNegCP24 = complex<double>(cos(delta3), -sin(delta3));
  
  complex<double> Umu2  =  c12 * c23 * c24
                        -  c24 * s12 * s13 * s23 * conj(expNegCP13)
                        -  c13 * s12 * s14 * s24 * expNegCP24 * conj(expNegCP14);

  complex<double> Umu3  =  c13 * c24 * s23
                        -  s13 * s14 * s24 * expNegCP13 * expNegCP24 * conj(expNegCP14);
  
  complex<double> Umu4  =  c14 * s24 * expNegCP24;


  complex<double> Ue2   =  c13 * c14 * s12;
  complex<double> Ue3   =  c14 * s13 * expNegCP13;
  complex<double> Ue4   =  s14 * expNegCP14;
  
  complex<double> i(0.0, 1.0);
  
  double DeltaM21 = (k1267 * dm221 * baseline / kKmUnits) / energy;
  double DeltaM31 = (k1267 * dm231 * baseline / kKmUnits) / energy;
  double DeltaM41 = (k1267 * dm241 * baseline / kKmUnits) / energy;
  
  complex<double> expDeltaM21 = complex<double>(cos(DeltaM21), -sin(DeltaM21));
  complex<double> expDeltaM31 = complex<double>(cos(DeltaM31), -sin(DeltaM31));
  complex<double> expDeltaM41 = complex<double>(cos(DeltaM41), -sin(DeltaM41));
  
  double oscProb = norm(-2.0 * i * conj(Umu2) * Ue2 * sin(DeltaM21) * expDeltaM21       
			-2.0 * i * conj(Umu3) * Ue3 * sin(DeltaM31) * expDeltaM31
			-2.0 * i * conj(Umu4) * Ue4 * sin(DeltaM41) * expDeltaM41);

  return oscProb;
}
//---------------------------------------------------------------------------------
double CalcChi2::FourFlavourNuMuToNuTauProbability
(const double energy, 
 double dm232, const double theta23, double dm221, 
 double dm243, const double theta12, 
 const double theta13, const double theta14,
 const double theta24, const double theta34,
 double delta1, double delta2, double delta3,
 const double baseline)
{

  // Calculate other mass splittings
  const double dm231 = dm221 + dm232;
  const double dm241 = dm231 + dm243;

  const double c12 = cos(theta12); const double s12 = sin(theta12);
  const double c13 = cos(theta13); const double s13 = sin(theta13);
  const double c14 = cos(theta14); const double s14 = sin(theta14);
  const double c23 = cos(theta23); const double s23 = sin(theta23);
  const double c24 = cos(theta24); const double s24 = sin(theta24);
  const double c34 = cos(theta34); const double s34 = sin(theta34);

  complex<double> expNegCP13 = complex<double>(cos(delta1), -sin(delta1));
  complex<double> expNegCP14 = complex<double>(cos(delta2), -sin(delta2));
  complex<double> expNegCP24 = complex<double>(cos(delta3), -sin(delta3));

  complex<double> Umu2  =  c12 * c23 * c24
                        -  c24 * s12 * s13 * s23 * conj(expNegCP13)
                        -  c13 * s12 * s14 * s24 * expNegCP24 * conj(expNegCP14);

  complex<double> Umu3  =  c13 * c24 * s23
                        -  s13 * s14 * s24 * expNegCP13 * expNegCP24 * conj(expNegCP14);
  
  complex<double> Umu4  =  c14 * s24 * expNegCP24;


  complex<double> Utau2 =  -c12 * c34 * s23
                           -c23 * c34 * s12 * s13 * conj(expNegCP13)
                           -c13 * c24 * s12 * s14 * s34 * conj(expNegCP14)
                           -c12 * c23 * s24 * s34 * conj(expNegCP24)
                           +s12 * s13 * s23 * s24 * s34 * conj(expNegCP13 * expNegCP24);

  complex<double> Utau3 =  c13 * c23 * c34
                        -  c24 * s13 * s14 * s34 * expNegCP13 * conj(expNegCP14)
                        -  c13 * s23 * s24 * s34 * conj(expNegCP24);

  complex<double> Utau4 =  c14 * c24 * s34;

  complex<double> i(0.0, 1.0);
  
  double DeltaM21 = (k1267 * dm221 * baseline / kKmUnits) / energy;
  double DeltaM31 = (k1267 * dm231 * baseline / kKmUnits) / energy;
  double DeltaM41 = (k1267 * dm241 * baseline / kKmUnits) / energy;
  
  complex<double> expDeltaM21 = complex<double>(cos(DeltaM21), -sin(DeltaM21));
  complex<double> expDeltaM31 = complex<double>(cos(DeltaM31), -sin(DeltaM31));
  complex<double> expDeltaM41 = complex<double>(cos(DeltaM41), -sin(DeltaM41));
  
  double oscProb =  norm(-2.0 * i * conj(Umu2) * Utau2 * sin(DeltaM21) * expDeltaM21     
			 -2.0 * i * conj(Umu3) * Utau3 * sin(DeltaM31) * expDeltaM31
			 -2.0 * i * conj(Umu4) * Utau4 * sin(DeltaM41) * expDeltaM41);

  return oscProb;
}
//---------------------------------------------------------------------------------
TH1D* CalcChi2::CreateSpectrumComponent(params my_pars, TString OscType, TH2D* oscDummy, Double_t baseline)
{
  TString name = oscDummy->GetName();
  name += "_";
  name += OscType;
  name += "_";
  name += baseline;
  name += "_py";
  
  TH1D* bintemplate = oscDummy->ProjectionY(name);
  bintemplate->Reset();
 
  const double k1267 = 1.26693276;

  // Loop over every true energy bin in the reco vs. true matrices, then loop over every reco energy in that bin                                                 
  // to calculate an oscillation weight for that reco energy based on the true energy. 
  TAxis *Yaxis = oscDummy->GetYaxis();
  TAxis *Xaxis = oscDummy->GetXaxis();

  // Define Dm243 such that its actually Dm241 being altered.
  //41 = 43 + 32 + 21 
  //43 = 41 - 32 - 21
  Double_t dm243 = 0.0;

  dm243 = my_pars.Dm241 - my_pars.Dm232 - my_pars.Dm221;

  for(Int_t x = 1; x <= Xaxis->GetNbins(); x++){
    Double_t OscWeight = 0.0;
    
    if(baseline > 0){
      
      // Default iterations (1 at bin center)
      Int_t n_LoverE = 1;
      Double_t LoverE[5];
      LoverE[0] = Xaxis->GetBinCenter(x);
      
      // This is averaging oscialltions in true energy bins - see Technical Note http://minos-docdb.fnal.gov/cgi-bin/RetrieveFile?docid=10203&version=2
      const Double_t W = Xaxis->GetBinWidth(x);
      const Double_t arg = k1267*dm243*W; // half-period of oscillation
      Double_t sample = W/2/sqrt(3);

      if(arg!=0) sample = TMath::ACos(TMath::Sin(arg)/arg)/arg*W/2;

      n_LoverE = 2;
      Double_t bc = LoverE[0]; // bin center
      LoverE[0] = bc - sample;
      LoverE[1] = bc + sample;

      Double_t E = 1.0;
      Double_t L;
      
      for(int i = 0; i < n_LoverE; i++){
	L = LoverE[i];
	E = 1.0;
	
	// each Osctype has a different probability function
	if(OscType == "TrueNC"){
	  
	  OscWeight += FourFlavourNuMuToNuSProbability( E, 
  							my_pars.Dm232, 
							my_pars.th23, 
							my_pars.Dm221, 
							dm243, 
							my_pars.th12, 
							my_pars.th13, 
							my_pars.th14, 
							my_pars.th24, 
							my_pars.th34, 
							my_pars.deltaCP,  
							my_pars.delta14, 
							my_pars.delta24,
							L*kKmUnits);
	}
	if(OscType == "NuMu"){

	  OscWeight += FourFlavourDisappearanceWeight( E, 
  							my_pars.Dm232, 
							my_pars.th23, 
							my_pars.Dm221, 
							dm243, 
							my_pars.th12, 
							my_pars.th13, 
							my_pars.th14, 
							my_pars.th24, 
							my_pars.th34, 
							my_pars.deltaCP,  
							my_pars.delta14, 
							my_pars.delta24,
							L*kKmUnits);
	}
	if(OscType == "BeamNue"){
	  
	  OscWeight += FourFlavourNuESurvivalProbability( E, 
  							my_pars.Dm232, 
							my_pars.th23, 
							my_pars.Dm221, 
							dm243, 
							my_pars.th12, 
							my_pars.th13, 
							my_pars.th14, 
							my_pars.th24, 
							my_pars.th34, 
							my_pars.deltaCP,  
							my_pars.delta14, 
							my_pars.delta24,
							L*kKmUnits);
	}
	if(OscType == "AppNue"){
	  
	  OscWeight += FourFlavourNuMuToNuEProbability( E, 
  							my_pars.Dm232, 
							my_pars.th23, 
							my_pars.Dm221, 
							dm243, 
							my_pars.th12, 
							my_pars.th13, 
							my_pars.th14, 
							my_pars.th24, 
							my_pars.th34, 
							my_pars.deltaCP,  
							my_pars.delta14, 
							my_pars.delta24,
							L*kKmUnits);
	}
	if(OscType == "AppNuTau"){

	  OscWeight += FourFlavourNuMuToNuTauProbability( E, 
  							my_pars.Dm232, 
							my_pars.th23, 
							my_pars.Dm221, 
							dm243, 
							my_pars.th12, 
							my_pars.th13, 
							my_pars.th14, 
							my_pars.th24, 
							my_pars.th34, 
							my_pars.deltaCP,  
							my_pars.delta14, 
							my_pars.delta24,
							L*kKmUnits);
	}
      }
      // Now average this
      OscWeight /= n_LoverE;
    }
    else // if baseline < 0
      {
	
        if(OscType == "TrueNC")   OscWeight = 0.0;
        if(OscType == "NuMu")     OscWeight = 1.0;
        if(OscType == "BeamNue")  OscWeight = 1.0;
        if(OscType == "AppNue")   OscWeight = 0.0;
        if(OscType == "AppNuTau") OscWeight = 0.0;
      }

    // using the oscillation weight, fill a 1d histogram for each type of event with the oscillated reco energy 
    for(Int_t y = 1; y <= Yaxis->GetNbins(); y++){
      
      Double_t sumWeights = 0;
      
      if(OscType == "TrueNC"){
	sumWeights += oscDummy->GetBinContent(x,y)*(1.0-OscWeight);
      }
      else{
	sumWeights += oscDummy->GetBinContent(x,y)*(OscWeight);
      }
      Double_t currBinContents = bintemplate->GetBinContent( y );
      bintemplate->SetBinContent( y, sumWeights + currBinContents);

    }

  }

  return bintemplate;
}
