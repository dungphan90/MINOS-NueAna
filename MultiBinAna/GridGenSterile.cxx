#define GridGen_C

#include "NueAna/MultiBinAna/GridGenSterile.h"
#include <cstdlib>

GridGenSterile::GridGenSterile() {
  SetOutputFile();
  SetNExperiments();

  SetNormalHierarchy(true);
  SetTheta12();
  SetTheta13();
  SetTheta23();
  SetAbsValDeltaMSq23();
  SetDeltaMSq12();

  Extrap.clear();
  return;
}

GridGenSterile::~GridGenSterile() {}

void GridGenSterile::AddExtrap(Extrapolate2D *E) {
  Extrap.push_back(E);
  return;
}

void GridGenSterile::SetTheta12(double val, double errup, double errdn) {
  Theta12 = val;
  dTheta12_up = errup;
  dTheta12_dn = errdn;
  return;
}

void GridGenSterile::SetTheta23(double val, double errup, double errdn) {
  Theta23 = val;
  dTheta23_up = errup;
  dTheta23_dn = errdn;
  return;
}

void GridGenSterile::SetAbsValDeltaMSq23(double val, double errup, double errdn) {
  DeltaMSq23 = TMath::Abs(val);
  dDeltaMSq23_up = errup;
  dDeltaMSq23_dn = errdn;
  return;
}

void GridGenSterile::SetDeltaMSq12(double val, double errup, double errdn) {
  DeltaMSq12 = val;
  dDeltaMSq12_up = errup;
  dDeltaMSq12_dn = errdn;
  return;
}

void GridGenSterile::SetTheta13(double val, double errup, double errdn) {
  Theta13 = val;
  dTheta13_up = errup;
  dTheta13_dn = errdn;
  return;
}

double GridGenSterile::AsymGaus(double sigma_minus, double sigma_plus) {
  double u = gRandom->Uniform();
  double g = gRandom->Gaus(0, 1);
  double p = (u > sigma_plus / (sigma_plus + sigma_minus)) ? -TMath::Abs(g * sigma_minus) : TMath::Abs(g * sigma_plus);

  return p;
}

void GridGenSterile::RunMultiBinOscParErrsSterileFit(string oscpardist_name) {
  if (Extrap.size() == 0) {
    std::cout << "No Extrapolate2D input.  Quitting..." << std::endl;
    return;
  }

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap.at(ie)->GetPrediction();
  }

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    Extrap.at(ie)->SetOscPar(OscPar::kTh12, Theta12);
    Extrap.at(ie)->SetOscPar(OscPar::kTh13, Theta13);
    Extrap.at(ie)->SetOscPar(OscPar::kTh23, Theta23);
    Extrap.at(ie)->SetOscPar(OscPar::kDeltaM12, DeltaMSq12);
    Extrap.at(ie)->SetOscPar(OscPar::kDeltaM23, DeltaMSq23);
    Extrap.at(ie)->SetOscPar(OscPar::kDelta, 0);
    // Sterile params
    Extrap.at(ie)->SetOscPar(OscPar::kTh14, 0);
    Extrap.at(ie)->SetOscPar(OscPar::kTh24, 0);
    Extrap.at(ie)->SetOscPar(OscPar::kTh34, 0);
    Extrap.at(ie)->SetOscPar(OscPar::kDm41, 0);
    Extrap.at(ie)->SetOscPar(OscPar::kDelta14, 0);
    Extrap.at(ie)->SetOscPar(OscPar::kDelta24, 0);
    if (!NormalHier) {
      Extrap.at(ie)->InvertMassHierarchy();
    }
    Extrap.at(ie)->OscillatePrediction();
  }

  int nbins = Extrap[0]->Pred_TotalBkgd_VsBinNumber->GetNbinsX();
  int nPID = Extrap[0]->GetNPID();
  std::vector<double> sig, bkgd;
  std::vector<std::vector<double> > nc, numucc, bnue, tau, nue;
  std::vector<double> oscparerr;
  std::vector<double> oscparerr_offdiag;

  for (int i = 0; i < nbins; i++) {
    sig.push_back(0);
    bkgd.push_back(0);
    oscparerr.push_back(0);
    nc.push_back(std::vector<double>());
    numucc.push_back(std::vector<double>());
    bnue.push_back(std::vector<double>());
    tau.push_back(std::vector<double>());
    nue.push_back(std::vector<double>());
    for (int k = i + 1; k < nbins; k++) {
      oscparerr_offdiag.push_back(0);
    }
    for (unsigned int j = 0; j < Extrap.size(); j++) {
      nc.at(i).push_back(0);
      numucc.at(i).push_back(0);
      bnue.at(i).push_back(0);
      tau.at(i).push_back(0);
      nue.at(i).push_back(0);
    }
  }

  double grid_dmsq41;
  double grid_ssqth14;
  double grid_ssqth24;
  double profiling_delta24;
  double profiling_delta14;
  double profiling_delta13;
  double profiling_th34;
  std::vector<TTree *> ftree;
  std::vector<std::vector<TTree *> > ftree2;
  int noff = 0;
  for (int i = 0; i < nbins; i++) {
    TTree* ttmp = new TTree(Form("Bin_%i", i), Form("Bin_%i", i));
    ttmp->Branch("Dmsq41",                      &grid_dmsq41,             "Dmsq41/D");
    ttmp->Branch("SinSqTh14",                   &grid_ssqth14,            "SinSqTh14/D");
    ttmp->Branch("SinSqTh24",                   &grid_ssqth24,            "SinSqTh24/D");
    ttmp->Branch("Delta13",                     &profiling_delta13,       "Delta13/D");
    ttmp->Branch("Delta14",                     &profiling_delta14,       "Delta14/D");
    ttmp->Branch("Delta24",                     &profiling_delta24,       "Delta24/D");
    ttmp->Branch("Theta34",                     &profiling_th34,          "Theta34/D");
    ttmp->Branch("Signal",                      &sig.at(i),               "Signal/D");
    ttmp->Branch("Background",                  &bkgd.at(i),              "Background/D");
    ttmp->Branch("DNExp_DOscPars",              &oscparerr.at(i),         "DNExp_DOscPars/D");
    for (int k = i + 1; k < nbins; k++) {
      ttmp->Branch(Form("Bin_%i_Bin_%i", i, k), &oscparerr_offdiag.at(noff), Form("Bin_%i_Bin_%i/D", i, k));
      noff++;
    }
    ftree.push_back(ttmp);

    ftree2.push_back(std::vector<TTree*>());
    for (unsigned int j = 0; j < Extrap.size(); j++) {
      TTree* ttmp2 = new TTree(Form("Bin_%i_Run_%i", i, j), Form("Bin_%i_Run_%i", i, j));
      ttmp2->Branch("Dmsq41",          &grid_dmsq41,            "Dmsq41/D");
      ttmp2->Branch("SinSqTh14",       &grid_ssqth14,           "SinSqTh14/D");
      ttmp2->Branch("SinSqTh24",       &grid_ssqth24,           "SinSqTh24/D");
      ttmp2->Branch("Delta13",         &profiling_delta13,      "Delta13/D");
      ttmp2->Branch("Delta14",         &profiling_delta14,      "Delta14/D");
      ttmp2->Branch("Delta24",         &profiling_delta24,      "Delta14/D");
      ttmp2->Branch("Theta34",         &profiling_th34,         "Theta34/D");
      ttmp2->Branch("Signal",          &nue.at(i).at(j),        "Signal/D");
      ttmp2->Branch("NC",              &nc.at(i).at(j),         "NC/D");
      ttmp2->Branch("NuMuCC",          &numucc.at(i).at(j),     "NuMuCC/D");
      ttmp2->Branch("BNueCC",          &bnue.at(i).at(j),       "BNueCC/D");
      ttmp2->Branch("NuTauCC",         &tau.at(i).at(j),        "NuTauCC/D");
      ftree2.at(i).push_back(ttmp2);
    }
  }

  std::vector<double> nexp, nobs;
  std::vector<TH1D *> delnexphist;
  std::vector<TH1D *> delnexphist_offdiag;
  for (int i = 0; i < nbins; i++) {
    delnexphist.push_back(new TH1D(Form("delnexphist_%i", i), "", 400, -1, 1));
    nexp.push_back(0);
    nobs.push_back(0);
    for (int k = i + 1; k < nbins; k++) {
      delnexphist_offdiag.push_back(new TH1D(Form("delnexphist_%i_%i", i, k), "", 400, -1, 1));
    }
  }

  gRandom->SetSeed(time(NULL));

  TFile *oscpardist_file = new TFile(gSystem->ExpandPathName(oscpardist_name.c_str()), "RECREATE");

  grid_dmsq41  = DMSQ41Grid;    // THOSE VALUES ARE FIXED ON THE GRID
  grid_ssqth14 = SinSqTh14Grid; // THOSE VALUES ARE FIXED ON THE GRID
  grid_ssqth24 = SinSqTh24Grid; // THOSE VALUES ARE FIXED ON THE GRID
  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    // Set these fixed parameters and not to change them in the rest of the function
    Extrap.at(ie)->SetOscPar(OscPar::kDm41, grid_dmsq41);
    Extrap.at(ie)->SetOscPar(OscPar::kTh14, TMath::ASin(TMath::Sqrt(grid_ssqth14)));
    Extrap.at(ie)->SetOscPar(OscPar::kTh24, TMath::ASin(TMath::Sqrt(grid_ssqth24)));
  }

  profiling_delta14 = (double)(((double)(rand() % 100)) / 100) * 2.0 * TMath::Pi(); // Randomize profiling parameters
  profiling_delta24 = (double)(((double)(rand() % 100)) / 100) * 2.0 * TMath::Pi(); // Randomize profiling parameters
  profiling_delta13 = (double)(((double)(rand() % 100)) / 100) * 2.0 * TMath::Pi(); // Randomize profiling parameters
  profiling_th34    = (double)(((double)(rand() % 100)) / 100) * 0.5 * TMath::Pi(); // Randomize profiling parameters

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    // Get nominal prediction
    Extrap.at(ie)->SetOscPar(OscPar::kTh34,     profiling_th34);
    Extrap.at(ie)->SetOscPar(OscPar::kDelta,    profiling_delta13);
    Extrap.at(ie)->SetOscPar(OscPar::kDelta14,  profiling_delta14);
    Extrap.at(ie)->SetOscPar(OscPar::kDelta24,  profiling_delta24);
    Extrap.at(ie)->SetOscPar(OscPar::kTh12,     Theta12);
    Extrap.at(ie)->SetOscPar(OscPar::kTh13,     Theta13);
    Extrap.at(ie)->SetOscPar(OscPar::kTh23,     Theta23);
    Extrap.at(ie)->SetOscPar(OscPar::kDeltaM12, DeltaMSq12);
    Extrap.at(ie)->SetOscPar(OscPar::kDeltaM23, DeltaMSq23);
    if (!NormalHier) { Extrap.at(ie)->InvertMassHierarchy(); }
    Extrap.at(ie)->OscillatePrediction();
  }

  // Here calculate the NExp due to the nominal prediction
  for (int i = 0; i < nbins; i++) {
    sig[i] = 0;
    bkgd[i] = 0;
    int ir = (int) (i / nPID);
    int ip = i % nPID;
    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      bkgd   .at(i)        += Extrap.at(ie)->Pred_TotalBkgd_VsBinNumber ->GetBinContent(i + 1);
      sig    .at(i)        += Extrap.at(ie)->Pred_Signal_VsBinNumber    ->GetBinContent(i + 1);
      nc     .at(i).at(ie)  = Extrap.at(ie)->Pred[Background::kNC]      ->GetBinContent(ip + 1, ir + 1);
      numucc .at(i).at(ie)  = Extrap.at(ie)->Pred[Background::kNuMuCC]  ->GetBinContent(ip + 1, ir + 1);
      bnue   .at(i).at(ie)  = Extrap.at(ie)->Pred[Background::kBNueCC]  ->GetBinContent(ip + 1, ir + 1);
      tau    .at(i).at(ie)  = Extrap.at(ie)->Pred[Background::kNuTauCC] ->GetBinContent(ip + 1, ir + 1);
      nue    .at(i).at(ie)  = Extrap.at(ie)->Pred[Background::kNueCC]   ->GetBinContent(ip + 1, ir + 1);
    }
    nexp[i] = sig[i] + bkgd[i];
  }

  // Generate pseudo experiments
  noff = 0;
  for (int i = 0; i < nbins; i++) {
    delnexphist.at(i)->Reset();
    delnexphist.at(i)->SetName(Form("DeltaNexp_Diag_%i", i));
    for (int k = i + 1; k < nbins; k++) {
      delnexphist_offdiag.at(noff)->Reset();
      delnexphist_offdiag.at(noff)->SetName(Form("DeltaNexp_OffDiag_%i_%i", i, k));
      noff++;
    }
  }

  for (int u = 0; u < NumExpts; u++) {
    // Small fluctuation of 3-flavor parameters
    double fluctuating_theta12 = Theta12    + AsymGaus(dTheta12_dn,    dTheta12_up);
    double fluctuating_theta13 = Theta13    + AsymGaus(dTheta13_dn,    dTheta13_up);
    double fluctuating_theta23 = Theta23    + AsymGaus(dTheta23_dn,    dTheta23_up);
    double fluctuating_dm21    = DeltaMSq12 + AsymGaus(dDeltaMSq12_dn, dDeltaMSq12_up);
    double fluctuating_dm32    = DeltaMSq23 + AsymGaus(dDeltaMSq23_dn, dDeltaMSq23_up);

    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      Extrap.at(ie)->SetOscPar(OscPar::kTh12,     fluctuating_theta12);
      Extrap.at(ie)->SetOscPar(OscPar::kTh13,     fluctuating_theta13);
      Extrap.at(ie)->SetOscPar(OscPar::kTh23,     fluctuating_theta23);
      Extrap.at(ie)->SetOscPar(OscPar::kDeltaM12, fluctuating_dm21);
      Extrap.at(ie)->SetOscPar(OscPar::kDeltaM23, fluctuating_dm32);
      if (!NormalHier) { Extrap.at(ie)->InvertMassHierarchy(); }
      Extrap.at(ie)->OscillatePrediction();
    }

    noff = 0;
    // Generate an NObs and "compare" this NObs to the NExp (calculated above, outside of the pseudo-experiment loop)
    for (int i = 0; i < nbins; i++) {
        nobs.at(i) = 0;
        for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
          nobs.at(i) += (Extrap.at(ie)->Pred_TotalBkgd_VsBinNumber->GetBinContent(i + 1));
          nobs.at(i) += (Extrap.at(ie)->Pred_Signal_VsBinNumber->GetBinContent(i + 1));
        }
        delnexphist.at(i)->Fill((nobs.at(i) - nexp.at(i)) / (nexp.at(i)));
    }
    // Bin-to-bin shift (due to the fluctuation of the 3-flavor parameters)
    for (int i = 0; i < nbins; i++) {
        for (int k = i + 1; k < nbins; k++) {
          delnexphist_offdiag.at(noff)->Fill((nobs.at(i) - nexp.at(i)) * (nobs.at(k) - nexp.at(k)) / (nexp.at(i) * nexp.at(k)));
          noff++;
        }
    }
  } // End of pseudo-experiment loop

  noff = 0;
  // Save the compare histogram
  for (int i = 0; i < nbins; i++) {
    oscparerr.at(i) = delnexphist.at(i)->GetRMS();
    delnexphist.at(i)->Write();
    for (int k = i + 1; k < nbins; k++) {
      oscparerr_offdiag.at(noff) = delnexphist_offdiag.at(noff)->GetRMS();
      if (delnexphist_offdiag.at(noff)->GetMean() < 0) oscparerr_offdiag.at(noff) = -1. * oscparerr_offdiag.at(noff);
      delnexphist_offdiag.at(noff)->Write();
      noff++;
    }
  }

  // Close the oscpardist file
  oscpardist_file->Close();

  gROOT->cd("/");

  for (int i = 0; i < nbins; i++) {
    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      ftree2.at(i).at(ie)->Fill();
    }
    ftree.at(i)->Fill();
  }

  double nPOTNear, nPOTFar;

  TTree *paramtree = new TTree("paramtree", "paramtree");
  paramtree->Branch("nearPOT",     &nPOTNear,        "nearPOT/D");
  paramtree->Branch("farPOT",      &nPOTFar,         "farPOT/D");
  paramtree->Branch("Theta12",     &Theta12,         "Theta12/D");
  paramtree->Branch("Theta13",     &Theta13,         "Theta13/D");
  paramtree->Branch("Theta23",     &Theta23,         "Theta23/D");
  paramtree->Branch("DeltaMSq23",  &DeltaMSq23,      "DeltaMSq23/D");
  paramtree->Branch("DeltaMSq12",  &DeltaMSq12,      "DeltaMSq12/D");
  paramtree->Branch("SinSqTh14",   &grid_ssqth14,    "SinSqTh14/D");
  paramtree->Branch("SinSqTh24",   &grid_ssqth24,    "SinSqTh24/D");
  paramtree->Branch("DeltaMSq41",  &grid_dmsq41,     "DeltaMSq41/D");

  if (!NormalHier) {
    DeltaMSq23 = -1. * DeltaMSq23;
  }

  for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
    nPOTNear = Extrap.at(ie)->GetNearPOT();
    nPOTFar = Extrap.at(ie)->GetFarPOT();
    paramtree->Fill();
  }

  TFile *fout = new TFile(gSystem->ExpandPathName(outFileName.c_str()), "RECREATE");
  for (int i = 0; i < nbins; i++) {
    ftree.at(i)->Write();
    for (unsigned int ie = 0; ie < Extrap.size(); ie++) {
      ftree2.at(i).at(ie)->Write();
    }
  }
  paramtree->Write();
  fout->Close();

  return;
}